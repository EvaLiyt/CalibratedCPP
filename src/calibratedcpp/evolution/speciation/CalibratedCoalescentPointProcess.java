package calibratedcpp.evolution.speciation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

import calibratedcpp.BirthDeathModel;
import calibratedcpp.CoalescentPointProcessModel;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Collections;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Marcus Overwater
 */

@Description("A general class of birth-death processes with incomplete extant sampling and conditioning on clade calibrations")
public class CalibratedCoalescentPointProcess extends SpeciesTreeDistribution {
    public Input<CoalescentPointProcessModel> cppModelInput =
            new Input<>("treeModel", "The tree model", (CoalescentPointProcessModel) null);

    public Input<List<CalibrationPoint>> calibrationsInput =
            new Input<>("calibrations","Clade calibrations", (List<CalibrationPoint>) null);

    public Input<Boolean> conditionOnCalibrationInput =
            new Input<>("conditionOnCalibrations","Boolean if the likelihood is conditioned on the clade calibrations (Default: true). For large trees with many calibrations it is recommended to set this to false and use the exchange operator.", true);

    protected TreeInterface tree;

    protected CoalescentPointProcessModel model;
    protected List<CalibrationPoint> calibrations = new ArrayList<>();
    protected boolean conditionOnCalibrations;

    protected Double origin;
    protected Double rootAge;

    protected boolean conditionOnRoot;
    protected double maxTime;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();
        model = cppModelInput.get();
        calibrations = calibrationsInput.get();
        conditionOnCalibrations = (calibrations != null) ? conditionOnCalibrationInput.get() : false;

        if (conditionOnCalibrations) {
            calibrations = postOrderTopologicalSort(tree, calibrations);
        }

        origin = model.getOrigin();

        conditionOnRoot = model.isConditionOnRoot();

        rootAge = tree.getRoot().getHeight();

        if (!conditionOnRoot && origin < rootAge) {
            throw new RuntimeException("rootAge (" + rootAge + ") > origin (" + origin + "): Origin must be greater than root age.");
        }

        maxTime = conditionOnRoot ? rootAge : origin;
    }

    public double calculateUnConditionedTreeLogLikelihood(TreeInterface tree) {

        double logP = Math.log1p(-Math.exp(model.calculateLogCDF(maxTime)));

        if (conditionOnRoot) {
            logP += logP - model.calculateLogDensity(maxTime);
        }

        for (Node node : tree.getInternalNodes()){
            double age = node.getHeight();
            logP += model.calculateLogDensity(age);
        }
        return logP;
    }

    public double calculateLogMarginalDensityOfCalibrations(TreeInterface tree, List<CalibrationPoint> calibrations) {
        double marginalDensity = Math.log1p(-Math.exp(model.calculateLogCDF(maxTime)));

        if (conditionOnRoot) {
            marginalDensity += marginalDensity;
        }

        int numTaxa = tree.getLeafNodeCount();
        int numCalibrations = calibrations.size();

        double[] calibrationAges = new double[numCalibrations];
        int[] cladeSizes = new int[numCalibrations];
        int sumCladeSizes = 0;

        double logQt = model.calculateLogCDF(maxTime);

        double[] logQti = new double[numCalibrations];
        double[] logDiff = new double[numCalibrations];

        for (int i = 0 ; i < numCalibrations ; i++) {
            CalibrationPoint calibration = calibrations.get(i);

            calibrationAges[i] = getMRCA(tree, calibration.taxa().asStringList()).getHeight();
            cladeSizes[i] = calibration.taxa().getTaxonSet().size();
            sumCladeSizes += cladeSizes[i];

            logQti[i] = model.calculateLogCDF(calibrationAges[i]);
            logDiff[i] = logDiffExp(logQt, logQti[i]); // Precompute log difference of Q(t) - Q(t_i)

            marginalDensity += model.calculateLogDensity(calibrationAges[i]) + (cladeSizes[i] - 2) * logQti[i];
        }

        marginalDensity += calculateLogSumOfPermutations(numCalibrations, sumCladeSizes, numTaxa, logQt, logDiff);

        return marginalDensity;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        logP = calculateUnConditionedTreeLogLikelihood(tree);

        if (conditionOnCalibrations) {
            logP -= calculateLogMarginalDensityOfCalibrations(tree, calibrations);
        }

        return logP;
    }

    private double calculateLogSumOfPermutations(int numCalibrations, int sumCladeSizes, int numTaxa, double logQ_t, double[] logDiff) {

        int m = numTaxa - sumCladeSizes;

        int totalElements = m + numCalibrations;

        List<List<Integer>> permutations = new ArrayList<>();
        generatePermutations(totalElements, numCalibrations, new ArrayList<>(), new boolean[totalElements], permutations);

        List<Double> logTerms = new ArrayList<>();

        double logTerm = 0;

        for (List<Integer> perm : permutations) {
            int sum_s = 0;
            int[] s = new int[numCalibrations];

            for (int i = 0; i < numCalibrations; i++) {
                int ell_i = perm.get(i);
                int countAdj = 0;
                for (int j = 0; j < i; j++) {
                    int ell_j = perm.get(j);
                    if (ell_j == ell_i - 1 || ell_j == ell_i + 1) {
                        countAdj++;
                    }
                }
                s[i] = (ell_i == 0 ? 1 : 0) + (ell_i == totalElements - 1 ? 1 : 0) + countAdj;
                sum_s += s[i];

                logTerm += (2 - s[i]) * logDiff[i]; // Compute log term for this permutation
            }

            logTerm += (m - numCalibrations - 1 + sum_s) * logQ_t;

            logTerms.add(logTerm);
        }
        // Use logSumExp to sum all terms in log space
        return logSumExp(logTerms);
    }

    private Node getMRCA(TreeInterface tree, List<String> taxonIDs) {
        // Get leaf nodes for the given taxa
        List<Node> nodes = new ArrayList<>();
        for (String id : taxonIDs) {
            int index = tree.getTaxonset().getTaxonIndex(id);
            nodes.add(tree.getNode(index));
        }

        return findMRCA(nodes);
    }

    private Node findMRCA(List<Node> nodes) {
        if (nodes == null || nodes.isEmpty()) return null;
        if (nodes.size() == 1) return nodes.get(0);

        Node mrca = nodes.get(0);

        for (int i = 1; i < nodes.size(); i++) {
            mrca = findMRCA(mrca, nodes.get(i));
        }

        return mrca;
    }

    private Node findMRCA(Node a, Node b) {
        // Collect ancestors of a
        List<Node> ancestorsA = new ArrayList<>();
        while (a != null) {
            ancestorsA.add(a);
            a = a.getParent();
        }

        // Walk b's ancestry and return first shared node
        while (b != null) {
            for (Node ancestor : ancestorsA) {
                if (b == ancestor) { // use reference equality!
                    return b;
                }
            }
            b = b.getParent();
        }

        return null; // Should never happen if both are in same tree
    }

    private void generatePermutations(int n, int k, List<Integer> current, boolean[] used, List<List<Integer>> result) {
        if (current.size() == k) {
            result.add(new ArrayList<>(current));
            return;
        }

        for (int i = 0; i < n; i++) {
            if (!used[i]) {
                used[i] = true;
                current.add(i);
                generatePermutations(n, k, current, used, result);
                current.remove(current.size() - 1);
                used[i] = false;
            }
        }
    }

    private double logSumExp(List<Double> logValues) {
        if (logValues.isEmpty()) return Double.NEGATIVE_INFINITY;
        double max = Collections.max(logValues);
        double sum = 0.0;
        for (double val : logValues) sum += Math.exp(val - max);
        return max + Math.log(sum);
    }

    private double logDiffExp(double a, double b) {
        if (b > a) throw new IllegalArgumentException("logDiffExp: b must be <= a");
        if (a == b) return Double.NEGATIVE_INFINITY;
        return a + Math.log1p(-Math.exp(b - a));
    }

    private boolean isNested(CalibrationPoint inner, CalibrationPoint outer, TreeInterface tree) {
        Node innerMRCA = getMRCA(tree, inner.taxa().asStringList());
        Node outerMRCA = getMRCA(tree, outer.taxa().asStringList());
        return isDescendant(innerMRCA, outerMRCA);
    }

    private boolean isDescendant(Node node, Node ancestor) {
        while (node != null) {
            if (node == ancestor) return true;
            node = node.getParent();
        }
        return false;
    }

    private Map<CalibrationPoint, List<CalibrationPoint>> buildNestingDAG(List<CalibrationPoint> calibrations, TreeInterface tree) {
        Map<CalibrationPoint, List<CalibrationPoint>> graph = new HashMap<>();
        for (CalibrationPoint c1 : calibrations) {
            graph.putIfAbsent(c1, new ArrayList<>());
            for (CalibrationPoint c2 : calibrations) {
                if (c1 == c2) continue;
                if (isNested(c2, c1, tree)) {
                    graph.get(c1).add(c2);  // c1 contains c2
                }
            }
        }
        return graph;
    }

    private List<CalibrationPoint> postOrderTopologicalSort(TreeInterface tree, List<CalibrationPoint> calibrations) {
        Map<CalibrationPoint, List<CalibrationPoint>> graph = buildNestingDAG(calibrations, tree);
        List<CalibrationPoint> sorted = new ArrayList<>();
        Set<CalibrationPoint> visited = new HashSet<>();

        for (CalibrationPoint c : calibrations) {
            postOrderDFS(c, graph, visited, sorted, tree);
        }

        return sorted;
    }

    private void postOrderDFS(CalibrationPoint current, Map<CalibrationPoint, List<CalibrationPoint>> graph,
                              Set<CalibrationPoint> visited, List<CalibrationPoint> result, TreeInterface tree) {
        if (visited.contains(current)) return;
        visited.add(current);

        List<CalibrationPoint> children = graph.getOrDefault(current, new ArrayList<>());
        children.sort(Comparator.comparingDouble(
                (CalibrationPoint c) -> getMRCA(tree, c.taxa().asStringList()).getHeight()
        ).reversed()); // Oldest to youngest

        for (CalibrationPoint child : children) {
            postOrderDFS(child, graph, visited, result, tree);
        }

        result.add(current); // Post-order: visit after children
    }

    public static void main(String[] args){
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        BirthDeathModel birthDeath = new BirthDeathModel();

        birthDeath.initByName("birthRate", new RealParameter("2.0"),
                "deathRate", new RealParameter("2.0"),
                "rho", new RealParameter("0.1"),
                "conditionOnRoot", true,
                "origin", new RealParameter("3.0")
        );

        boolean b;
        b = birthDeath.birthRateInput.get().getValue() - birthDeath.deathRateInput.get().getValue() < 1e-10;

        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeath
                );
        System.out.println("tree = " + tree);
        System.out.println("isCritical = " + b);
        System.out.println("logP = " + cpp.calculateLogP());
    }
}
