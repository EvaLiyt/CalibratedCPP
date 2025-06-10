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
            new Input<>("conditionOnCalibrations","Boolean if the likelihood is conditioned on the clade calibrations (Default: true)", true);

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

        tree = treeInput.get();
        model = cppModelInput.get();
        calibrations = calibrationsInput.get();
        conditionOnCalibrations = (calibrations != null) ? conditionOnCalibrationInput.get() : false;

        if (conditionOnCalibrations) {
            calibrations.sort(Comparator.comparingDouble((CalibrationPoint calibration) ->
                    getMRCA(tree, calibration.taxa().asStringList()).getHeight()
            ).reversed());
        }

        origin = model.getOrigin();

        conditionOnRoot = model.isConditionOnRoot();

        rootAge = tree.getRoot().getHeight();

        if (!conditionOnRoot && origin < rootAge) {
            throw new RuntimeException("rootAge(" + rootAge + ") > origin(" + origin + "): Origin must be greater than root age.");
        }
        if(conditionOnRoot && origin != null){
            System.err.println("WARNING: Overriding origin when conditionOnRoot is true");
        }

        maxTime = conditionOnRoot ? rootAge : origin;

        super.initAndValidate();
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

        int numCalibrations = calibrations.size();

        int[] cladeSize = new int[numCalibrations];

        double[] calibrationAges = new double[numCalibrations];

        double logQt = model.calculateLogCDF(maxTime);

        int numTaxa = tree.getLeafNodeCount();

        for (int i = 0 ; i < numCalibrations ; i++) {
            CalibrationPoint calibration = calibrations.get(i);

            calibrationAges[i] = getMRCA(tree, calibration.taxa().asStringList()).getHeight();

            cladeSize[i] = calibration.taxa().getTaxonSet().size();

            marginalDensity += model.calculateLogDensity(calibrationAges[i]) + (cladeSize[i] - 2) * model.calculateLogCDF(calibrationAges[i]);
        }

        marginalDensity += calculateLogSumOfPermutations(calibrationAges, numCalibrations, cladeSize, numTaxa, logQt);

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

    private double calculateLogSumOfPermutations(double[] cladeAges, int numCalibrations, int[] cladeSizes, int numTaxa, double logQ_t) {
        int c = 0;
        for (int size : cladeSizes) c += size;
        int m = numTaxa - c;

        // Precompute log Q(t_i) and log q(t_i)
        double[] logQ_ti = new double[numCalibrations];
        double[] logDiff = new double[numCalibrations];

        for (int i = 0; i < numCalibrations; i++) {
            logQ_ti[i] = model.calculateLogCDF(cladeAges[i]);
            logDiff[i] = logDiffExp(logQ_t, logQ_ti[i]); // Precompute log difference of Q(t) - Q(t_i)
        }

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

// Same generatePermutations and logSumExp as before

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

    public static void main(String[] args){
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        BirthDeathModel birthDeath = new BirthDeathModel();

        birthDeath.initByName("birthRate", new RealParameter("2.0"),
                "deathRate", new RealParameter("2.0"),
                "rho", new RealParameter("0.1"),
                "conditionOnRootAge", true,
                "origin", new RealParameter("0.0")
        );
        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeath
                );
        System.out.println(tree);
        System.out.println("logP=" + cpp.calculateLogP());
    }
}
