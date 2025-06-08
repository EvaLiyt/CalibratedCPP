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
import java.util.List;
import java.util.Collections;

/**
 * @author Marcus Overwater
 */

@Description("A general class of birth-death processes with homochronous sampling conditioning on clade calibrations")
public class CalibratedCoalescentPointProcess extends SpeciesTreeDistribution {
    public Input<RealParameter> originInput =
            new Input<>("origin", "Age of the origin", (RealParameter) null);

    public Input<RealParameter> rootInput =
            new Input<>("rootAge", "Optional, if provided, the tree density is conditioned on the root age", (RealParameter) null);

    public Input<CoalescentPointProcessModel> coalescentDensityInput =
            new Input<>("treeModel", "coalescentDensity", (CoalescentPointProcessModel) null);

    public Input<List<CalibrationPoint>> calibrationsInput =
            new Input<>("calibrations","Clade calibrations", (List<CalibrationPoint>) null);

    public Input<Boolean> conditionOnCalibrationInput =
            new Input<>("conditionOnCalibrations","Boolean if the likelihood is conditioned on the clade calibrations, default true", true);

    protected CoalescentPointProcessModel q;
    protected Double origin;
    protected Double rootAge;
    protected List<CalibrationPoint> calibrations = new ArrayList<>();
    protected boolean conditionOnRoot;
    protected TreeInterface tree;
    protected double maxTime;
    protected boolean conditionOnCalibrations;

    @Override
    public void initAndValidate() {

        q = coalescentDensityInput.get();
        origin = originInput.get().getValue();
        rootAge = rootInput.get().getValue();
        calibrations = calibrationsInput.get();
        tree = treeInput.get();

        conditionOnCalibrations = (calibrations != null) ? conditionOnCalibrationInput.get() : false;

        conditionOnRoot = (rootAge != null);

        if (origin == null && rootAge == null) {
            throw new IllegalArgumentException("Either root age or origin must be specified.");
        }

        maxTime = conditionOnRoot ? rootAge : origin;

        super.initAndValidate();
    }

    public double calculateUnConditionedTreeLogLikelihood(TreeInterface tree) {

        double logP = Math.log1p(-Math.exp(q.calculateLogCDF(maxTime)));

        if (conditionOnRoot) {
            logP += logP - q.calculateLogDensity(maxTime);
        }

        for (Node node : tree.getInternalNodes()){
            double age = node.getHeight();
            logP += q.calculateLogDensity(age);
        }

        return logP;
    }

    public double calculateLogMarginalDensityOfCalibrations(TreeInterface tree, List<CalibrationPoint> calibrations) {
        double marginalDensity = Math.log1p(-Math.exp(q.calculateLogCDF(maxTime)));

        if (conditionOnRoot) {
            marginalDensity *= 2;
        }

        int numCalibrations = calibrations.size();

        int[] cladeSize = new int[numCalibrations];

        double[] calibrationAges = new double[numCalibrations];

        int numTaxa = tree.getLeafNodeCount();

        for (int i = 0 ; i < numCalibrations ; i++) {
            CalibrationPoint calibration = calibrations.get(i);

            calibrationAges[i] = getMRCA(tree, calibration.taxa().asStringList()).getHeight();

            cladeSize[i] = calibration.taxa().getTaxonSet().size();

            marginalDensity += q.calculateLogDensity(calibrationAges[i]) + (cladeSize[i] - 2) * q.calculateLogCDF(calibrationAges[i]);
        }

        marginalDensity += calculateSumOfPermutations(calibrationAges, numCalibrations, cladeSize, numTaxa);

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

    private double calculateSumOfPermutations(double[] cladeAges, int numCalibrations, int[] cladeSizes, int numTaxa) {
        int totalConditionedTaxa = 0;
        for (int c : cladeSizes) {
            totalConditionedTaxa += c;
        }
        int m = numTaxa - totalConditionedTaxa;

        return 0.0;
    }

    private double logSumExp(List<Double> logValues) {
        if (logValues.isEmpty()) return Double.NEGATIVE_INFINITY;
        double max = Collections.max(logValues);
        double sum = 0.0;
        for (double val : logValues) sum += Math.exp(val - max);
        return max + Math.log(sum);
    }

    public static void main(String[] args){
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        BirthDeathModel birthDeath = new BirthDeathModel();

        birthDeath.initByName("deathRate", new RealParameter("1.0"),
                "diversificationRate", new RealParameter("0.0"),
                "rho", new RealParameter("0.1")
        );

        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeath,
                "origin", new RealParameter("5.0")
                );

        System.out.println(tree);
        System.out.println("logP=" + cpp.calculateLogP());
    }
}