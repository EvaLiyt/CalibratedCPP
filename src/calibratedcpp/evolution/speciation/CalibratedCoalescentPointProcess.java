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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marcus Overwater
 */

@Description("A general class of birth-death processes with homochronous sampling clade calibrations")

public class CalibratedCoalescentPointProcess extends SpeciesTreeDistribution {
    public Input<RealParameter> originInput =
            new Input<RealParameter>("origin","Age of the origin",(RealParameter) null);
    public Input<RealParameter> rootInput =
            new Input<RealParameter>("rootAge","Age of the root",(RealParameter) null);
    public Input<CoalescentDistribution> coalescentDensityInput =
            new Input<CoalescentDistribution>("treeModel", "coalescentDensity", (CoalescentDistribution) null);
    public Input<List<CalibrationPoint>> calibrationsInput =
            new Input<List<CalibrationPoint>>("calibrations","Clade calibrations", (List<CalibrationPoint>) null);
    public Input<Boolean> conditionOnRootInput =
            new Input<Boolean>("conditionOnRoot","Condition on the root age", false);

    protected CoalescentDistribution q;
    protected Double origin;
    protected Double rootAge;
    protected List<CalibrationPoint> calibrations = new ArrayList<>();
    protected Boolean conditionOnRoot;
    protected TreeInterface tree;
    protected double maxTime;

    @Override
    public void initAndValidate() {
        q = coalescentDensityInput.get();
        origin = originInput.get().getValue();
        calibrations = calibrationsInput.get();
        conditionOnRoot = conditionOnRootInput.get();
        tree = treeInput.get();
        maxTime = origin;

        if (conditionOnRoot) {
            rootAge = rootInput.get().getValue();
            maxTime = rootAge;
        }

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

    public double calculateLogMarginalDensityOfCalibrations() {
        double marginalDensity = Math.log1p(-Math.exp(q.calculateLogCDF(maxTime)));

        if (conditionOnRoot)
            marginalDensity *= 2;

        for (CalibrationPoint calibration : calibrations) {
            double calibrationAge = getMRCA(tree, calibration.taxa().asStringList()).getHeight();
            int calibrationSize = calibration.taxa().getTaxonSet().size();
            marginalDensity += q.calculateLogDensity(calibrationAge) + (calibrationSize - 2) * q.calculateLogCDF(calibrationAge);
        }

        return marginalDensity;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        logP = calculateUnConditionedTreeLogLikelihood(tree);

        if (calibrations!=null) {
            logP -= calculateLogMarginalDensityOfCalibrations();
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

    public static void main(String[] args){
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        BirthDeathCoalescentDistribution birthDeath = new BirthDeathCoalescentDistribution();

        birthDeath.initByName("birthRate", new RealParameter("2.0"),
                "deathRate", new RealParameter("1.0"),
                "rho", new RealParameter("0.1")
        );

        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "model", birthDeath,
                "origin", new RealParameter("5.0")
                );

        System.out.println(tree);
        System.out.println("logP=" + cpp.calculateLogP());
        System.out.println("origin=" + cpp.origin);
    }
}