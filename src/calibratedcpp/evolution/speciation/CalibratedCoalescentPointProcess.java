package coalescentcpp.evolution.speciation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marcus Overwater
 */

@Description("A general class of birth-death processes with homochronous sampling clade calibrations")

public class CalibratedCoalescentPointProcess extends SpeciesTreeDistribution {
    Input<RealParameter> originInput =
            new Input<RealParameter>("origin","Age of the origin",(RealParameter) null);
    Input<RealParameter> rootInput =
            new Input<RealParameter>("rootAge","Age of the root",(RealParameter) null);
    Input<Distribution> coalescentDensityInput =
            new Input<Distribution>("model", "coalescentDensity", Distribution.class);
    Input<List<CalibrationPoint>> calibrationsInput =
            new Input<List<CalibrationPoint>>("calibrations","Clade calibrations", CalibrationPoint.class);
    Input<Boolean> conditionOnRootInput =
            new Input<Boolean>("conditionOnRoot","Condition on the root age",(Boolean) null);

    protected CoalescentDensity q;
    protected Double origin;
    protected Double rootAge;
    protected List<CalibrationPoint> calibrations = new ArrayList<>();
    protected Boolean conditionOnRoot;
    protected TreeInterface tree;

    @Override
    public void initAndValidate() {
        q = (CoalescentDensity) coalescentDensityInput.get();
        origin = originInput.get().getValue();
        rootAge = rootInput.get().getValue();
        calibrations = calibrationsInput.get();
        conditionOnRoot = conditionOnRootInput.get();
        tree = treeInput.get();

        super.initAndValidate();
    }

    public double calculateUnConditionedTreeLogLikelihood(TreeInterface tree) {
        double logP = q.calculateLogCDF(origin);
        for (Node node : tree.getInternalNodes()){
            double age = node.getHeight();
            logP += q.calculateLogDensity(age);
        }
        return logP;
    }

    public double calculateLogMarginalDensityOfCalibrations(){
        for (CalibrationPoint calibrationPoint : calibrations){
            double mrca = tree.getCo
        }
        return 0;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        logP = calculateUnConditionedTreeLogLikelihood(tree) - calculateLogMarginalDensityOfCalibrations();
        return logP;
    }

    public static void main(String[] args){

    }
}