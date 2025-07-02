package calibratedcpp;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Abstract class for the distribution of node ages in a CPP")
public abstract class CoalescentPointProcessModel extends BEASTObject {

    public Input<RealParameter> originInput =
            new Input<>("origin", "Age of the origin (time of process start)", (RealParameter) null);

    public Input<Boolean> conditionOnRootInput =
            new Input<>("conditionOnRoot", "Whether the model is conditioned on the root age (default: false)", false);

    public abstract double calculateLogDensity(double time);
    public abstract double calculateLogCDF(double time);

    public Double origin;
    public boolean conditionOnRoot;

    @Override
    public void initAndValidate() {
        RealParameter originParam = originInput.get();
        origin = (originParam != null) ? originParam.getValue() : null;
        conditionOnRoot = conditionOnRootInput.get();

        if (origin == null && !conditionOnRoot) {
            throw new IllegalArgumentException("You must either provide an origin age or set conditionOnRoot=true.");
        }
    }

    public Double getOrigin() {
        return origin;
    }

    public boolean isConditionedOnRoot() {
        return conditionOnRoot;
    }
}