package calibratedcpp;

import beast.base.core.BEASTObject;
import beast.base.core.Description;

/**
 * @author Marcus Overwater
 */

@Description("Abstract class for the distribution of node ages in a CPP")
public abstract class CoalescentPointProcessModel extends BEASTObject {
    public abstract double calculateLogDensity(double time);
    public abstract double calculateLogCDF(double time);
}