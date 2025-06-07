package calibratedcpp.evolution.speciation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Node age distribution for the CPP representation of the birth-death process")
public class BirthDeathCoalescentDistribution extends CoalescentDistribution {
    public Input<RealParameter> birthRateInput =
            new Input<RealParameter>("birthRate","the birth rate",(RealParameter)null);
    public Input<RealParameter> deathRateInput =
            new Input<RealParameter>("deathRate","the death rate",(RealParameter)null);
    public Input<RealParameter> rhoInput =
            new Input<RealParameter>("rho","Probability with which each individual in the population is sampled",(RealParameter)null);

    protected Double birthRate;
    protected Double deathRate;
    protected Double rho;
    protected Double diversificationRate;

    protected Double logBirthRate;
    protected Double logDeathRate;
    protected Double logRho;
    protected Double logDiversificationRate;

    @Override
    public void initAndValidate() {
        birthRate = birthRateInput.get().getValue();
        deathRate = deathRateInput.get().getValue();
        rho = rhoInput.get().getValue();
        diversificationRate = birthRate - deathRate;

        logBirthRate = Math.log(birthRate);
        logDeathRate = Math.log(deathRate);
        logRho = Math.log(rho);
        logDiversificationRate = Math.log(diversificationRate);
    }

    @Override
    public double calculateLogDensity(double time){
        return logRho + logBirthRate + logDiversificationRate - diversificationRate * time
        - 2 * Math.log(rho * birthRate + (birthRate * (1 - rho) - deathRate) * Math.exp(- diversificationRate * time));
    }

    @Override
    public double calculateLogCDF(double time){
        return logRho + logBirthRate + Math.log(1 - Math.exp(-diversificationRate * time)) -
                Math.log(rho * birthRate + (birthRate * (1 - rho) - deathRate) * Math.exp(-diversificationRate * time));
    }
}
