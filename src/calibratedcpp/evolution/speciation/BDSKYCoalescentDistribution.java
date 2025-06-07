package calibratedcpp.evolution.speciation;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

public class BDSKYCoalescentDensity extends CoalescentDistribution {
    Input<RealParameter> birthRateInput =
            new Input<RealParameter>("birthRate","the birth rate",(RealParameter)null);
    Input<RealParameter> deathRateInput =
            new Input<RealParameter>("deathRate","the death rate",(RealParameter)null);
    Input<RealParameter> rhoInput =
            new Input<RealParameter>("rho","the probability with which each individual in the total population is sampled",(RealParameter)null);
    Input<RealParameter> birthRateChangeTimesInput =
            new Input<RealParameter>("birthRateChangeTimes","the birth rate change times",(RealParameter)null);
    Input<RealParameter> deathRateChangeTimesInput =
            new Input<RealParameter>("deathRateChangeTimes","the death rate change times",(RealParameter)null);

    protected Double birthRate;
    protected Double deathRate;
    protected Double rho;
    protected Double birthRateChangeTimes;
    protected Double deathRateChangeTimes;

    @Override
    public void initAndValidate() {

        birthRate = birthRateInput.get().getValue();
        deathRate = deathRateInput.get().getValue();
        rho = rhoInput.get().getValue();
        birthRateChangeTimes = birthRateChangeTimesInput.get().getValue();
        deathRateChangeTimes = deathRateChangeTimesInput.get().getValue();

        super.initAndValidate();
    }

    @Override
    public double calculateLogDensity(double time){
        return 0;
    }

    @Override
    public double calculateLogCDF(double time){
        return 0;
    }

    @Override
    public double calculateLogP(){
        return logP;
    }
}
