package calibratedcpp.evolution.speciation;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

public class BirthDeathCoalescentDensity extends CoalescentDensity {
    Input<RealParameter> birthRateInput =
            new Input<>("birthRate","the birth rate",(RealParameter)null);
    Input<RealParameter> deathRateInput =
            new Input<>("deathRate","the death rate",(RealParameter)null);
    Input<RealParameter> rhoInput =
            new Input<>("rho","Probability with which each individual in the population is sampled",(RealParameter)null);

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

        super.initAndValidate();
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
