package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Node age distribution for the CPP representation of the birth-death process")
public class BirthDeathModel extends CoalescentPointProcessModel {
    public Input<RealParameter> birthRateInput =
            new Input<>("birthRate","the birth rate", (RealParameter) null);

    public Input<RealParameter> deathRateInput =
            new Input<>("deathRate","the death rate", (RealParameter) null);

    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "Diversification rate lambda - mu", (RealParameter) null);

    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "Reproductive number lambda / mu", (RealParameter) null);

    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "Turnover mu / lambda", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho","Probability with which each individual in the population is sampled",(RealParameter)null);

    protected Double birthRate;
    protected Double deathRate;
    protected Double diversificationRate;
    protected Double reproductiveNumber;
    protected Double turnover;
    protected Double rho;

    protected double logBirthRate;
    protected double logDeathRate;
    protected double logRho;
    protected double logDiversificationRate;

    protected boolean isCritical;

    @Override
    public void initAndValidate() {
        rho = rhoInput.get().getValue();

        birthRate = birthRateInput.get().getValue();
        deathRate = deathRateInput.get().getValue();
        diversificationRate = diversificationRateInput.get().getValue();
        reproductiveNumber = reproductiveNumberInput.get().getValue();
        turnover = turnoverInput.get().getValue();

        int specified = countNonNull(birthRate, deathRate, diversificationRate, reproductiveNumber, turnover);

        if (specified != 2) {
            throw new IllegalArgumentException("Exactly TWO of {birthRate, deathRate, diversificationRate, reproductiveNumber, turnover} must be specified.");
        }

        // disallow repNumber + turnover
        if (reproductiveNumber != null && turnover != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover together.");
        }

        // determine parameters from the valid combinations
       if (birthRate != null && diversificationRate != null) {
            deathRate = birthRate - diversificationRate;
        } else if (deathRate != null && diversificationRate != null) {
            birthRate = deathRate + diversificationRate;
        } else if (birthRate != null && reproductiveNumber != null) {
            deathRate = birthRate / reproductiveNumber;
        } else if (deathRate != null && reproductiveNumber != null) {
            birthRate = deathRate * reproductiveNumber;
        } else if (birthRate != null && turnover != null) {
            deathRate = birthRate * turnover;
        } else if (deathRate != null && turnover != null) {
            birthRate = deathRate / turnover;
        } else if (diversificationRate != null && reproductiveNumber != null) {
            deathRate = diversificationRate / (reproductiveNumber - 1);
            birthRate = deathRate * reproductiveNumber;
        } else if (diversificationRate != null && turnover != null) {
            birthRate = diversificationRate / (1 - turnover);
            deathRate = birthRate * turnover;
        } else {
            throw new IllegalArgumentException("Unsupported parameter combination.");
        }

        // Validation
        if (birthRate <= 0.0 || deathRate < 0.0 || rho < 0.0 || rho > 1.0) {
            throw new IllegalArgumentException("birthRate must be > 0, deathRate must be >= 0, and rho must be between 0 and 1.");
        }

        isCritical = birthRate.equals(deathRate);

        diversificationRate = birthRate - deathRate;

        logBirthRate = Math.log(birthRate);
        logDeathRate = Math.log(deathRate);
        logRho = Math.log(rho);
        logDiversificationRate = Math.log(diversificationRate);
    }

    @Override
    public double calculateLogDensity(double time){
        double logDensity;
        if (isCritical) {
            logDensity = logRho + logBirthRate - 2 * Math.log(1 + rho * birthRate * time);
        }
        else {
            logDensity = logRho + logBirthRate + 2 * logDiversificationRate - diversificationRate * time
                    - 2 * Math.log(rho * birthRate + (birthRate * (1 - rho) - deathRate) * Math.exp(-diversificationRate * time));
        }
        return logDensity;
    }

    @Override
    public double calculateLogCDF(double time){
        double logCDF;
        if (isCritical) {
            logCDF = logRho + logBirthRate + Math.log(time) - Math.log(1 + rho * birthRate * time);
        }
        else{
            logCDF = logRho + logBirthRate + Math.log(1 - Math.exp(-diversificationRate * time)) -
                    Math.log(rho * birthRate + (birthRate * (1 - rho) - deathRate) * Math.exp(-diversificationRate * time));;
        }
        return logCDF;
    }

    private int countNonNull(Double... values) {
        int count = 0;
        for (Double v : values) {
            if (v != null) count++;
        }
        return count;
    }
}
