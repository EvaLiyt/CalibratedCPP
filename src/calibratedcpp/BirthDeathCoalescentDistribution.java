package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Node age distribution for the CPP representation of the birth-death process")
public class BirthDeathCoalescentDistribution extends CoalescentDistribution {
    public Input<RealParameter> birthRateInput =
            new Input<RealParameter>("birthRate","the birth rate", (RealParameter) null);

    public Input<RealParameter> deathRateInput =
            new Input<RealParameter>("deathRate","the death rate", (RealParameter) null);

    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "Diversification rate lambda - mu", (RealParameter) null);

    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "Reproductive number lambda / mu", (RealParameter) null);

    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "Turnover mu / lambda", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<RealParameter>("rho","Probability with which each individual in the population is sampled",(RealParameter)null);

    protected Double birthRate;
    protected Double deathRate;
    protected Double diversificationRate;
    protected Double reproductiveNumber;
    protected Double turnover;
    protected Double rho;

    protected Double logBirthRate;
    protected Double logDeathRate;
    protected Double logRho;
    protected Double logDiversificationRate;

    @Override
    public void initAndValidate() {
        rho = rhoInput.get().getValue();

        Double birth = getValue(birthRateInput);
        Double death = getValue(deathRateInput);
        Double div = getValue(diversificationRateInput);
        Double rep = getValue(reproductiveNumberInput);
        Double turn = getValue(turnoverInput);

        int specified = countNonNull(birth, death, div, rep, turn);

        if (specified != 2) {
            throw new IllegalArgumentException("Exactly TWO of {birthRate, deathRate, diversificationRate, reproductiveNumber, turnover} must be specified.");
        }

        // disallow repNumber + turnover
        if (rep != null && turn != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover together.");
        }

        // determine parameters from the valid combinations
        if (birth != null && death != null) {
            birthRate = birth;
            deathRate = death;
        } else if (birth != null && div != null) {
            birthRate = birth;
            deathRate = birth - div;
        } else if (death != null && div != null) {
            deathRate = death;
            birthRate = death + div;
        } else if (birth != null && rep != null) {
            birthRate = birth;
            deathRate = birth / rep;
        } else if (death != null && rep != null) {
            deathRate = death;
            birthRate = death * rep;
        } else if (birth != null && turn != null) {
            birthRate = birth;
            deathRate = birth * turn;
        } else if (death != null && turn != null) {
            deathRate = death;
            birthRate = death / turn;
        } else if (div != null && rep != null) {
            deathRate = div / (rep - 1);
            birthRate = deathRate * rep;
        } else if (div != null && turn != null) {
            birthRate = div / (1 - turn);
            deathRate = birthRate * turn;
        } else {
            throw new IllegalArgumentException("Unsupported parameter combination.");
        }

        // Validation
        if (birthRate <= 0.0 || deathRate < 0.0) {
            throw new IllegalArgumentException("birthRate must be > 0 and deathRate must be >= 0.");
        }

        diversificationRate = birthRate - deathRate;
        reproductiveNumber = birthRate / deathRate;
        turnover = deathRate / birthRate;

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

    private Double getValue(Input<RealParameter> input) {
        return input.get() != null ? input.get().getValue() : null;
    }

    private int countNonNull(Double... values) {
        int count = 0;
        for (Double v : values) {
            if (v != null) count++;
        }
        return count;
    }
}
