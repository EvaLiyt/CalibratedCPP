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
            new Input<>("birthRate","The birth rate (lambda)", (RealParameter) null);

    public Input<RealParameter> deathRateInput =
            new Input<>("deathRate","The death rate (mu)", (RealParameter) null);

    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "Diversification rate (lambda - mu)", (RealParameter) null);

    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "Reproductive number (lambda / mu)", (RealParameter) null);

    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "Turnover (mu / lambda)", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho","Probability with which each individual in the population is sampled.", (RealParameter) null );

    protected Double birthRate;
    protected Double deathRate;
    protected Double diversificationRate;
    protected Double reproductiveNumber;
    protected Double turnover;
    protected Double rho;

    protected double logBirthRate;
    protected double logDeathRate;
    protected double logDiversificationRate;
    protected double logRho;

    protected double A;
    protected double B;

    @Override
    public void initAndValidate() {
        birthRate = safeGet(birthRateInput);
        deathRate = safeGet(deathRateInput);
        diversificationRate = safeGet(diversificationRateInput);
        reproductiveNumber = safeGet(reproductiveNumberInput);
        turnover = safeGet(turnoverInput);
        rho = safeGet(rhoInput);

        int specified = 0;
        for (Double i : new Double[] {birthRate, deathRate, diversificationRate, reproductiveNumber, turnover}) {
            if (i != null) specified++;
        }

        if (rho == null) {
           throw new IllegalArgumentException("rho parameter must be specified.");
        }

        if (specified != 2) {
            throw new IllegalArgumentException("Exactly TWO of {birthRate, deathRate, diversificationRate, reproductiveNumber, turnover} must be specified.");
        }

        // disallow repNumber + turnover
        if (reproductiveNumber != null && turnover != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover together.");
        }

        // determine birth and death rate parameters from the valid combinations
        if (birthRate != null && deathRate != null) {
            diversificationRate = birthRate - deathRate;
        } else if (birthRate != null && diversificationRate != null) {
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
        if (birthRate <= 0.0 || deathRate < 0.0 || rho > 1.0 || rho <= 0.0) {
            throw new IllegalArgumentException("birthRate (" + birthRate + ") must be > 0, deathRate (" + deathRate + ") must be >= 0, AND rho (" + rho + ") must be between 0.0 and 1.0.");
        }

        A = rho * birthRate;
        B = birthRate * (1 - rho) - deathRate;

        diversificationRate = birthRate - deathRate;

        logBirthRate = Math.log(birthRate);
        logDeathRate = Math.log(deathRate);
        logRho = Math.log(rho);
        logDiversificationRate = Math.log(Math.abs(diversificationRate));

        super.initAndValidate();
    }

    @Override
    public double calculateLogDensity(double time){
        double logDensity;
        double rt = diversificationRate * time;

        if (diversificationRate == 0.0) {
            // Critical case
            logDensity = logRho + logBirthRate - 2 * Math.log(1 + rho * birthRate * time);
        }
        else if (diversificationRate < 0) {
                // Sub-critical case: use stable form with exp(r * t)
                logDensity = logRho + logBirthRate + 2 * logDiversificationRate + rt
                        - 2 * Math.log(Math.abs(A * Math.exp(rt) + B));
            } else {
                // Supercritical case: formula with exp(-r * t)
                logDensity = logRho + logBirthRate + 2 * logDiversificationRate - rt
                        - 2 * Math.log(A + B * Math.exp(-rt));
            }
        return logDensity;
    }

    @Override
    public double calculateLogCDF(double time){
        double logCDF;

        if (diversificationRate == 0.0) {
            // Critical case
            logCDF = logRho + logBirthRate + Math.log(time) - Math.log1p(rho * birthRate * time);
        }
        else if (diversificationRate < 0) {
                // Sub-critical case
                double exp_rt = Math.exp(diversificationRate * time); // decays, stable

                logCDF = logRho + logBirthRate
                        + Math.log1p(-exp_rt)     // log(1 - exp(r * t)) stable for r<0
                        - Math.log(- A * exp_rt - B);
            } else {
                // Supercritical case
                double exp_neg_rt = Math.exp(-diversificationRate * time);

                logCDF = logRho + logBirthRate
                        + Math.log1p(-exp_neg_rt)  // log(1 - exp(-r * t)) stable for r>=0
                        - Math.log(A + B * exp_neg_rt);
            }
        return logCDF;
    }

    private Double safeGet(Input<RealParameter> input) {
        RealParameter param = input.get();
        return (param != null) ? param.getValue() : null;
    }
}
