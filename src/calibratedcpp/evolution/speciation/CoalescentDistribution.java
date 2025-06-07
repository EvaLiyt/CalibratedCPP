package calibratedcpp.evolution.speciation;

import beast.base.inference.Distribution;
import beast.base.inference.State;

import java.util.List;
import java.util.Random;

public abstract class CoalescentDensity extends Distribution implements CoalescentDensityInterface {

    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public double calculateLogP(){
        return logP;
    }

    public double calculateLogDensity(double time) {
        return 0;
    }

    public double calculateLogCDF(double time) {
        return 0;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
