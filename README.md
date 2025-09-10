# The Calibrated Coalescent Point Process
This is an implementation of the Calibrated Coalescent Point Process (Calibrated CPP) in BEAST2.

The CPP is a model of ultrametric trees where node ages are i.i.d. random variables. This captures a general class of birth-death processes with time-dependent birth and death rates (as in BDSKY) and age dependent death rates.
- Individuals die at some rate $\mu(x,t)$ that depends on the age of the individual $x$ and the time $t$.
- Individuals give birth to new individuals at a rate $\lambda(t)$ depending only on time.

[Lambert and Stadler (2012)](https://doi.org/10.1016/j.tpb.2013.10.002) show that these models give a uniform distribution over ranked labelled (or oriented) tree topoligies, and that the node ages are i.i.d. random variables. The node age density where $\lambda=\lambda(t)$ and $\mu=\mu(x,t)$ and each individual is sampled with probability $rho$ is,
$$q(t) = \frac{\rho\lambda (\lambda-\mu)e^{-(\lambda-\mu)}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}$$
the cumulative distribution function is,
$$\frac{\rho\lambda(1-e^{-(\lambda-\mu)t})}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{(\lambda-\mu)t}}.$$

## Structure

The abstract class CoalescentPointProcessModel has abstract methods for the density and of the node age.

BirthDeathModel extends CoalescentPointProcessModel with node age density and CDF for the birth-death process.

CalibratedCoalescentPointProcess extends SpeciesTreeDistribution and takes a CoalescentPointProcessModel and a list of calibrations, and the origin age OR conditionOnRoot as input.

