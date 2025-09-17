# The Calibrated Coalescent Point Process
[![Unit/integration tests](https://github.com/moverwater/CalibratedCPP/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/moverwater/CalibratedCPP/actions/workflows/main.yml)

This is an implementation of the Calibrated Coalescent Point Process (Calibrated CPP) in BEAST2.

The CPP is a model of ultrametric trees where node ages are i.i.d. random variables. This captures a general class of birth-death processes with time-dependent birth and death rates (as in BDSKY) and age dependent death rates.
- Individuals die at some rate $\mu(x,t)$ that depends on the age of the individual $x$ and the time $t$.
- Individuals give birth to new individuals at a rate $\lambda(t)$ depending only on time.

[Lambert and Stadler (2012)](https://doi.org/10.1016/j.tpb.2013.10.002) show that these models give a uniform distribution over ranked labelled (or oriented) tree topologies, and that the node ages are i.i.d. random variables. The node age density where $\lambda=\lambda(t)$ and $\mu=\mu(x,t)$ and each individual is sampled with probability $\rho$ is,

$$q(t) = \frac{\rho\lambda (\lambda-\mu)e^{-(\lambda-\mu)}}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}$$

the cumulative distribution function is,

$$Q(t) = \frac{\rho\lambda(1-e^{-(\lambda-\mu)t})}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}.$$

The Calibrated CPP is a calibrated tree prior using the CPP. Calibrated tree priors are used for molecular clock dating by conditioning on the existence and ages of the most recent common ancestors of monophyletic clades.

## Structure

The implementation has the following structure:
- The abstract class CoalescentPointProcessModel has abstract methods for the density and CDF of the node age.
- BirthDeathModel extends CoalescentPointProcessModel with node age density and CDF for the constant rate birth-death process.
- CalibratedCoalescentPointProcess extends SpeciesTreeDistribution and takes a CoalescentPointProcessModel, a list of calibrations, and the origin age OR conditionOnRoot as inputs.

## License

CalibratedCPP is free software.  It is distributed under the terms of version 3 of the GNU General Public License.  A copy of this license should be found in the file [COPYING](./COPYING) located in the root directory of this repository. If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.










