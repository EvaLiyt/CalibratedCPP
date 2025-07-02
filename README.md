# CalibratedCPP_BEAST
Implementation of the Calibrated Coalescent Point Process (Calibrated CPP) in BEAST2.

The CPP is a model of ultrametric trees where node ages are i.i.d. random variables. This captures a general class of birth-death processes with time-dependent birth and death rates (as in BDSKY) and age dependent death rates.

## Structure

The abstract class CoalescentPointProcessModel has abstract methods for the density and of the node age.

BirthDeathModel extends CoalescentPointProcessModel with node age density and CDF for the birth-death process.

CalibratedCoalescentPointProcess extends SpeciesTreeDistribution and takes a CoalescentPointProcessModel and a list of calibrations, and the origin age OR conditionOnRoot as input.
