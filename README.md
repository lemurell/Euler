# Experimental L-function computations
 
This repository contains Mathematica code to compute experimental L-functions.

The file Lfunctions.m contains all auxiliary functions that are called from the other files.

The main file to run is Total.m. It tries to compute all L-functions of a specified type with parameters in a specified region. This is done in three steps:
* Step 1: Scan the complete region by evaluations in a symmetric grid. Identify candidates for L-functions.
* Step 2: Try to zoom in on each of the candiates found in Step 1. Remove candidates that don't improve and remove duplicates. (Files used: Step 2, ZoomOne)
* Step 3: Zoom in further on the candidates that passed Step 2. (Files used: Step 2, ZoomOne)

The results of the computations are stored in a set of files. The ones that are estimated to be approximations of L-functions are collected in the file that ends in XXSuccess. More detailed information is stored in the files that end in XXZoom?? where ?? takes the values 1,2,3, ... for each of the candidates. The file that ends in candidates.txt contains all the candidates in Steps 2 and 3.

Before running the file Total.m a number of variables need to be set. The file LongRun.m contains an example on how to set all the necessary variables and then run Total.m.

## Variables to specify to run Total.m

The variables can be divided into three different subsets.

1. Variables specifying the type of the L-function
    * Ltype: Should be one of the values "Maass", "Holomorph", "GL3", "GL3Holo", "SP4", "GL4", "HoloxMaass", "SiegelPara", "SiegelSpin", "SP4Std", "SP6", see more inforamtion below. 
    * level: The conductor of the L-function
    * charvalue: Set to "No" if trivial character
    * parity: Speciefies the (fixed) real parts of the spectral parameters
    * OMEGA: The sign of the functional equation
    * realOrImaginary: 1 = real coefficients (self dual), 0 = complex coefficients
    * knownCoef: Normally equal to {} (if some coefficients are known beforehand, then they can be specified here)
2. Variables specifying the region to search (we suppose that the number of unknown parameters is d)
    * Rtuple: A list of length d specifying the smallest tuple of the d-dimensional box to search
    * sideLengths: A list of length d specifying the length of the sides of the d-dimensional box to search
    * RstepList: A list of length d specifying the step size in each of the d dimensions
    * sideCounts: A list of length d specifying the number of values in each dimensions, computed by getSideCounts[sideLengths, RstepList] (should be put into Total.m)
 3. Variables specifying parameters for the algorithm and where to store the results
    * testfunctiontypeSweep: Which type of testfunctions to use in the initial scan, could be either "DS" (recommended) or "Classic"
    * testfunctiontypeZoom: Which type of testfunctions to use when zooming in, could be either "DS" (recommended) or "Classic"
    * bwidthStart: Specifies range of parameters in the testfunctions. Recommended values are 1 in degree 3 and 3/2 in degree 4.
    * maxPower: Specifies the maximal power of a prime to use for the multiplicative relations between coefficients. The recommended and maximal value is 4. This specifies for example that a_16 (16=2^4) is not an unknown but evaluated from a_2 (and a_4) but a_32 (32=2^5) is an unknown.
    * If using "Classic" then the following parameters need to be set (with recommended values): initNN = 20, minWidth = 40 / (realOrImaginary + 1), minMult = 1/200
    * fileNameBase: This specifies the base of the name of the files to store the results. It's relative to the current directory. For example "Runs/Test1" would put the results in the folder "Runs" int the current directory and all the produced files would start with "Test1".

## The different types of L-functions

The type of the L-function is specified by the variable **Ltype** and the meaning of the different values are as follows:

1. **Degree 2**
   * "Maass": L-function of a degree 2 Maass form. One unknown parameter. Set parity = {0} for even Maass forms, parity = {1} for odd Maass forms.
   * "Holomorph": L-function of a classical modular form. No unknown parameter. Set parity = {k} for a weight k modular form.
2. **Degree 3**
   * "GL3": L-function with three Gamma-factors. Two unknown parameters. Set parity = {0} for r0r0r0 (spherical Maass forms), parity = {1} for r1r1r0. (Other two cases not implemented yet.)
   * "GL3Holo": L-functions with two Gamma factors. One unknown parameter. Set parity = {d, k} for rdc(k-1) (so k is a positive even integer)
3. **Degree 4**
   * "SP4": Self dual L-function with four Gamma-factors. Two unknown parameters. Set parity = {d, e} for rdrerdre.
   * "GL4": L-functions corresponding to GL4 Maass forms. Three unknown parameters. Only implemented with "trivial parity", parity = {}
   * "HoloxMaass": Self dual L-functions with three Gamma-factors. One unknown parameter. Set parity = {d, k} for rdrdc(k-1) (so k is a positive even integer)
   * "SiegelPara": TBA
   * "SiegelSpin": TBA
4. **Degree 5**
   * "SP4Std": TBA
5. **Degree 6**
   * "SP6": Self dual L-function with six Gamma-factors. Three unknown parameters. Only implemented with "trivial parity", parity = {}.

## Structure of the file containing all the found L-functions (name ends in XXSuccess)

The file that ends in XXSuccess contains a list of all results that are estimated to be an experimental L-function. Each item in the list corresponds to one L-function and is a list with the following structure:

{ {Ltype, level, charvalue, OMEGA, parity}, Rtuple, {{primeCoef} , {{a2, a3}, 1}, {list of errors in indicators}}, RMargin, coefMargin}

* The variables in the first element are explained above.
* Rtuple: This is the approximation of the imaginary parts of the spectral parameters. The normalization is so that it's half of the values in the paper, so e.g. the first GL3 Maass form will have Rtuple approximately equal to {6.8, 2.4}.
* primeCoef: This will contain approximatons of all coefficients that can't be evaluated from the other ones. So for degrees 2 and 3 it will be primes ({a2, a3, a5, a7, ...}), for degree 4 and 5 it will also contain squares of primes ({a2, a3, a4, a5, a7, a9, ...}) and for degree 6 (and 7) it will also contain cubes of primes ({a2, a3, a4, a5, a7, a8, a9, ...}).
* The list of errors in the indicators could be ignored.
* Rmargin: Gives a list of estimated error margins of the values in Rtuple. This should only be regarded as an estimation.
* coefMargin: Gives a list of estimated error margins of the values in primeCoef. This should only be regarded as an estimation.

## The file LongRun.m

This file contains an example on how to set the necessary parameters and then run the file Total.m. It's possible to do several different boxes in one run. The actual content of the file will compute the self dual L-function of type r0r0r0r0 with the smallest larger spectral parameter, i.e. this one https://www.lmfdb.org/L/4/1/1.1/r0e4/c4.72c12.47/0.

Before running the file all the Mathematica files in this folder should be downloaded and put in a local folder. In this folder one should create a subfolder Runs where generated files will end up. 

It takes approxiamately 20 minutes to run on an average computer.
