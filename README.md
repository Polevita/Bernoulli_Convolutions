# Bernoulli_Convolutions
Code in fortran for our joint paper with M.Pollicott and V. Kleptsyn "Uniform lower bounds on the dimension of Bernoulli convolutions", Advances in Mathematics, Advances in Mathematics, Volume 395, 2022, paper No 108090, https://doi.org/10.1016/j.aim.2021.108090.
(https://www.sciencedirect.com/science/article/pii/S0001870821005296). 

The repository contains the set of Fortran programs used to compute lower bounds
for the correlation dimension and regularity exponents of the invariant measure of
an iterated function scheme of similarities. There are four routines for 
Bernoulli convolutions and one routine adapted to the {0,1,3}-system. 

===========

Bernoulli Convolutions

===========

The program dimcor.f90 is used to improve an existing lower bound
on the correlation dimension of Bernoulli convolution for a large interval
of parameter values, using a given partition. 

It should be compiled with -fopenmp option. It will then use as many threads
as available. 

The upper bound on the dimension and regularity exponent for the Bernoulli 
convolution routines is assumed to be 1. 

When run, it prompts for the following parameters 
- theta for the operator \widehat A_\theta
- the total number lam_N of values of lambda and corresponding guesses to be read from the file 
- the number of first end point lambda to consider lam_0
- the number of the last end point lambda to consider lam_k, k <= n
- the number of partition intervals for the test function
- the number of iterations of the diffusion operator 
- the refinement parameter epsilon
- the file name for the output data
- the file name of the file with the initial guesses 

A sample input file: (lambda      dimension guess; dim.dat)

0.5709848000        0.9950000000
0.5709696000        0.9950000000 
0.5709544000        0.9950000000
0.5709392000        0.9950000000
0.5709240000        0.9950000000
0.5709088000        0.9950000000

It is assumed that the values of lambda are the end points of the partition in Lambda. 

The output allows the user to follow the progress of the computation. 
The result is written in the file specified by the user in the following format 
#; value of lambda; the number of refinements in alpha, the number of iterations 
of the diffusion operator, and the dimension value. 

------

The program dimcorsingle.f90 is used to compute a lower bound on correlation dimension
for the Bernoulli convolution measure for a single parameter value. 
When run, it prompts for the following parameters: 
 
- theta for the operator \widehat A_\theta 
- the length of the interval Lambda containing the parameter value lambda
- the value of lambda 
- the dimension guess 
- the number of partition intervals for the test function
- the number of iterations of the diffusion operator 
- the refinement parameter epsilon
- the file name for the output data

The output file contains the numerical values of the piecewise constant function $\psi$,
which is the image of the indicator function of an admissible interval J_\lamdba 

-------

The program regularity.f90 is used to compute a lower bound on the regularity exponent
of the invariant measure for Bernoulli convolutions. When run, it prompts the user for the same 
parameters as the function dimcor.f90. The sample input file is nums.dat.

-------

The program regularitysingle.f90 is used to compute a lower bound on the regularity exponent
for a single value of parameter lambda.  It requires the similar parameters as the function 
dimcorsingle.f90. 

==============

{0,1,3} - system

==============

The program 013systemcor.f90 is used to compute lower bound on the correlation dimension 
of the invariant measure in the {0,1,3} system. It requires the same user input as the program 
dimcor.f90, but it additionally needs the file to provide an initial guess
for the both  lower and upper bounds on the correlation dimension. 
It also makes no assumption on the length of Lambda-intervals, 
and requires the user to provide it as a parameter (the same value for all values of lambda). 

Sample input file (lambda, lower bound, upper bound -- guesses.dat): 

0.3333333333333      0.750000000       1.0000000000000000
0.3027756377319      0.801451700       0.9195230257947472
0.3049692805112      0.803376500       0.9251127356890137
0.3221853546260      0.800489200       0.9699672204962546
0.3225141365024      0.808122700       0.9708414823016051
0.3294085281925      0.754574200       0.9893338621886850
0.3294521704854      0.750000000       0.9894519043256954
0.3319890295845      0.758827500       0.9963351443598648
The archive contains the set of Fortran programs used to compute lower bounds
for the correlation dimension and regularity exponents of the invariant measure of
an iterated function scheme of similarities. There are four routines for 
Bernoulli convolutions and one routine adapted to the {0,1,3}-system. 

===========

Bernoulli Convolutions

===========

The program dimcor.f90 is used to improve an existing lower bound
on the correlation dimension of Bernoulli convolution for a large interval
of parameter values, using a given partition. 

It should be compiled with -fopenmp option. It will then use as many threads
as available. 

The upper bound on the dimension and regularity exponent for the Bernoulli 
convolution routines is assumed to be 1. 

When run, it prompts for the following parameters 
- theta for the operator \widehat A_\theta
- the total number lam_N of values of lambda and corresponding guesses to be read from the file 
- the number of first end point lambda to consider lam_0
- the number of the last end point lambda to consider lam_k, k <= n
- the number of partition intervals for the test function
- the number of iterations of the diffusion operator 
- the refinement parameter epsilon
- the file name for the output data
- the file name of the file with the initial guesses 

A sample input file: (lambda      dimension guess; dim.dat)

0.5709848000        0.9950000000
0.5709696000        0.9950000000 
0.5709544000        0.9950000000
0.5709392000        0.9950000000
0.5709240000        0.9950000000
0.5709088000        0.9950000000

It is assumed that the values of lambda are the end points of the partition in Lambda. 

The output allows the user to follow the progress of the computation. 
The result is written in the file specified by the user in the following format 
#; value of lambda; the number of refinements in alpha, the number of iterations 
of the diffusion operator, and the dimension value. 

------

The program dimcorsingle.f90 is used to compute a lower bound on correlation dimension
for the Bernoulli convolution measure for a single parameter value. 
When run, it prompts for the following parameters: 
 
- theta for the operator \widehat A_\theta 
- the length of the interval Lambda containing the parameter value lambda
- the value of lambda 
- the dimension guess 
- the number of partition intervals for the test function
- the number of iterations of the diffusion operator 
- the refinement parameter epsilon
- the file name for the output data

The output file contains the numerical values of the piecewise constant function $\psi$,
which is the image of the indicator function of an admissible interval J_\lamdba 

-------

The program regularity.f90 is used to compute a lower bound on the regularity exponent
of the invariant measure for Bernoulli convolutions. When run, it prompts the user for the same 
parameters as the function dimcor.f90. The sample input file is nums.dat.

-------

The program regularitysingle.f90 is used to compute a lower bound on the regularity exponent
for a single value of parameter lambda.  It requires the similar parameters as the function 
dimcorsingle.f90. 

==============

{0,1,3} - system

==============

The program 013systemcor.f90 is used to compute lower bound on the correlation dimension 
of the invariant measure in the {0,1,3} system. It requires the same user input as the program 
dimcor.f90, but it additionally needs the file to provide an initial guess
for the both  lower and upper bounds on the correlation dimension. 
It also makes no assumption on the length of Lambda-intervals, 
and requires the user to provide it as a parameter (the same value for all values of lambda). 

Sample input file (lambda, lower bound, upper bound -- guesses.dat): 

0.3333333333333      0.750000000       1.0000000000000000
0.3027756377319      0.801451700       0.9195230257947472
0.3049692805112      0.803376500       0.9251127356890137
0.3221853546260      0.800489200       0.9699672204962546
0.3225141365024      0.808122700       0.9708414823016051
0.3294085281925      0.754574200       0.9893338621886850
0.3294521704854      0.750000000       0.9894519043256954
0.3319890295845      0.758827500       0.9963351443598648

