# Adaptive Time-Steppers for the Poisson-Nernst-Planck Equations

This repo contains the MATLAB code used in the papers [Adaptive Time-stepping Schemes for the Solution of the Poisson-Nernst-Planck Equations](https://arxiv.org/abs/1703.10297) and [A Study of the Numerical Stability of an ImEx 
Scheme with Application to the Poisson-Nernst-Planck Equations](https://arxiv.org/abs/1905.01368). The same time-stepper and solver was used to study voltammetry in [this paper](https://arxiv.org/abs/1608.07004), and in general is intended to solve the PNP-gFBV equations efficiently.

Unfortunately, as this work was done as part of my doctoral work which has since wrapped up, I am no longer actively maintaining it, and it is somewhat disorganized as a result of modifications over the course of my PhD. Though the code should be relatively straightforward to modify if you wish to run your own custom simulations, feel free to ask me any questions about the papers or the code by submitting an issue. Specifics for parts of the code are documented below.

## Notes for code pertaining to **Adaptive Time-stepping Schemes for the Solution of the Poisson-Nernst-Planck Equations**

The following entry points are set up to create the results in Tables 2 and 3 of [Adaptive Time-stepping Schemes for the Solution of the Poisson-Nernst-Planck Equations](https://arxiv.org/abs/1703.10297): 

- ```dy_run_semi_implicit.m``` is the entry point to run the semi-implicit (VSSBDF2) time-stepper
- ```dy_run_fully_implicit_newton.m``` is the entry point to run the fully implicit (VSBDF2) time-stepper with our hand-coded Newton-Raphson method
- ```dy_run_fully_implicit_fsolve.m``` is the entry point to run the fully implicit (VSBDF2) time-stepper with MATLAB's fsolve

For each of these scripts, you can change the value of epsilon (the nondimensionalized Debye length) at the beginning of the scripts to run with a different value, as these are intended to be sweeps in epsilon. You can turn on per-step console output by uncommenting the various ```display(sprintf('n=%d,err=%.3g,time=%.3g,dt=%.3g,c=%d',n,err,t(n),dt,c));``` lines in ```run.m```, ```run_fully_implicit.m```, and ```run_fully_implicit_fsolve.m```.

In order to recreate the convergence test results in Table 1 of [Adaptive Time-stepping Schemes for the Solution of the Poisson-Nernst-Planck Equations](https://arxiv.org/abs/1703.10297), run ```convergence_test.m```. This runs the convergence test for the SBDF2 stepping scheme (but can be modified to call ```run_fully_implicit.m``` instead of ```run.m``` to run the convergence test on BDF2). In order to use one of the "direct" boundary conditions on concentration instead of ghost points, you will have to go into ```concentration.m``` and manually comment out lines ```38-49``` and uncomment either lines ```59-67``` for a two-point stencil direct method, or lines ```71-73```, ```76-86``` and ```89-95``` for a three-point stencil.

In order to perform simulations with custom parameters, see ```parameters.m```, which contains all of the physical constants (epsilon is commented out because it is being set elsewhere). You will also need to set the simulation hyperparameters (spatial and temporal discretization parameters), and set initial conditions on c\_p and c\_m. Examples of how to do this are in the existing code. Note that smaller epsilons require a larger number of spatial mesh points to resolve.

Set ```bc=1``` for voltage boundary conditions (voltage as a function of time is set in ```voltage.m```) or ```bc=0``` for current boundary conditions (current as a function of time is set in ```cur.m```)

Finally, Richardson Extrapolation is turned off (i.e. commented out) in this code because we were testing for numerical properties rather than running optimal simulations, so the time step update is directly equal to the "coarse" values, instead of blending "coarse" and "fine" values. If you want to run this code to actually simulate the PNP-FBV equations, Richardson Extrapolation should be turned on -- extrapolation will result in higher order accuracy for the local truncation error, and ultimately faster runtimes. You can find the commented out code in the various ```run_xxx.m``` files.

## Notes for code pertaining to **A Study of the Numerical Stability of an ImEx Scheme with Application to the Poisson-Nernst-Planck Equations**

To create a collection of steady state solutions, one for each value
of epsilon: run the script ```find_steady_state.m```  This creates a file
```steady_state_data.mat``` It includes an array EPSILON and various arrays
ending in '\_ss'.  These are the steady states, for each value of
epsilon in EPSILON there's an accompanying steady state solution.

The current script creates EPSILON of length 51 and has N=91 nodes in
the interval [0,1]; the steady states are of size 91x51 and to see a
sample steady state you can ```plot(x,phi_ss(:,12))```.  To see the
electric field at the right endpoint and how it depends on epsilon you
can ```plot(EPSILON,phi_x_right_ss)```

For a fixed value of epsilon, we linearize the time-stepping scheme
about the steady-state solution.

The linearized scheme has one free parameter: the time-step size dt.
For each fixed dt, the matrix has 4N eigenvalues where N is the number
of nodes in the interval.  These eigenvalues can be real or complex
and so for each value of dt, we plot the magnitudes of the
eigenvalues.  If these are plotted for a range of dt, we get curves
showing us how the magnitudes of particular eigenvalues depend on dt.
This is one figure per epsilon.  If we then view these figures as
epsilon increases, we see which eigenvalues are the ones to determine
the numerical stability time-step constraint and can understand why
the stability diagram can have jumpes and corners in it.

To make such a movie, first create the data using
```generate_movie_data.m```.  You can then watch the movie using
```watch_movie.m``` or you can create an avi file of the movie using
```make_movie.m``` Creating the data takes over 12 hours and so the data
file ```movie_data.mat``` can be downloaded [here](https://drive.google.com/file/d/1cNDTi17iHKDwnTBmtA2n1ctaQ8pBBKS9/view?usp=sharing), sparing you of needing to use
```generate_movie_data.m```.

For each value of epsilon, there's a critical value of dt at which one
of the eigenvalues has its magnitude equal to 1, with its magnitude
being greater than 1 for larger values of dt.  This is the linear
stability constraint, dt^\*, for that value of epsilon.  To find dt^\*
for each epsilon value in EPSILON, run ```sweep.m```.  ```sweep.m``` loads
steady_state_data.mat and computes dt_thresh and then saves the data
back to steady_state_data.mat

The data file provided, steady_state_data.mat has all the steady states
and the dt_threshold already computed.  Only 51 values of epsilon are
considered in the file --- in the article we used a larger range of
epsilon values and made sure to refine the epsilon values near corners
and cusps.
