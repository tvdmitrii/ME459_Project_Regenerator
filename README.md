# Nonlinear Simultaneous Equations Solver Based on Broyden’s Method

As a research assistant at Solar Energy Lab, I worked on regenerative heat exchanger design code
that calculates the physical dimensions and performance metrics of the regenerator based on its
thermodynamic model and target parameters. Design code needs to solve a system of nonlinear
non-differentiable simultaneous equations. Equation solver was originally implemented using
nested monotonic equation solvers, which can handle only one function of one variable. System of
equations was split into five separate equations of one variable by picking single variable that
seemed to play the largest role. Those equations were then solved in a hierarchy. First, last four
variables would be fixed, and first variable would be changed until the first equation was solved.
Then, first variable was fixed, and second variable was changed once to try and get the value of a
second function closer to its desired value. After that last four variables would be fixed again, and
first variable would be adjusted until first equation was solved. The process would then slowly
propagate, starting to affect third, fourth and fifth variables, until the whole system was solved.
First two issues are speed and scalability. Nested structure of monotonic equation solvers is
inefficient, because solution is discarded 99% of the time. If a value of the fifth variable is wrong,
all the extensive work of precisely solving four equations down the hierarchy is unnecessary, yet it
is required. Moreover, monotonic equation solver requires two initial guess values to start the
iterative solution process, which means that thermodynamic model must be called twice to evaluate
all the functions. For every guess value, the solver one step down in the hierarchy needs two guess
values of its own and so on. Just this preparation work requires 2 N calls to thermodynamic model,
where N is the number of nested solvers. This has significant impact on speed and limits the
practical number of equations in a system, which was the main problem that I ran into, as
complexity of the model increased.
Second issue is robustness to initial guess values. Since at any given point only one equation is
being solved and only one variable can be adjusted, incorrect guess values of one solver make it
impossible to solve the equations down the hierarchy, because of physical limits of the variables.
For example, diameter can be set so large, that even the smallest reasonable length produces
conductance value that is too large.
Third issue is difficulty of code maintenance due to large overhead. Five solvers calling each other
in a specific order, spread out variables and functions evaluations, ten guess values – all this makes
it difficult to change order of the solvers or to introduce new ones.
To resolve problems of the original approach, algorithm using Broyden’s method was
implemented. Broyden’s method provides an iterative way for finding roots of a system of
nonlinear non-differentiable simultaneous equations. Thermodynamic model is complicated and
relies on interpolated values from lookup tables, which makes obtaining a closed-form expression
for a derivative of a function impossible. Broyden’s method uses Jacobian to iteratively achieve
successively better approximations of the roots. The functions can be of multiple variables, which
would improve robustness to initial guess values. Moreover, Jacobian is not calculated at every
iteration, but rather is successively improved, which cuts down on calls to thermodynamic model
and improves the performance of the code.
