Vaccine Bounds Project

This repository contains code for computing non parametric bounds for vaccine effects (VE) under broken blinding scenarios. The code includes linear programming (LP) based bounds, non linear optimization based bounds and monotonicity-based bounds, with options to impose additional assumptions such as M monotonicity.


The repository includes:

1.LP based bounds: Functions for computing upper and lower bounds for VE using linear programming.

2.Monotonicity based bounds: Code for incorporating monotonicity assumptions to obtain sharper bounds.

3.Non linear optimization: R functions that directly maximize or minimize the VE expressions under user specified constraints.

4.Simulation code: Scripts for generating data under different causal data generating mechanisms (DGMs) to evaluate the bounds.

