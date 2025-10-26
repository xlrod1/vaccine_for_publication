## Vaccine Bounds Project

This repository contains the code used for the partial identification of vaccine effects (VE) under broken blinding scenarios. It includes functions for deriving non parametric bounds for VE using both linear programming (LP) and non linear optimization approaches, along with numerical illustrations and a real data example.

### Folder Structure

- **LP_based_bounds_formulas**  
  Contains the derivations and functions for computing VE bounds using linear programming. These bounds follow the approach described in the manuscript.  
  Note: this folder includes the required functions adapted from the causalOptim package. The full package does not need to be installed in order to run these scripts.

- **Non_linear_optimazation**  
  Contains functions that compute VE bounds using non linear optimization methods. These functions allow for flexible incorporation of user defined constraints such as monotonicity or restricted differences.

- **Numerical_example**  
  Includes reproducible code for the numerical examples presented in Section 8 of the manuscript.

- **Real_data_example**  
  Includes the code used to analyse the ENSEMBLE2 study. A mock dataset is included for demonstration. To run the full real data analysis, you will need access to the original dataset.
