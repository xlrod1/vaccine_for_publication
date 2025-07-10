#Upload all functions
source("G:\\My Drive\\PhD3\\vaccine_for_publications\\LP_based_bounds_formulas\\figure_2\\functions.R")

# VE(0)
compute_bounds("Y_10",B =  0,"min",constraint = NULL)
compute_bounds("Y_10",B =  1,"min",constraint = NULL)
compute_bounds("Y_00",B =  0,"max",constraint = NULL)
compute_bounds("Y_00",B =  1,"max",constraint = NULL)

#VE(1)
