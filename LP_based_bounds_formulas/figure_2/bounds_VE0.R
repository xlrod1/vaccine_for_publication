#Upload all functions
source("figure 2 functions path.R")

# VE(0)
compute_bounds("Y_10",B =  0,"min",constraint = NULL)
compute_bounds("Y_10",B =  1,"min",constraint = NULL)
compute_bounds("Y_00",B =  0,"max",constraint = NULL)
compute_bounds("Y_00",B =  1,"max",constraint = NULL)

#VE(1)
