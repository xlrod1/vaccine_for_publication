set.seed(2025)

#Source function
folder_path <- "G:\\My Drive\\PhD3\\vaccine_for_publications\\Non_linear_optimazation\\functions"
# List all R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file
sapply(r_files, source)


# Load your dataset
data <- read.csv("G:\\My Drive\\PhD3\\vaccine_for_publications\\Non_linear_optimazation\\example_data.csv")


#Without S
#calculate the LP bounds without monotonically_constraint
result_LP_no_S<-get_all_LP_bounds_summary_no_s (data,  monotonicity_constraint ="NULL")
#calculate the LP bounds with monotonicity_constraint
result_LP_M_no_S<-get_all_LP_bounds_summary_no_s (data,  monotonicity_constraint = "positive_M")



#With S
#calculate the LP bounds without monotonically_constraint
result_LP<-get_all_LP_bounds_summary(data,  monotonicity_constraint ="NULL")
#calculate the LP bounds with monotonicity_constraint
result_LP_M<-get_all_LP_bounds_summary(data,  monotonicity_constraint = "positive_M")


