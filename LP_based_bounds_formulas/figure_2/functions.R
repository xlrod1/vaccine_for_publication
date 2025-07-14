library(rcdd)
linear_term<-function (number, string) 
{
  if (number == 0) 
    return("")
  if (number == 1) 
    return(paste0(" + ", string))
  if (number == -1) 
    return(paste0(" - ", string))
  if (number > 0) 
    return(paste0(" + ", number, string))
  return(paste0(" - ", abs(number), string))
}

linear_expression<-function (numbers, strings) 
{
  paste0(mapply(linear_term, numbers, strings), collapse = "")
}

constant_term<-function (numbers1, numbers2) 
{
  as.character(sum(numbers1 * numbers2))
}
evaluate_objective<-function (c1_num, p, y) 
{
  m1 <- nrow(c1_num)
  m <- length(y)
  number_indices <- 1:m1
  parameter_indices <- (m1 + 1):m
  y1 <- y[number_indices]
  y2 <- y[parameter_indices]
  const_term <- constant_term(numbers1 = c1_num, numbers2 = y1)
  lin_expr <- linear_expression(numbers = y2, strings = p)
  aff_expr <- paste0(as.character(const_term), lin_expr)
  if (const_term == "0" && lin_expr != "") {
    sign_of_first_lin_term <- substr(lin_expr, start = 2, 
                                     stop = 2)
    aff_expr <- substr(aff_expr, start = 5, stop = nchar(aff_expr))
    if (sign_of_first_lin_term == "-") 
      aff_expr <- paste0("-", aff_expr)
  }
  aff_expr <- paste0("  ", aff_expr)
  return(aff_expr)
}
compute_bounds <- function(target_param, B, opt = "max", constraint = NULL) {
  # Determine 'a' and 'm' based on the target parameter
  if (target_param == "Y_11") {
    a <- 1; m <- 1
    c0 <- matrix(c(0,0,0,0,1,0,0,1,0,1,1,0,1,1,1,1), byrow = TRUE)
  } else if (target_param == "Y_10") {
    a <- 1; m <- 0
    c0 <- matrix(c(0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,1), byrow= TRUE)
  } else if (target_param == "Y_01") {
    a <- 0; m <- 1
    c0 <- matrix(c(0,0,1,0,0,1,0,0,1,1,0,1,1,0,1,1), byrow = TRUE)
  } else if (target_param == "Y_00") {
    a <- 0; m <- 0
    c0 <- matrix(c(0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1), byrow = TRUE)
  } else {
    stop("Invalid target parameter. Choose from 'Y_11', 'Y_10', 'Y_01', 'Y_00'.")
  }
  
  n <- nrow(c0)
  
  # Define constraints (c3 and c4) based on B
  if (B == 0) {
    if (a == 0) {
      c4 <- c(0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1)  # (0,-1) column
    } else {
      c4 <- c(0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,1)  # (1,-1) column
    }
  } else if (B == 1) {
    if (a == 0) {
      c4 <- c(0,0,1,0,0,1,0,0,1,1,0,1,1,0,1,1)  # (0,-1) column from second table
    } else {
      c4 <- c(0,0,0,0,1,0,0,1,0,1,1,0,1,1,1,1)  # (1,-1) column from second table
    }
  } else {
    stop("Invalid B value. Choose 0 or 1.")
  }
  
  c3 <- 1 - c4  # Ensure sum is 1 in each column
  
  # Define the probability vector p
  p<- c(paste0("p0", B, ".", a, "_"), paste0("p1", B, ".", a, "_"))
  
  
  
  
  A_e <- matrix(c(c(rep(1,16)),c3,c4),nrow = 3,ncol=16,byrow = TRUE)#matrix of constarins
  #p <- c("p00.1_" ,"p10.1_")
  m_e <- nrow(A_e)
  A_l <- NULL# extra constraint , for example:monotonicity
  b_l <- NULL
  # Apply constraints to q (some q's must be zero)
  if (!is.null(constraint)) {
    indices_to_zero <- NULL
    
   if (constraint == "positive_a") {
    # Non-negative effect of A: Y^{a=0,m} ≤ Y^{a=1,m} for all m
    # Set q_{1.ab}, q_{2.ab}, q_{5.ab}, q_{7.ab}, q_{8.ab}, q_{11.ab}, q_{12.ab} = 0
    indices_to_zero <- c(1, 2, 5, 7, 8, 11, 12)
  } else if (constraint == "negative_a") {
    # Non-positive effect of A: Y^{a=1,m} ≤ Y^{a=0,m} for all m
    # Set q_{3.ab}, q_{4.ab}, q_{7.ab}, q_{8.ab}, q_{10.ab}, q_{13.ab}, q_{14.ab} = 0
    indices_to_zero <- c(3, 4, 7, 8, 10, 13, 14)
  } else if (constraint == "positive_m") {
    # Non-negative effect of m: Y^{a,m=0} ≤ Y^{a,m=1} for all a
    # Set q_{1.ab}, q_{3.ab}, q_{6.ab}, q_{7.ab}, q_{8.ab}, q_{11.ab}, q_{13.ab} = 0
    indices_to_zero <- c(1, 3, 6, 7, 8, 11, 13)
  } else if (constraint == "negative_m") {
    # Non-positive effect of m: Y^{a,m=1} ≤ Y^{a,m=0} for all a
    # Set q_{2.ab}, q_{4.ab}, q_{7.ab}, q_{8.ab}, q_{9.ab}, q_{12.ab}, q_{14.ab} = 0
    indices_to_zero <- c(2, 4, 7, 8, 9, 12, 14)
  }
    
    if (!is.null(indices_to_zero)) {
      A_l <- rbind(A_l, diag(n)[indices_to_zero, , drop = FALSE])
      data<-c(b_l, rep(0, length(indices_to_zero)))
      b_l <- matrix(data=data, nrow = length(data), ncol = 1)
    }
  }
  
  
  if (is.null(A_l)) 
    A_l <- matrix(data = 0, nrow = 0, ncol = n)
  #b_l <- NULL
  if (is.null(b_l)) 
    b_l <- matrix(data = 0, nrow = 0, ncol = 1)
  if (!is.numeric(A_l) || !is.numeric(b_l)) 
    stop("Inequality entries must be numeric.")
  if (ncol(A_l) != n) 
    stop("Dimension mismatch of inequality constaints.")
  m_l <- nrow(A_l)
  if (nrow(b_l) != m_l) 
    stop("Dimension mismatch of inequality constaints.")
  m <- m_l + m_e
  a1 <- rbind(cbind(t(A_l), t(A_e)), cbind(diag(x = 1, nrow = m_l, 
                                                ncol = m_l), matrix(data = 0, nrow = m_l, ncol = m_e)))
  b1 <- rbind(c0, matrix(data = 0, nrow = m_l, ncol = 1))
  
  
  
  #opt="max"
  a1 <- rbind(cbind(t(A_l), t(A_e)), cbind(diag(x = 1, nrow = m_l, 
                                                ncol = m_l), matrix(data = 0, nrow = m_l, ncol = m_e)))
  b1 <- rbind(c0, matrix(data = 0, nrow = m_l, ncol = 1))
  
  if (opt == "max") {
    a1 <- -a1
    b1 <- -b1
  }
  
  
  hrep <- makeH(a1 = a1, b1 = b1)
  vrep <- scdd(input = hrep, adjacency = TRUE, inputadjacency = TRUE, 
               incidence = TRUE, inputincidence = TRUE)
  matrix_of_vrep <- vrep$output
  indices_of_vertices <- matrix_of_vrep[, 1] == 0 & matrix_of_vrep[, 
                                                                   2] == 1
  vertices <- matrix_of_vrep[indices_of_vertices, -c(1, 2), 
                             drop = FALSE]
  c1_num <- rbind(b_l, 1)
  expressions <- apply(vertices, 1, function(y) evaluate_objective(c1_num = c1_num, 
                                                                   p = p, y = y))
  elements <- paste(expressions, sep = ",", collapse = ",\n")
  opt_bound <- paste0(if (opt == "min") 
    "\nMAX {\n"
    else "\nMIN {\n", elements, "\n}\n")
  opt_bound <- structure(list(expr = opt_bound, type = if (opt == 
                                                           "min") "lower" else "upper", dual_vertices = vertices, 
                              dual_vrep = vrep), class = "optbound")
  #opt_bound$expr
  #opt_bound$type
  return(list(expr = opt_bound$expr, type = opt_bound$type))
}
