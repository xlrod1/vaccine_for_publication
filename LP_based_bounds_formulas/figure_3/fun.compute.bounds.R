
comute.bounds<-function(c0,c1,c2,p,nrow,ncol,opt){
  c0<-c0
  A_e <- matrix(c(c(rep(1,ncol)),c1,c2),nrow = nrow,,ncol=ncol,byrow = TRUE)#matrix of constarins
  p <- p
  m_e <- nrow(A_e)
  A_l <- NULL# extra constraint , for example:monotonicity
  n <- nrow(c0)
  if (is.null(A_l)) 
    A_l <- matrix(data = 0, nrow = 0, ncol = n)
  b_l <- NULL
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
  
  
  
  #opt="min"
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
  opt_bound$expr
  #opt_bound$type
}

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