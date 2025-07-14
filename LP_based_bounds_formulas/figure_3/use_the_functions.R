#library(causaloptim)
library(rcdd)
source("G:/My Drive/PhD3/bounds/bounds.DAG2d/new.R")
source("G:/My Drive/PhD3/bounds/bounds.DAG2d/fun.compute.bounds.R")
#Y^{11}
compute.bounds.DAG_YS(A = 0,B = 0,S=0,"110")
compute.bounds.DAG_YS(A = 1,B = 0,S=0,"110")
compute.bounds.DAG_YS(A = 0,B = 0,S=1,"111")
compute.bounds.DAG_YS(A = 1,B = 0,S=1,"111")

compute.bounds.DAG_YS(A = 0,B = 1,S=0,"110")
compute.bounds.DAG_YS(B = 1,A = 1,S=0,"110")
compute.bounds.DAG_YS(A = 0,B = 1,S=1,"111")
compute.bounds.DAG_YS(B = 1,A = 1,S=1,"111")


#Y^{1,0}
compute.bounds.DAG_YS(A = 0,B = 0,S=0,"100")
compute.bounds.DAG_YS(B = 0,A = 1,S=0,"101")
compute.bounds.DAG_YS(A = 0,B = 0,S=1,"101")
compute.bounds.DAG_YS(B = 0,A = 1,S=1,"101")

compute.bounds.DAG_YS(A = 0,B = 1,S=0,"100")
compute.bounds.DAG_YS(B = 1,A = 1,S=0,"100")
compute.bounds.DAG_YS(A = 0,B = 1,S=0,"100")
compute.bounds.DAG_YS(B = 1,A = 1,S=1,"101")


#Y^{0,1}
compute.bounds.DAG_YS(A = 0,B = 0,S=0,"010")
compute.bounds.DAG_YS(B = 0,A = 1,S=0,"011")
compute.bounds.DAG_YS(A = 0,B = 0,S=1,"011")
compute.bounds.DAG_YS(B = 0,A = 1,S=1,"011")

compute.bounds.DAG_YS(A = 0,B = 1,S=0,"010")
compute.bounds.DAG_YS(B = 1,A = 1,S=0,"010")
compute.bounds.DAG_YS(A = 0,B = 1,S=1,"011")
compute.bounds.DAG_YS(B = 1,A = 1,S=1,"011")
