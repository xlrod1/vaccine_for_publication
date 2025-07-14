compute.bounds.DAG_YS<-function(B,A,S,ATE){
  
  #crating the Y table
  
  y.ams<-expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
  names(y.ams)<-c("000","001","010","011","100","101","110","111")
  
  
  y.ams.b0<-y.ams
  y.ams.b1<-y.ams
  
  y.ams.b0$'0-10'<-y.ams$'000'
  y.ams.b0$'0-11'<-y.ams$'001'
  y.ams.b0$'1-10'<-y.ams$'100'
  y.ams.b0$'1-11'<-y.ams$'101'
  
  
  y.ams.b1$'0-10'<-y.ams$'010'
  y.ams.b1$'0-11'<-y.ams$'011'
  y.ams.b1$'1-10'<-y.ams$'110'
  y.ams.b1$'1-11'<-y.ams$'111'  
  
  cs<-paste(A,-1,S,sep = "")
  if (B==0){
    c1<-y.ams.b0[,cs]
    c2<-1-c1
    p1<-paste("p","1",A,B,S,sep = "")
    p2<-paste("p","0",A,B,S,sep = "")
    p <- c(p1,p2)
  }
  if (B==1){
    c1<-y.ams.b1[,cs]
    c2<-1-c1
    p1<-paste("p","1",A,B,S,sep = "")
    p2<-paste("p","0",A,B,S,sep = "")
    p <- c(p1,p2)
  }
  
  
  #The average treatment effcet
  c0<-matrix(c(y.ams[,ATE]),byrow = TRUE)
  
  nrow = length(p)+1
  ncol<-nrow(y.ams)
  
  
  
  
  min<-comute.bounds(c0,c1,c2,p,nrow,ncol,"min")
  max<-comute.bounds(c0,c1,c2,p,nrow,ncol,"max")
  
  return(list(min=min,max=max))
}
