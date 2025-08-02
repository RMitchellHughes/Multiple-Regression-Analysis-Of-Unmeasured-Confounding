cnfd.int.raw <- function(y, x, w, Uy, Ux) {
  # Description of Variables
  # - y:  vector (or nx1 matrix) of response variable
  # - x:  vector (or nx1 matrix)  of explanatory variable
  # - w:  data frame or matrix of covariate variables
  # - Uy: number in range (0,1) that bounds the coefficient of determination R^2_{w,u;y}
  # - Ux: number in range (0,1) that bounds the coefficient of determination R^2_{w,u;x}
  
  # Verify Parameters
  # y:
  if (!is.vector(y) && !(is.matrix(y) && ncol(y) == 1)) {
    stop("Parameter 'y' must be a vector or an nx1 matrix.")
  }
  if (any(is.na(y))) {
    stop("Parameter 'y' contains NA/NaN values.")
  }
  if (!is.numeric(y)) {
    stop("Parameter 'y' must be numeric.")
  }
  
  # x:
  if (!is.vector(x) && !(is.matrix(x) && ncol(x) == 1)) {
    stop("Parameter 'x' must be a vector or an nx1 matrix.")
  }
  if (any(is.na(x))) {
    stop("Parameter 'x' contains NA/NaN values.")
  }
  if (!is.numeric(x)) {
    stop("Parameter 'x' must be numeric.")
  }
  
  # w:
  if (!is.data.frame(w) && !is.matrix(w)) {
    stop("Parameter 'w' must be a data frame or a matrix.")
  }
  if (any(is.na(w))) {
    stop("Parameter 'w' contains NA/NaN values.")
  }
  if (!is.numeric(as.matrix(w))) {
    stop("Parameter 'w' must contain only numeric values.")
  }
  
  # Uy:
  if (!is.numeric(Uy) || length(Uy) != 1 || Uy <= 0 || Uy >= 1) {
    stop("Parameter 'Uy' must be a single numeric value between 0 and 1 (exclusive).")
  }
  if (is.na(Uy)) {
    stop("Parameter 'Uy' cannot be NA.")
  }
  
  # Ux:
  if (!is.numeric(Ux) || length(Ux) != 1 || Ux <= 0 || Ux >= 1) {
    stop("Parameter 'Ux' must be a single numeric value between 0 and 1 (exclusive).")
  }
  if (is.na(Ux)) {
    stop("Parameter 'Ux' cannot be NA.")
  }
  
  # Preprocess function parameters
  y <- as.vector(y)
  x <- as.vector(x)
  w <- as.data.frame(w)
  
  # Consistency in n
  if (length(y) != length(x) || length(y) != nrow(w)) {
    stop("The number of elements/rows in 'y', 'x', and 'w' must be consistent.")
  }
  if (length(y) < 2) {
    stop("The number of elements/rows in 'y', 'x', and 'w' must be greater than 1.")
  }
  
  # Function to calculate confounding interval
  # From "Regression Analysis of Unmeasured Confounding" by Knaeble, Osting, Abramson
  # Function was written by Abramson and supplied by Knaeble. See link below:
  # https://github.com/bknaeble/ConfoundingIntervals/tree/master
  f=function(p,sr,lx2,ux2,ly2,uy2,lxy,uxy) {
    tol=.00001
    lx=sqrt(lx2)
    ux=sqrt(ux2)
    ly=sqrt(ly2)
    uy=sqrt(uy2)
    v=numeric(88*3)
    M=matrix(v,ncol=3)
    bx=c(lx,ux) 
    by=c(ly,uy)
    bxy=c(lxy,uxy)
    s1=function(bx,by,bxy) c(((-2*p+sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
    s2=function(bx,by,bxy) c(((-2*p-sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
    s3=function(bx,by,bxy) c((p+1)/(bxy+1),(p+1)/(bxy+1),bxy)
    s4=function(bx,by,bxy) c((p-1)/(bxy-1),(p-1)/(bxy-1),bxy)
    s5=function(bx,by,bxy) c(bx^2,by^2,bxy)
    s6=function(bx,by,bxy) c(bx,by,(p+sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
    s7=function(bx,by,bxy) c(bx,by,(p-sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
    s8=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
    s9=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
    s10=function(bx,by,bxy) c(((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
    s11=function(bx,by,bxy) c(((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
    w=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11)
    for (i in 1:2) {
      for (j in 1:2) {
        for (k in 1:2) {
          for (l in 1:11) {
            r=((l-1)*8)+(4*(i-1)+2*(j-1)+1*(k-1)+1)
            M[r,]=w[[l]](bx[i],by[j],bxy[k])
          }}}}
    M
    feas=function(arg) 
    {
      arg[1]>=(lx2-tol) & 
        arg[1]<=(ux2+tol) &
        arg[2]>=(ly2-tol) &
        arg[2]<=(uy2+tol) &
        arg[3]>=(lxy-tol) &
        arg[3]<=(uxy+tol) 
    }
    fs=apply(M,1,feas)
    N=M[which(fs==TRUE),]
    obj=function(arg) sr*(p-sqrt(arg[1])*sqrt(arg[2])*arg[3])/(1-arg[1])
    l=min(apply(N,1,obj), na.rm = TRUE)
    u=max(apply(N,1,obj), na.rm = TRUE)
    return(c(l,u))
  }
  
  # This is where the methodology of "Partial Identification For Causal Inference
  # With Regression Models" by Knaeble and Hughes begins
  
  # Regress y and x onto w
  model.yw <- lm(y ~ ., data = w)
  model.xw <- lm(x ~ ., data = w)
  
  # Get residuals from models
  res.yw <- model.yw$residuals
  res.xw <- model.xw$residuals
  
  # Calculate R^2 values from models
  R2.wy <- summary(model.yw)$r.squared
  R2.wx <- summary(model.xw)$r.squared
  
  # Calculate upper bounds uy and ux
  uy <- (Uy - R2.wy) / (1 - R2.wy)
  ux <- (Ux - R2.wx) / (1 - R2.wx)
  
  # Apply algorithm of "Regression Analysis of Unmeasured Confounding" to residual vectors
  confounding.interval <- suppressWarnings(f(cor(res.yw, res.xw), sd(res.yw) / sd(res.xw), 0, ux, 0, uy, -1, 1))
  
  return(confounding.interval)
}