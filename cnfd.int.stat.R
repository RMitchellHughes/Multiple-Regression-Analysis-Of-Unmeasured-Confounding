cnfd.int.stat <- function(cor.res, sd.ratio.res, R2.wy, R2.wx, Uy, Ux) {
  # Description of Variables: For all of the following variables, res.yw is defined as the
  # residuals from regressing y onto w and res.xw is the residuals of regressing x onto w
  
  # - cor.res: number in range [-1,1] representing the correlation between res.yw  and res.xw
  # - sd.ratio.res: non-negative number representing the standard deviation of res.yw divided by standard deviation of res.xw
  # - R2.wy: number in [0,1) that represents the coefficient of determination R^2_{w;y}
  # - R2.wx: number in [0,1) that represents the coefficient of determination R^2_{w;x}
  # - Uy: number in range (0,1) that bounds the coefficient of determination R^2_{w,u;y}
  # - Ux: number in range (0,1) that bounds the coefficient of determination R^2_{w,u;x}
  
  # Verify Parameters
  # cor.res:
  # Must be a single numeric value.
  if (!is.numeric(cor.res) || length(cor.res) != 1) {
    stop("Parameter 'cor.res' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(cor.res)) {
    stop("Parameter 'cor.res' cannot be NA/NaN.")
  }
  # Must be within the range [-1, 1].
  if (cor.res < -1 || cor.res > 1) {
    stop("Parameter 'cor.res' must be in the range [0, 1].")
  }
  
  # sd.ratio.res:
  # Must be a single numeric value.
  if (!is.numeric(sd.ratio.res) || length(sd.ratio.res) != 1) {
    stop("Parameter 'sd.ratio.res' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(sd.ratio.res)) {
    stop("Parameter 'sd.ratio.res' cannot be NA/NaN.")
  }
  # Must bve non-negative
  if (sd.ratio.res < 0) {
    stop("Parameter 'sd.ratio.res' must be non-negative.")
  }
  
  # R2.wy:
  # Must be a single numeric value.
  if (!is.numeric(R2.wy) || length(R2.wy) != 1) {
    stop("Parameter 'R2.wy' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(R2.wy)) {
    stop("Parameter 'R2.wy' cannot be NA/NaN.")
  }
  # Must be in range [0,1)
  if (R2.wy < 0 || R2.wy >= 1) {
    stop("Parameter 'R2.wy' must be a single numeric value in range [0,1).")
  }
  
  # R2.wx:
  # Must be a single numeric value.
  if (!is.numeric(R2.wx) || length(R2.wx) != 1) {
    stop("Parameter 'R2.wx' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(R2.wx)) {
    stop("Parameter 'R2.wx' cannot be NA/NaN.")
  }
  # Must be in range [0,1)
  if (R2.wx < 0 || R2.wx >= 1) {
    stop("Parameter 'R2.wx' must be a single numeric value in range [0,1).")
  }
  
  # Uy:
  # Must be a single numeric value.
  if (!is.numeric(Uy) || length(Uy) != 1) {
    stop("Parameter 'Uy' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(Uy)) {
    stop("Parameter 'Uy' cannot be NA/NaN.")
  }
  # Must be strictly between 0 and 1 (exclusive).
  if (Uy <= 0 || Uy >= 1) {
    stop("Parameter 'Uy' must be a single numeric value between 0 and 1 (exclusive).")
  }
  
  # Ux:
  # Must be a single numeric value.
  if (!is.numeric(Ux) || length(Ux) != 1) {
    stop("Parameter 'Ux' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(Ux)) {
    stop("Parameter 'Ux' cannot be NA/NaN.")
  }
  # Must be strictly between 0 and 1 (exclusive).
  if (Ux <= 0 || Ux >= 1) {
    stop("Parameter 'Ux' must be a single numeric value between 0 and 1 (exclusive).")
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
  
  # Calculate upper bounds uy and ux
  uy <- (Uy - R2.wy) / (1 - R2.wy)
  ux <- (Ux - R2.wx) / (1 - R2.wx)
  
  # Apply algorithm of "Regression Analysis of Unmeasured Confounding" to residual vectors
  confounding.interval <- suppressWarnings(f(cor.res, sd.ratio.res, 0, ux, 0, uy, -1, 1))
  
  return(confounding.interval)
}