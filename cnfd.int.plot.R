cnfd.int.plot <- function(cor.res, sd.ratio.res, R2.wy, R2.wx, granularity = 50, axis.labels = c('y', 'x')) {
  # Description of Variables: For all of the following variables, res.yw is defined as the
  # residuals from regressing y onto w and res.xw is the residuals of regressing x onto w
  
  # - cor.res: number in range [-1,1] representing the correlation between res.yw  and res.xw
  # - sd.ratio.res: non-negative number representing the standard deviation of res.yw divided by standard deviation of res.xw
  # - R2.wy: number in [0,1) that represents the coefficient of determination R^2_{w;y}
  # - R2.wx: number in [0,1) that represents the coefficient of determination R^2_{w;x}
  # - granularity: integer in [50,300] with larger number producing more refined plots
  
  
  # Verify plotly installed
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "The 'plotly' package is required to run this function.\n",
      "Please install it using: install.packages('plotly')"
    )
  }
  
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop(
      "The 'latex2exp' package is required to run this function.\n",
      "Please install it using: install.packages('latex2exp')"
    )
  }
  
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
  
  # granularity:
  # Must be a single integer
  if (!(round(granularity) == granularity) || length(granularity) != 1) {
    stop("Parameter 'granularity' must be a single integer value.")
  }
  # Must not be NA/NaN.
  if (is.na(granularity)) {
    stop("Parameter 'granularity' cannot be NA/NaN.")
  }
  # Must be in range [50,300]
  if (granularity < 50 || granularity > 300) {
    stop("Parameter 'granularity' must be in range [50,300].")
  }
  
  # axis.labels:
  # Must be a character vector.
  if (!is.character(axis.labels)) {
    stop("Parameter 'axis.labels' must be a vector of strings.")
  }
  # Must be length 2.
  if (length(axis.labels) != 2) {
    stop("Parameter 'axis.labels' must be a vector of length 2.")
  }
  # Must not contain NA values.
  if (any(is.na(axis.labels))) {
    stop("Parameter 'axis.labels' cannot contain NA values.")
  }
  # Must not contain empty strings.
  if (any(axis.labels == "")) {
    stop("Parameter 'axis.labels' cannot contain empty strings.")
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
  
  # Fine grid for values of Bx,By
  # Bx_vals <- seq(R2.wx, 1, length.out = granularity+1)[1:granularity]
  # By_vals <- seq(R2.wy, 1, length.out = granularity+1)[1:granularity]
  Bx_vals <- seq(R2.wx, .96, length.out = granularity)
  By_vals <- seq(R2.wy, .96, length.out = granularity)
  
  # Calculate upper and lower bounds U(Bx,By) and L(Bx,By) on Bx,By grid
  U.BxBy <- matrix(ncol = granularity, nrow = granularity)
  L.BxBy <- matrix(ncol = granularity, nrow = granularity)
  for (i in 1:granularity) {
    for (j in 1:granularity) {
      # Calculate Upper Limits From Equation 12
      upper_limit_y <- (By_vals[j] - R2.wy) / (1 - R2.wy)
      upper_limit_x <- (Bx_vals[i] - R2.wx) / (1 - R2.wx)
      
      U.BxBy[j,i] <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x, 0, upper_limit_y, -1, 1)[2])
      L.BxBy[j,i] <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x, 0, upper_limit_y, -1, 1)[1])
    }
  }
  
  # Remove Inf values
  U.BxBy[U.BxBy == Inf] <- NA
  U.BxBy[U.BxBy == -Inf] <- NA
  L.BxBy[L.BxBy == Inf] <- NA
  L.BxBy[L.BxBy == -Inf] <- NA
  
  # Define 3d plot "fig"
  fig <- plot_ly()
  
  # Define tooltips
  upper.tooltip <- matrix(
    paste0(
      "B<sub>", axis.labels[2], "</sub>: ", round(Bx_vals[col(U.BxBy)], 5), "<br>",
      "B<sub>", axis.labels[1], "</sub>: ", round(By_vals[row(U.BxBy)], 5), "<br>",
      "U(B<sub>", axis.labels[2], "</sub>,B<sub>", axis.labels[1], "</sub>) = ", round(U.BxBy, 5)
    ),
    nrow = nrow(U.BxBy),
    ncol = ncol(U.BxBy)
  )
  
  lower.tooltip <- matrix(
    paste0(
      "B<sub>", axis.labels[2], "</sub>: ", round(Bx_vals[col(L.BxBy)], 5), "<br>",
      "B<sub>", axis.labels[1], "</sub>: ", round(By_vals[row(L.BxBy)], 5), "<br>",
      "L(B<sub>", axis.labels[2], "</sub>,B<sub>", axis.labels[1], "</sub>) = ", round(L.BxBy, 5)
    ),
    nrow = nrow(L.BxBy),
    ncol = ncol(L.BxBy)
  )
  
  # Add upper bound U(Bx,By)
  # U(Bx,By) >= 0 ==> light blue
  # U(Bx,By) <  0 ==> dark blue
  if (all(U.BxBy >= 0, na.rm = TRUE)) {
    fig <- fig %>% add_surface(
      x = Bx_vals,
      y = By_vals,
      z = U.BxBy,
      colorscale = list(c(0, 1), c("skyblue", "skyblue")),
      name = "U(Bx,By)",
      showscale = FALSE,
      showlegend = TRUE,
      text = upper.tooltip,
      hoverinfo = "text"
    )
  } else {
    U.BxBy_color <- U.BxBy
    U.BxBy_color[U.BxBy >= 0] <- 1
    U.BxBy_color[U.BxBy < 0] <- 0
    
    fig <- fig %>% add_surface(
      x = Bx_vals,
      y = By_vals,
      z = U.BxBy,
      surfacecolor = U.BxBy_color,
      colorscale = list(c(0, 0.5, 0.5, 1), c("darkblue", "darkblue", "skyblue", "skyblue")),
      name = "U(Bx,By)",
      showscale = FALSE,
      showlegend = TRUE,
      text = upper.tooltip,
      hoverinfo = "text"
    )
  }
  
  # Add lower bound L(Bx,By)
  # L(Bx,By) >= 0 ==> red
  # L(Bx,By) <  0 ==> dark red
  if (all(L.BxBy < 0, na.rm = TRUE)) {
    fig <- fig %>% add_surface(
      x = Bx_vals,
      y = By_vals,
      z = L.BxBy,
      colorscale = list(c(0, 1), c("darkred", "darkred")),
      name = "L(Bx,By)",
      showscale = FALSE,
      showlegend = FALSE,
      text = lower.tooltip,
      hoverinfo = "text"
    )
  } else {
    L.BxBy_color <- L.BxBy
    L.BxBy_color[L.BxBy >= 0] <- 1
    L.BxBy_color[L.BxBy < 0] <- 0
    
    fig <- fig %>% add_surface(
      x = Bx_vals,
      y = By_vals,
      z = L.BxBy,
      surfacecolor = L.BxBy_color,
      colorscale = list(c(0, 0.5, 0.5, 1), c("darkred", "darkred", "red", "red")),
      name = "L(Bx,By)",
      showscale = FALSE,
      showlegend = FALSE,
      text = lower.tooltip,
      hoverinfo = "text"
    )
  }
  
  # Add title, axis labels, setup camera, etc.
  fig <- fig %>% config(
    displayModeBar = FALSE
  )
  
  fig <- fig %>% config(
    mathjax = "cdn"
  )
  
  fig <- fig %>% layout(
    scene = list(
      camera = list(
        eye = list(x = -1.25, y = -1.25, z = 1.25),
        center = list(x = 0, y = 0, z = 0),
        up = list(x = 0, y = 0, z = 1)
      ),
      xaxis = list(title = paste0("B<sub>", axis.labels[2], "</sub>")),
      yaxis = list(title = paste0("B<sub>", axis.labels[1], "</sub>")),
      # zaxis = list(title = plotly::TeX(paste0("\\beta_{B_{", axis.labels[2], "}. B_{", axis.labels[1], "}|w,u}")))
      zaxis = list(title = paste0("\u03b2<sub>", axis.labels[2], ".", axis.labels[1], "|w,u</sub>"))
    ),
    # title = "Confounding Interval of beta_{log(SD).log(BMI)|w,u}"
    title = list(text = plotly::TeX(paste0("\\text{Confounding Interval of } \\beta_{\\text{", axis.labels[2], "} . \\text{", axis.labels[1], "}|w,u}")))
  )
  
  return(fig)
}
