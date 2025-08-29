cnfd.int.contour <- function(cor.res, sd.ratio.res, R2.wx, R2.wy, granularity = 50, contour.config = NA, axis.labels = c('x', 'y')) {
  # Description of Variables: For all of the following variables, res.xw is defined as the
  # residuals from regressing x onto w and res.yw is the residuals of regressing y onto w
  
  # - cor.res: number in range [-1,1] representing the correlation between res.xw and res.yw
  # - sd.ratio.res: non-negative number representing the standard deviation of res.yw divided by standard deviation of res.xw
  # - R2.wx: number in [0,1) that represents the coefficient of determination R^2_{w;x}
  # - R2.wy: number in [0,1) that represents the coefficient of determination R^2_{w;y}
  # - practical.significance: number representing a threshold of practical significance
  # - contour.config: a numeric vector of length three c(start, end, step), or NA to use default contours values, note 
  #   - start: the level for the first contour
  #   - end: the level for the last contour
  #   - step: the interval between contour levels
  # - granularity: integer in [50,300] with larger number producing more refined plots
  # - axis.labels: vector of strings with the names of the x and y variables in that order
  
  # Verify plotly installed
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "The 'plotly' package is required to run this function.\n",
      "Please install it using: install.packages('plotly')"
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
  
  # contour.config:
  # Must be a numeric vector of length 3, or a single NA value.
  if (! (any(is.na(contour.config)) && length(contour.config) == 1) &&
      ! (is.numeric(contour.config) && length(contour.config) == 3) ) {
    stop("Parameter 'contour.config' must be a numeric vector of length three or a single NA value.")
  }

  # If numeric, must not contain NA/NaN values.
  if (is.numeric(contour.config) && any(is.na(contour.config))) {
    stop("Parameter 'contour.config' cannot contain NA/NaN values when specified as a vector.")
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
  Bx_vals <- seq(R2.wx, .99, length.out = granularity)
  By_vals <- seq(R2.wy, .99, length.out = granularity)
  
  # Calculate upper and lower bounds U(Bx,By) and L(Bx,By) on Bx,By grid
  U.BxBy <- matrix(ncol = granularity, nrow = granularity)
  L.BxBy <- matrix(ncol = granularity, nrow = granularity)
  for (i in 1:granularity) {
    upper_limit_x <- (Bx_vals[i] - R2.wx) / (1 - R2.wx)
    
    for (j in 1:granularity) {
      upper_limit_y <- (By_vals[j] - R2.wy) / (1 - R2.wy)
      result <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x, 0, upper_limit_y, -1, 1))
      U.BxBy[j,i] <- result[2]
      L.BxBy[j,i] <- result[1]
    }
  }
  
  # Remove Inf values
  U.BxBy[U.BxBy == Inf] <- NA
  U.BxBy[U.BxBy == -Inf] <- NA
  L.BxBy[L.BxBy == Inf] <- NA
  L.BxBy[L.BxBy == -Inf] <- NA
  
  # Get Initial Slope
  initial.slope <- suppressWarnings(f(cor.res, sd.ratio.res, 0, 0, 0, 0, -1, 1)[1])
  
  # Define tooltiop
  tooltip <- matrix(
    paste0(
      "B<sub>", axis.labels[1], "</sub> = ", round(Bx_vals[col(L.BxBy)], 5), "<br>",
      "B<sub>", axis.labels[2], "</sub> = ", round(By_vals[row(L.BxBy)], 5), "<br>",
      "U(B<sub>", axis.labels[1], "</sub>,B<sub>", axis.labels[2], "</sub>) = ", round(U.BxBy, 5), "<br>",
      "L(B<sub>", axis.labels[1], "</sub>,B<sub>", axis.labels[2], "</sub>) = ", round(L.BxBy, 5)
    ),
    nrow = nrow(L.BxBy),
    ncol = ncol(L.BxBy)
  )
  
  # Adjust max/min if needed so contour lines appear correctly
  if (!any(is.na(contour.config))) {
    contour.config <- c(min(contour.config[1:2]), max(contour.config[1:2]), contour.config[3])
    if (initial.slope > 0) {
      contour.config[1] <- contour.config[1] + (contour.config[2] - contour.config[1]) %% contour.config[3]
      }
    } else {
      contour.config[2] <- contour.config[2] - (contour.config[2] - contour.config[1]) %% contour.config[3]
    }
  }
  
  # Settings for contour lines
  if (any(is.na(contour.config))) {
    contour.settings <- list(
      type = 'levels',
      showlabels = TRUE
    )
  } else {
    contour.settings <- list(
      type = 'levels',
      start = contour.config[1],
      end = contour.config[2],
      size = contour.config[3],
      showlabels = TRUE
    )
  }
  
  # Create the Plot
  fig <- plot_ly()
  
  if (initial.slope > 0) {
    # Use L(Bx,By) if initial.slope is positive
    fig <- fig %>% add_contour(
      x = Bx_vals,
      y = By_vals,
      z = L.BxBy,
      type = "contour",
      autocontour = is.na(contour.config),
      contours = contour.settings,
      colorscale = list(c(0, "rgb(50,50,50)"), c(1, "rgb(200,200,200)")),
      line = list(color = "black"),
      text = tooltip,
      hoverinfo = "text"
    )
    
    line.color <- 'black'
  } else {
    # Use U(Bx,By) if initial.slope is non-positive
    fig <- fig %>% add_contour(
      x = Bx_vals,
      y = By_vals,
      z = U.BxBy,
      type = "contour",
      autocontour = is.na(contour.config),
      contours = contour.settings,
      colorscale = list(c(0, "rgb(50,50,50)"), c(1, "rgb(200,200,200)")),
      line = list(color = "white"),
      text = tooltip,
      hoverinfo = "text"
    )
    
    line.color <- 'gray'
  }
  
  # Add line for R2.wx
  fig <- fig %>% add_trace(
    x = rep(R2.wx,200),
    y = seq(0,1,length.out=200),
    type = "scatter",
    mode = "lines",
    line = list(color = line.color, width = 5),
    showlegend = FALSE,
    name = paste0('R<sup>2</sup><sub>wx</sub> = ', round(R2.wx,3))
  )
  
  # Add line for R2.wy
  fig <- fig %>% add_trace(
    x = seq(0,1,length.out=200),
    y = rep(R2.wy,200),
    type = "scatter",
    mode = "lines",
    line = list(color = line.color, width = 5),
    showlegend = FALSE,
    name = paste0('R<sup>2</sup><sub>wy</sub> = ', R2.wy)
  )
  
  # Add title, axis labels,
  fig <- fig %>% config(
    displayModeBar = FALSE
  )
  
  fig <- fig %>% layout(
    title = list(text = plotly::TeX(paste0("\\text{Contour Plot of } L(B_{\\text{", axis.labels[1], "}}, B_{\\text{", axis.labels[2], "}})"))),
    xaxis = list(title = plotly::TeX(paste0("B_{\\text{", axis.labels[1], "}}")),
                 range = c(0,1)),
    yaxis = list(title = plotly::TeX(paste0("B_{\\text{", axis.labels[2], "}}")),
                 range = c(0,1)),
    margin = list(t = 70, l = 70, r = 70, b = 70)
  )
  
  fig <- fig %>% config(
    mathjax = "cdn"
  )
  
  return(fig)
}
