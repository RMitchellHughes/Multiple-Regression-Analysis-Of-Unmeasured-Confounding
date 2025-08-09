cnfd.int.prac <- function(cor.res, sd.ratio.res, R2.wx, R2.wy, practical.significance, granularity = 100, axis.labels = c('x', 'y')) {
  # Description of Variables: For all of the following variables, res.yw is defined as the
  # residuals from regressing y onto w and res.xw is the residuals of regressing x onto w
  
  # - cor.res: number in range [-1,1] representing the correlation between res.yw  and res.xw
  # - sd.ratio.res: non-negative number representing the standard deviation of res.yw divided by standard deviation of res.xw
  # - R2.wy: number in [0,1) that represents the coefficient of determination R^2_{w;y}
  # - R2.wx: number in [0,1) that represents the coefficient of determination R^2_{w;x}
  # - practical.significance: number that represents practical significance for beta_{x.y}, should be between beta_{x.y} and 0 for best results
  # - granularity: integer in [100,300] with larger number producing more refined plots
  # - axis.labels: a vector of strings of length 2 that represent the labels for the x and y axes (in that order)
  
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
  
  # practical.significance:
  # Must be a single numeric value.
  if (!is.numeric(practical.significance) || length(practical.significance) != 1) {
    stop("Parameter 'practical.significance' must be a single numeric value.")
  }
  # Must not be NA/NaN.
  if (is.na(practical.significance)) {
    stop("Parameter 'practical.significance' cannot be NA/NaN.")
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
  # Must be in range [100,300]
  if (granularity < 100 || granularity > 300) {
    stop("Parameter 'granularity' must be in range [100,300].")
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
  Bx_vals <- seq(R2.wx, 1, length.out = granularity+1)[1:granularity]
  By_vals <- seq(R2.wy, 1, length.out = granularity+1)[1:granularity]
  
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
  
  # Define tooltip
  tooltip <- matrix(
    paste0(
      "B<sub>", axis.labels[1], "</sub> = ", round(Bx_vals[col(U.BxBy)], 5), "<br>",
      "B<sub>", axis.labels[2], "</sub> = ", round(By_vals[row(U.BxBy)], 5), "<br>",
      "U(B<sub>", axis.labels[1], "</sub>,B<sub>", axis.labels[2], "</sub>) = ", round(U.BxBy, 5), "<br>",
      "L(B<sub>", axis.labels[1], "</sub>,B<sub>", axis.labels[2], "</sub>) = ", round(L.BxBy, 5)
    ),
    nrow = nrow(U.BxBy),
    ncol = ncol(U.BxBy)
  )
  
  # Determine direction of inequality
  initial.slope <- suppressWarnings(f(cor.res, sd.ratio.res, 0, 0, 0, 0, -1, 1)[1])
  if (practical.significance > 0 | (practical.significance == 0 & initial.slope > 0)) {
    bound.index <- 1
    plot.title <- plotly::TeX(paste0('\\large \\text{Region where } \\beta_{\\text{', axis.labels[1], '}. \\text{', axis.labels[2], '}|w,u} \\ge', practical.significance))
  } else {
    bound.index <- 2
    plot.title <- plotly::TeX(paste0('\\large \\text{Region where } \\beta_{\\text{', axis.labels[1], '}. \\text{', axis.labels[2], '}|w,u} \\le', practical.significance))
  }
  
  # Initialize Plot
  fig <- plot_ly() %>%
    layout(
      title = list(text = plot.title),
      xaxis = list(title = plotly::TeX(paste0("\\Large R^2_{w,u;\\text{", axis.labels[1], "}}")),
                   range = c(0,1)),
      yaxis = list(title = plotly::TeX(paste0("\\Large R^2_{w,u;\\text{", axis.labels[2], "}}")),
                   range = c(0,1)),
      margin = list(t = 100, l = 100, r = 100, b = 100)
    )%>%
    config(
      mathjax = "cdn",
      displayModeBar = FALSE
    )
  
  # Add line for R2.wx
  fig <- fig %>%
    add_trace(
      x = rep(R2.wx,200),
      y = seq(0,1,length.out=200),
      type = "scatter",
      mode = "lines",
      line = list(color = 'black', width = 2),
      showlegend = FALSE,
      name = paste0('R<sup>2</sup><sub>wx</sub> = ', round(R2.wx,3))
    )
  
  # Add line for R2.wy
  fig <- fig %>%
    add_trace(
      x = seq(0,1,length.out=200),
      y = rep(R2.wy,200),
      type = "scatter",
      mode = "lines",
      line = list(color = 'black', width = 2),
      showlegend = FALSE,
      name = paste0('R<sup>2</sup><sub>wy</sub> = ', round(R2.wy,3))
    )
  
  # Add tooltip
  fig <- fig %>%
    add_trace(
      type = 'heatmap',
      x = Bx_vals,
      y = By_vals,
      z = U.BxBy,
      text = tooltip,
      hoverinfo = 'text',
      colorscale = list(list(0, 'rgba(0,0,0,0)'), list(1, 'rgba(0,0,0,0)')),
      showscale = FALSE
    )
  
  # Handle practical significance not between initial slope and 0
  if (initial.slope > 0) {
    if (practical.significance > initial.slope) {
      warning("Parameter 'practical.significance' is greater than positive initial slope. Thus, there is is no region of significance.")
      return(fig)
    } else if (practical.significance < 0) {
      stop("Parameter 'practical.significance' is negative, but initial slope is positive. 'practical.significance' should have the same sign as the initial slope.")
    }
  } else if (initial.slope < 0) {
    if (practical.significance < initial.slope) {
      warning("Parameter 'practical.significance' is less than negative initial slope. Thus, there is is no region of significance.")
      return(fig)
    } else if (practical.significance > 0) {
      stop("Parameter 'practical.significance' is positive, but initial slope is negative. 'practical.significance' should have the same sign as the initial slope.")
    }
  }
  
  # Find points along practical.significance contour line
  contour.points <- NULL
  
  if (bound.index == 1) {
    for (Bx in Bx_vals) {
      # Set point to NA if practical.significance not achieved in widest interval
      widest.interval <- suppressWarnings(f(cor.res, sd.ratio.res, 0, (Bx - R2.wx) / (1 - R2.wx), 0, 1, -1, 1))
      if (widest.interval[bound.index] > practical.significance) {
        contour.points <- c(contour.points, NA)
        next
      }
      
      # Initialize values for this Bx
      interval <- c(R2.wy,1)
      upper_limit_x <- (Bx - R2.wx) / (1 - R2.wx)
      while (TRUE) {
        # Calculate estimate
        upper_limit_y <- (mean(interval) - R2.wy) / (1 - R2.wy)
        contour.estimate <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x, 0, upper_limit_y, -1, 1)[bound.index])
        
        # Break if boundary estimate found
        if (practical.significance - contour.estimate < .00001 & practical.significance - contour.estimate > -.00001) {
          contour.points <- c(contour.points, mean(interval))
          break
        }
        
        # Eliminate half of the interval for next run
        if (practical.significance - contour.estimate < 0) {
          interval <- c(mean(interval), interval[2])
        } else {
          interval <- c(interval[1], mean(interval))
        }
      }
    }
  } else {
    for (Bx in Bx_vals) {
      # Set point to NA if practical.significance not achieved in widest interval
      widest.interval <- suppressWarnings(f(cor.res, sd.ratio.res, 0, (Bx - R2.wx) / (1 - R2.wx), 0, 1, -1, 1))
      if (widest.interval[bound.index] < practical.significance) {
        contour.points <- c(contour.points, NA)
        next
      }
      
      # Initialize values for this Bx
      interval <- c(R2.wy,1)
      upper_limit_x <- (Bx - R2.wx) / (1 - R2.wx)
      while (TRUE) {
        # Calculate estimate
        upper_limit_y <- (mean(interval) - R2.wy) / (1 - R2.wy)
        contour.estimate <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x, 0, upper_limit_y, -1, 1)[bound.index])
        
        # Break if boundary estimate found
        if (practical.significance - contour.estimate < .00001 & practical.significance - contour.estimate > -.00001) {
          contour.points <- c(contour.points, mean(interval))
          break
        }
        
        # Eliminate half of the interval for next run
        if (practical.significance - contour.estimate > 0) {
          interval <- c(mean(interval), interval[2])
        } else {
          interval <- c(interval[1], mean(interval))
        }
      }
    }
  }
  
  # Find end points of practical.significance contour line
  interval.top <- c(R2.wx,1)
  upper_limit_y.top <- 1
  
  interval.right <- c(R2.wy,1)
  upper_limit_x.right <- max(Bx_vals) * .99 + .01
  
  while(TRUE) {
    # Calculate estimates
    upper_limit_x.top <- (mean(interval.top) - R2.wx) / (1 - R2.wx)
    upper_limit_y.right <- (mean(interval.right) - R2.wy) / (1 - R2.wy)
    estimate.top <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x.top, 0, upper_limit_y.top, -1, 1)[bound.index])
    estimate.right <- suppressWarnings(f(cor.res, sd.ratio.res, 0, upper_limit_x.right, 0, upper_limit_y.right, -1, 1)[bound.index])
    
    # Break if boundary estimate found
    if (practical.significance - estimate.top < .00001 & practical.significance - estimate.top > -.00001 &
        practical.significance - estimate.right < .00001 & practical.significance - estimate.right > -.00001) {
      break
    }
    
    # Eliminate half of the interval for next run
    if (bound.index == 1) {
      # Interval for top
      if (practical.significance - estimate.top < 0) {
        interval.top <- c(mean(interval.top), interval.top[2])
      } else {
        interval.top <- c(interval.top[1], mean(interval.top))
      }
      # Interval for right
      if (practical.significance - estimate.right < 0) {
        interval.right <- c(mean(interval.right), interval.right[2])
      } else {
        interval.right <- c(interval.right[1], mean(interval.right))
      }
    } else {
      # Interval for top
      if (practical.significance - estimate.top > 0) {
        interval.top <- c(mean(interval.top), interval.top[2])
      } else {
        interval.top <- c(interval.top[1], mean(interval.top))
      }
      # Interval for right
      if (practical.significance - estimate.right > 0) {
        interval.right <- c(mean(interval.right), interval.right[2])
      } else {
        interval.right <- c(interval.right[1], mean(interval.right))
      }
    }
  }
  
  # Use only those points that are not NA
  index <- which(!is.na(contour.points))
  
  # Add shaded region
  fig <- fig %>%
    add_trace(
      x = c(R2.wx, R2.wx, mean(interval.top), Bx_vals[index], 1, 1, R2.wx),
      y = c(R2.wy, 1, 1, contour.points[index], mean(interval.right), R2.wy, R2.wy),
      type = 'scatter',
      mode = 'lines',
      fill = 'toself',
      fillcolor = 'rgba(125, 125, 125, 0.5)',
      line = list(color = 'transparent'),
      name = "Region of Practical Significance",
      showlegend = FALSE
    )
  
  # Add tooltip
  fig <- fig %>%
    add_trace(
      type = 'heatmap',
      x = Bx_vals,
      y = By_vals,
      z = U.BxBy,
      text = tooltip,
      hoverinfo = 'text',
      colorscale = list(list(0, 'rgba(0,0,0,0)'), list(1, 'rgba(0,0,0,0)')),
      showscale = FALSE
    )
  return(fig)
}
