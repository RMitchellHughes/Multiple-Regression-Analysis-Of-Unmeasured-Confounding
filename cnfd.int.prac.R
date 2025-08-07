cnfd.int.prac <- function(cor.res, sd.ratio.res, R2.wy, R2.wx, practical.significance, granularity = 100, axis.labels = c('y', 'x')) {
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
  
  # Define tooltip
  tooltip <- matrix(
    paste0(
      "B<sub>", axis.labels[2], "</sub>: ", round(Bx_vals[col(U.BxBy)], 5), "<br>",
      "B<sub>", axis.labels[1], "</sub>: ", round(By_vals[row(U.BxBy)], 5), "<br>",
      "U(B<sub>", axis.labels[2], "</sub>,B<sub>", axis.labels[1], "</sub>) = ", round(U.BxBy, 5), "<br>",
      "L(B<sub>", axis.labels[2], "</sub>,B<sub>", axis.labels[1], "</sub>) = ", round(L.BxBy, 5)
    ),
    nrow = nrow(U.BxBy),
    ncol = ncol(U.BxBy)
  )
  
  initial.slope <- suppressWarnings(f(cor.res, sd.ratio.res, 0, 0, 0, 0, -1, 1)[1])
  Bx_vals <- seq(R2.wx, 1, length.out = granularity + 1)[1:granularity]
  contour.points <- NULL
  
  # Find points along practical.significance contour line
  if (practical.significance > 0) {
    bound.index <- 1
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
    plot.title <- plotly::TeX(paste0('\\large \\text{Region where } \\beta_{\\text{', axis.labels[2], '}. \\text{', axis.labels[1], '}|w,u} \\ge', practical.significance))
  } else {
    # User enters practical.significance = -4 > 0
    # What happens if practical significance is negative? => Use U(Bx,By)
    # # If initial.slope = -5 < practical significance, then line shows up
    # # If initial.slope = -3 > practical significance, then line doesn't show up (warning) (empty plot)
    # How do I know if contour line is in this column?
    # # If max(U(Bx,By)) < practical.significance, set to NA
    
    bound.index <- 2
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
    plot.title <- plotly::TeX(paste0('\\large \\text{Region where } \\beta_{\\text{', axis.labels[2], '}. \\text{', axis.labels[1], '}|w,u} \\le', practical.significance))
  }
  
  # Find end points of practical.significance contour line
  interval.top <- c(R2.wx,1)
  upper_limit_y.top <- 1
  
  interval.right <- c(R2.wy,1)
  upper_limit_x.right <- max(Bx_vals)
  
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
    if (practical.significance > 0) {
      if (practical.significance - estimate.top < 0) {
        interval.top <- c(mean(interval.top), interval.top[2])
      } else {
        interval.top <- c(interval.top[1], mean(interval.top))
      }
      if (practical.significance - estimate.right < 0) {
        interval.right <- c(mean(interval.right), interval.right[2])
      } else {
        interval.right <- c(interval.right[1], mean(interval.right))
      }
    } else {
      if (practical.significance - estimate.top > 0) {
        interval.top <- c(mean(interval.top), interval.top[2])
      } else {
        interval.top <- c(interval.top[1], mean(interval.top))
      }
      if (practical.significance - estimate.right > 0) {
        interval.right <- c(mean(interval.right), interval.right[2])
      } else {
        interval.right <- c(interval.right[1], mean(interval.right))
      }
    }
  }
  
  index <- which(!is.na(contour.points))
  
  fig <- plot_ly() %>%
    layout(
      title = list(text = plot.title),
      xaxis = list(title = plotly::TeX(paste0("\\Large R^2_{w,u;\\text{", axis.labels[2], "}}")),
                   range = c(0,1)),
      yaxis = list(title = plotly::TeX(paste0("\\Large R^2_{w,u;\\text{", axis.labels[1], "}}")),
                   range = c(0,1)),
      margin = list(t = 100, l = 100, r = 100, b = 100)
    )%>%
    config(
      mathjax = "cdn",
      displayModeBar = FALSE
    )
  
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
  print(practical.significance)
  return(fig)
}