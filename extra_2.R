library(ggplot2) 

#%% functions f ####### 

f <- list(
  list(
    f = function (x, y) {
      ( 4 - 2.1 * x**2 + x**4 / 3 ) * x**2 + x * y + 4 * ( y**2 - 1 ) * y**2 
    }, 
    xmin = -3, 
    xmax = 3, 
    ymin = -2, 
    ymax = 2 
  ),
  list(
    f = function (x, y) {
      exp(sin(50 * x)) + sin(60 * exp(y)) + sin(70 * sin(x)) + 
          sin(sin(80 * y)) - sin(10 * (x + y)) + (x * x + y * y) / 4 
    }, 
    xmin = -pi / 2, 
    xmax = pi / 2, 
    ymin = -pi / 2, 
    ymax = pi / 2 
  )
)

#%% Simulated Annealing algorithm ####### 

# we use the geometric schedule and Metropolis update with Cauchy proposal as it has heavier tails 
# since the region is restricted, a proposal can be outside the region, in this case there are two ways 
# to proceed: either correct it to the boundary of the region, or reject it 
# correcting to the boundary has some advantage if the minimum is at the boundary, however, this would mean 
# that the distribution of the proposal is not symmetric anymore; rejecting such proposals is equivalent to 
# assuming that the function f is +oo outside the region, so we keep the symmetry of the proposal 
# accordingly, we choose the second method, though it might potentially lead to longer convergence 
# the two parameters to tune are the multiplicator of the schedule a and the scale of Cauchy distribution s 
sa <- function (f, xmin, xmax, ymin, ymax, x0, y0, b0, a, s, n) {
  fabsmin <- f(x0, y0) 
  fabsmax <- fabsmin 
  xabsmin <- x0 
  yabsmin <- y0 
  xyf <- matrix(c(x0, y0, fabsmin, rep(0, n * 3)), ncol = 3, byrow = TRUE) 
  cs <- matrix(rcauchy(n * 2, scale = s), ncol = 2, byrow = TRUE) 
  b <- b0 
  for (i in 1:n) {
    b <- a * b 
    x <- xyf[i, 1] + cs[i, 1] 
    y <- xyf[i, 2] + cs[i, 2] 
    if (x < xmin || x > xmax || y < ymin || y > ymax) {
      reject = TRUE 
    }
    else {
      fxy <- f(x, y) 
      reject = runif(1) > exp(-b * (fxy - xyf[i, 3])) 
    }
    if (reject) {
      xyf[i + 1, ] = xyf[i, ] 
    }
    else {
      xyf[i + 1, ] = c(x, y, fxy) 
      if (fxy < fabsmin) {
        fabsmin <- fxy 
        xabsmin <- x 
        yabsmin <- y 
      } 
      else if (fxy > fabsmax) { 
        fabsmax <- fxy 
      } 
    }
  }
  xyf = data.frame(xyf) 
  colnames(xyf) <- c("x", "y", "f") 
  list(xyf = xyf, xabsmin = xabsmin, yabsmin = yabsmin, fabsmin = fabsmin, fabsmax = fabsmax) 
}

#%% main code ####### 

graphs <- function (fn, a, s, n) {
  saresult <- sa(
    f[[fn]]$f, f[[fn]]$xmin, f[[fn]]$xmax, f[[fn]]$ymin, f[[fn]]$ymax, 
    x0 = 0, y0 = 0, b0 = 1, 
    a = a, s = s, n = n
  )

  cat("The process led to : \n") 
  print(saresult$xyf[(n-1):(n+1), ]) 
  cat(sprintf("During the process the range of values was from %f (at %f, %f) to %f.\n", 
          saresult$fabsmin, saresult$xabsmin, saresult$yabsmin, saresult$fabsmax)) 
  
  midpoint = (saresult$fabsmax + saresult$fabsmin) / 2 
  list( 
    ggplot(saresult$xyf, aes(x = x, y = y, color = f)) + geom_path() + 
        scale_color_gradient2(low = "darkgreen", mid = "yellow", high = "red", midpoint = midpoint), 
    ggplot(saresult$xyf, aes(x = as.numeric(row.names(saresult$xyf)), y = f, color = f)) + geom_path() + 
        scale_color_gradient2(low = "darkgreen", mid = "yellow", high = "red", midpoint = midpoint) 
  ) 
}

graphs(2, 1.001, 1/3, 10000) 
