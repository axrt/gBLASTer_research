getMLEandLoglike<-function (data, maxSteps = 50, weight = NULL) 
{
  if (missing(data)) 
    stop("A valid data set is required.")
  #Takes a transposed taxonomy data.frame and a transposed vector gStar
  #
  rowMatch <- function(A, B) {
    #Helper that simply pastes anything provided with ":"
    f <- function(...) paste(..., sep = ":")
    #Convert to matrix if not yet a matrix, which it is not, cuz it's a vector (presumably)
    if (!is.matrix(B)){
      B <- matrix(B, 1, length(B))
    } 
    a <- do.call("f", as.data.frame(A))
    b <- do.call("f", as.data.frame(B))
    match(b, a, nomatch = 0)
  }
  
  taxa <- rownames(data)
  #gStar is estimated as a mean of taxa
  if (is.null(weight)) {
    gStar <- apply(data, 1, mean)
  }
  else {
    gStar <- apply(data, 1, function(x) {
      p = (x %*% weight)/sum(weight)
    })
  }
  #iteratively move away the root from the mean
  while (rowMatch(t(data), t(gStar))[1] != 0){
    gStar <- gStar + 0.01
  } 
  #50 rows by default, 2 columns
  ret <- data.frame(matrix(0, nrow = maxSteps, ncol = 2))
  names(ret) <- c("f", "deltaf")
  
  #calculates the distance between all points and the gstar (euclid)
  #calc.f is an absolute sum of all distances in other word
  calc.f <- sum(apply(data, 2, function(x, gstart) {
    sqrt(sum((x - gstart)^2))
  }, gstart = gStar))
  
  
  delta.f <- 1
  count <- 0
  
  #delta is given as 10e-6, in other word, the algorithm stops if the difference betweent the previous and next step is smaller than delta
  #otherwise it stops after the number of counts runs out
  while ((count < maxSteps) && (delta.f > 10^(-6))) {
    count <- count + 1
    #ret number count is calculated value (funciton so far) and delta
    ret[count, ] <- c(calc.f, delta.f)
    
    #if weights given - account for the weights when calculating the distance
    if (is.null(weight)) {
      gStar <- rowSums(apply(data, 2, 
                             function(x, gstart) {
                               x/sqrt(sum((x - gstart)^2))
                               }, gstart = gStar))/
        sum(apply(data, 2, function(x, gstart) {
          1/sqrt(sum((x - gstart)^2))
          }, gstart = gStar))
      
      calc.f <- sum(apply(data, 2, function(x, gstart) {
        sqrt(sum((x - gstart)^2))
      }, gstart = gStar))
    }
    else {
      gStar <- (apply(data, 2, function(x, gstart) {
        x/sqrt(sum((x - gstart)^2))
      }, gstart = gStar) %*% weight)/as.vector(apply(data, 
                                                     2, function(x, gstart) {
                                                       1/sqrt(sum((x - gstart)^2))
                                                     }, gstart = gStar) %*% weight)
      calc.f <- as.vector(apply(data, 2, function(x, gstart) {
        sqrt(sum((x - gstart)^2))
      }, gstart = gStar) %*% weight)
    }
    
    delta.f <- abs(ret[count, 1] - calc.f)
  }
  
  
  N <- ncol(data)
  p <- nrow(data)
  if (is.null(weight)) {
    tau <- N/calc.f
  }
  else {
    tau <- sum(weight)/calc.f
  }
  gStar <- as.data.frame(gStar)
  if (is.null(weight)) {
    logLik <- N * (lgamma(p/2) - lgamma(p) + log(tau) - p * 
                     log(2) - p/2 * log(pi)) - tau * calc.f
  }
  else {
    logLik <- sum(weight) * (lgamma(p/2) - lgamma(p) + log(tau) - 
                               p * log(2) - p/2 * log(pi)) - tau * calc.f
  }
  mle.fit <- list(count, logLik, tau, gStar)
  names(mle.fit) <- c("iters", "Loglik", "tau", "mleTree")
  return(mle.fit)
}