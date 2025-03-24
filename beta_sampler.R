# a = 10
# b = 0.2
# 
# x = rbeta(100, a, b)
# theta = beta_sample(x, 1e4)
# profvis::profvis({theta = beta_sample(x, 1e4)})
# 
# s = system.time({theta = beta_sample(x, 1e4)})
# 
# s2 = system.time({meep = beta.ab.test(x, 1e3, 1e4, ab = T)})

beta_sample <- function(x, nsamples, burnin = 1e3){
  
  # Sample a's and b's
  rv = list()
  theta_sample = matrix(0, nrow = nsamples, ncol = 2)
  a_accept_rate = 0
  b_accept_rate = 0
  
  a = 1
  b = 1
  
  n = length(x)
  eps = .Machine$double.eps
  Slogx  = sum(log(x+eps))
  Slogcx = sum(log(1-x+eps))
  
  # Initial coordinate-wise ML search
  for(i in 1:5){
    
    a_approx = beta_approx(Slogx, Slogcx, n, b, a)
    a = a_approx$a_mode
    
    b_approx = beta_approx(Slogcx, Slogx, n, a, b)
    b = b_approx$a_mode
    
  }
  
  ## Sampling loop
  iter = 1
  iter_sample = 1
  while(iter_sample <= nsamples){
    
    # Gibbs steps
    a_approx = beta_approx(Slogx, Slogcx, n, b, a, x)
    a = a_approx$a_sample
    #a_mode = a_approx$a_mode
    
    b_approx = beta_approx(Slogcx, Slogx, n, a, b)
    b = b_approx$a_sample
    #b_mode = b_approx$a_mode
    
    # Store samples
    if (iter > burnin) {
      theta_sample[iter_sample, ] = c(a, b)
      iter_sample = iter_sample + 1
      
      if(a_approx$accept){
        a_accept_rate = a_accept_rate + 1
      }
      
      if(b_approx$accept){
        b_accept_rate = b_accept_rate + 1
      }
    }
    
    iter = iter + 1
    
  }
  
  rv$a_accept_rate = a_accept_rate/nsamples
  rv$b_accept_rate = b_accept_rate/nsamples
  rv$theta_sample = theta_sample
  
  if(a_approx$max_reached){
    cat('a')
  }
  
  if(b_approx$max_reached){
    cat('b')
  }
  
  return(rv)
  
}

beta_approx <- function(Slogx, Slogcx, n, b, a_prev, x){
  
  rv = list()
  rv$accept = F
  
  ## Get ML estimate and Hessian and use as starting gamma approximation
  a_est = a_est_beta(Slogx,Slogcx,n,b)
  y = a_est$y
  
  ## Fit the gamma distribution to the observed conditional log-probabilities
  mean_vector <- colMeans(a_est$t)
  t = a_est$t - matrix(mean_vector, nrow(a_est$t), ncol(a_est$t), byrow = TRUE) # center the estimates (so no need to fit intercept)
  
  #t = scale(a_est$t, center = T, scale = F)
     
  #eta <- lm(y ~ t - 1)$coefficients
  #eta = qr.solve(t, y)
  
  XtX = crossprod(t)
  XtX_a = XtX[1,1]
  XtX_b = XtX[1,2]
  XtX_c = XtX[2,1]
  XtX_d = XtX[2,2]
  
  det_XtX <- ( XtX_a * XtX_d) - (XtX_b * XtX_c)
  
  if(is.na(det_XtX)){
    browser()
  }
  
  # Check if the determinant is non-zero
  if (det_XtX != 0) {
    # Calculate the inverse of the matrix
    inv_XtX <- matrix(c(XtX_d, -XtX_b, -XtX_c, XtX_a), nrow = 2, ncol = 2, byrow = TRUE)/ det_XtX
    
  } else {
    break
    print("The matrix is singular and does not have an inverse.")
  }
  
  eta = inv_XtX %*% crossprod(t,y)
  
  ## Return fitted gamma model parameters (transform the natural parameter back to the well known shape, k and scale, theta parameter)
  k_fit = eta[1] + 1
  theta_fit = -1 / eta[2]
  
  ## Metropolis-Hastings sampling step
  
  # generate a random sample from the fitted distribution (proposal)
  a_new = rgamma(1, shape = k_fit, scale = theta_fit)
  
  if(is.na(a_new)){
    print('here1')
  }
  
  #num
  #prod((x^(a_new-1))/beta(a_new,b))/(1+a_new^2)
  #a_prev^(k_fit - 1)*exp(-a_prev/theta_fit)
  #prod((x^(a_new-1))/beta(a_new,b))/(1+a_new^2) * a_prev^(k_fit - 1)*exp(-a_prev/theta_fit)
  
  #den
  #prod((x^(a_prev-1))/beta(a_prev,b))/(1+a_prev^2) * a_new^(k_fit - 1)*exp(-a_new/theta_fit)
  num = beta_L(Slogx, n, a_new, b) + gamNLL(a_prev, k_fit, theta_fit)
  den = beta_L(Slogx, n, a_prev, b) + gamNLL(a_new, k_fit, theta_fit)
  p = min(1, exp(-(num-den)))
  
  # Accept/reject?
  a_sample = a_prev
  if (runif(1) < p){
    a_sample = a_new
    rv$accept = T
  }
  
  rv$k_fit = k_fit
  rv$theta_fit = theta_fit
  rv$a_sample = a_sample
  rv$a_mode = a_est$a
  rv$max_reached = a_est$max_reached
  #rv$iter = a_est$iter
  #cat(a_est$iter, "; ")
  
  return(rv)
}


a_est_beta <- function(Slogx, Slogcx, n, b, max.iter = 30){
  
  rv = list()
  rv$max_reached = F
  
  ## Initial guess
  Gx = Slogx/n
  Gcx = Slogcx/n
  a = 1/2 + Gx/2/(1-Gx-Gcx)
  a = max(a,1e-4)
  
  iter = 0
  c = 0.9
  kappa = 1
  
  while (T) {
    
    ## Newton-Raphson update
    g = ((digamma(a) - digamma(a + b)) * n - Slogx + 2 * a / (a^2 + 1))
    
    # Compute the Hessian H
    H = (n * (trigamma(a) - trigamma(a + b)) + 2 / (a^2 + 1) - 4 * a^2 / (a^2 + 1)^2)
    
    a_new = a - kappa*g/H
    
    if(a_new < 0){
      kappa = kappa * c
      next
    }
    
    if( (abs(a - a_new) / abs(a) < 1e-5)){
      break
    }
    
    if(iter > max.iter){
      rv$max_reached = T
      break
    }
    
    a = a_new
    iter = iter + 1
  }
  
  ## Compute the Hessian at the maximum
  H =  n * (trigamma(a) - trigamma(a + b))  # + 2/(a^2+1) - 4*a^2/(a^2+1)^2
  
  ## Create design points
  a_d = c(max(a - sqrt(2/H), a/2), a, a + sqrt(2/H))
  rv$t = matrix(c(log(a_d), a_d), ncol = 2)    # sufficient statistic
  rv$y = -beta_L(Slogx, n, a_d, b)     # compute the conditionals up to a constants at the 3 points
  
  rv$H = H
  rv$a = a
  rv$iter = iter
  #cat(iter, "; ")
  
  return(rv)
  
}

# Beta neg-log-posterior
beta_L <- function(Slogx, n, a, b){
  
  return( -(a-1)*Slogx + n*lbeta(a,b) + log(1+a^2))
  
}

# Gamma neg-log probability (proposal)
gamNLL <- function(x, k, theta){
  
  return(-(k-1)*log(x) + x/theta)
}



