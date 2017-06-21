# Define the function with input parameters:
# nt = number of iterations; nburn = burn-in

nt = 10
Rdippercode <- function(nt,nburn){
  # Define the parameter values:
  # ni = number of release years; nj = number of recapture years
  # nparam = maximum number of parameters
  ni = 6
  nj = 6
  nparam = 2
  # Read in the data:
  data <- matrix(c(
    11,2,0,0,0,0,9,
    0,24,1,0,0,0,35,
    0,0,34,2,0,0,42,
    0,0,0,45,1,2,32,
    0,0,0,0,51,0,37,
    0,0,0,0,0,52,46),nrow=ni,byrow=T)
  # Read in the priors:
  alpha <- c(1,1)
  beta <- c(1,1)
  # Parameters for MH updates (standard deviation for normal random walk):
  delta <- c(0.05,0.05)

  > dip <- dippercodeMHnorm(10000,1000)

  # Set initial parameter values:
  param <- c(0.9,0.5)
  # Calculate log-likelihood for initial state using function
  # "calclikhood":
  likhood <- calclikhood(ni, nj, data, param)
  # Define itns - array to store sample from posterior distribution;
  # output - vector for parameter values and associated log-likelihood:
  itns <- array(0, dim=c(nt, nparam))
  output <- dim(nparam+1)
  # MCMC updates - MH algorithm - cycle through each iteration:
  for (t in 1:nt){
    # Update the parameters in the model using function "updateparam":
    output <- updateparam(nparam,param,ni,nj,data,likhood,
                          alpha,beta,delta)
    # Set parameter values and log-likelihood to be the output from
    # the MH step:
    param <- output[1:nparam]
    likhood <- output[nparam+1]
    for (i in 1:nparam) {
      itns[t,i] <- param[i]
    }
  }


  # Remove the burn-in from the simulations and calculate the mean and
  # standard deviation of the parameters
  subitns <- itns[(nburn+1):nt,]
  mn <- array(0,nparam)
  std <- array(0,nparam)
  for (i in 1:nparam) {
    mn[i] <- mean(subitns[,i])
    std[i] <- sd(subitns[,i]))
  }
  # Output the posterior mean and standard deviation of the parameters
  # following burn-in to the screen:
  cat("Posterior summary estimates for each model:
cat("\n")
cat("mean (SD)", "\n")
cat("p: ", "\n")
cat(mn[1], "
      (", std[1], ")", "\n")
cat("\n")
cat("phi: ", "\n")
cat(mn[2], "
      (", std[2], ")", "\n")
", "\n")
  # Output the sample from posterior distribution and close the function:
  itns
}


# This function calculates the log-likelihood of capture-recapture data
calclikhood <- function(ni, nj, data, param){
  # Set up the arrays for parameters (phi and p), cell probabilities (q)
  # Set parameter values for each year of the study likhood to be zero
  phi <- array(0,nj)
  p <- array(0,nj)
  q <- array(0,dim=c(ni,nj+1))
  for (i in 1:nj) {
    p[i] <- param[1]
    phi[i] <- param[2] }
  likhood <- 0
  # Calculate multinomial cell probabilities and log-likelihood value
  for (i in 1:ni){
    # For diagonal elements:
    q[i,i] <- phi[i]*p[i]
    likhood <- likhood + data[i,i]*log(q[i,i])
    # Calculate the elements above the diagonal:
    if (i <= (nj-1)) {
      for (j in (i+1):nj) {
        q[i,j] <- prod(phi[i:j])*prod(1-p[i:(j-1)])*p[j]
        likhood <- likhood + data[i,j]*log(q[i,j]) } }
    # Probability of an animal never being seen again
    q[i,nj+1] <- 1 - sum(q[i,i:nj])
    likhood <- likhood + data[i,nj+1]*log(q[i,nj+1])
  }
  # Output the log-likelihood value:
  likhood
}

likhood=0
i=1
# Function for performing single-update MH algorithm:
updateparam <- function(nparam,param,ni,nj,data,likhood,alpha,beta,
                        delta){
  for (i in 1:nparam) {
    # Keep a record of the current parameter value being updated
    # Propose a new value using a random walk from normal proposal
    oldparam <- param[i]
    param[i] <- rnorm(1, param[i], delta[i])
    # If proposed value is in [0,1] calculate acceptable probability
    if (param[i] >= 0 & param[i] <= 1) {
      newlikhood <- calclikhood(ni, nj, data, param)
      # In acceptance probability include likelihood and prior contributions
      # (proposal contributions cancel since it is a symmetric distribution)
      num <- newlikhood + log(dbeta(param[i],alpha[i],beta[i]))
      den <- likhood + log(dbeta(oldparam,alpha[i],beta[i]))
      A <- min(1,exp(num-den))
    }
    else { A <- 0 }
    # Simulate a random number in [0,1] and accept move with probability A;
    # else reject move and return parameter value to previous value
    u <- runif(1)
    if (u <= A) { likhood <- newlikhood }
    else { param[i] <- oldparam }
  }
  # Set the values to be outputted from the function to be the
  # parameter values and log-likelihood value. Then output the values.
  output <- c(param, likhood)
  output
}
