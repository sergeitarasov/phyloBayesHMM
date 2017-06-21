# if we use adaptive MCMC we might have to modify our parameters; from TESSS
if ( adaptive == TRUE ) {
  if ( i <= burnin ) {## !!! I need to use this line in my code
    if ( i %% OPTIMIZATIONS == 0 ) {
      for ( j in 1:length(parameters) ) {
        rate <- accepted[j] / tried[j]
        if ( rate > 0.44 ) {
          delta[j] <- delta[j] * (1.0 + ((rate-0.44)/0.56) )
        } else {
          delta[j] <- delta[j] / (2.0 - rate/0.44 )
        }
        tried[j] <- 0
        accepted[j] <- 0
      }
    }
  }
}

#########

#M=matrix(c(-.02,.05,.02,-.05),2,2)
#calclikhood=function(treeedge, treeNnode, data, edgelength, base_vector, rate_param){
#  M=matrix(c(-rate_param[1],rate_param[2],rate_param[1],-rate_param[2]),2,2)
#  likelihood_root(treeedge, treeNnode, data, edgelength, M, base_vector)}
#calclikhood(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length,  c(.5,.5), param)

likelihood_root(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length, M, c(.5,.5))
likhood <- likelihood_root(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length, M, c(.5,.5))
# Set initial parameter values:
# Read in the priors:
prior_r1=2
prior_r2=2
#curve(dexp(x, 2))
# Parameters for MH updates (standard deviation for normal random walk):
accepted=c(0,0)
#tried
###iterations
nt = 10
nparam = 2

deltalog=c() ##log of tuning parameter delta
delta <- c(0.5,0.5)
accepted=c(0,0)#increment log if parameter is accepted
param <- c(0.5,0.5)#initial values of parameters
adapt_check=300#tune adaptation afern specified number of generations
M=matrix(c(-param[1],param[2],param[1],-param[2]),2,2)



  nt=1000
  itns <- array(0, dim=c(nt, nparam+3))# define array to keep track of all parameters

   output <- dim(nparam+1) # output - vector for parameter values and  log-likelihood:


   # MCMC updates - MH algorithm - cycle through each iteration:
    for (t in 1:nt){
    # Update the parameters in the model using function "updateparam":
    output <-  updateparam(nparam,param,likhood,prior_r1,delta,my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length,  c(.5,.5), accepted)
    # Set parameter values and log-likelihood to be the output from
    # the MH step:
    param <- output[1:nparam]
    likhood <- output[nparam+1]
    accepted <- output[4:5]
    itns[t,] <- output
    #write(output,"param-test-wrting.txt", append=T)

## tuning of delta parameter
    if ( t %% adapt_check == 0 ) {
      for ( j in 1:nparam) {
        rate <- accepted[j] / adapt_check
        if ( rate > 0.44 ) {
          delta[j] <- delta[j] * (1.0 + ((rate-0.44)/0.56) )
        } else {
          delta[j] <- delta[j] / (2.0 - rate/0.44 )
        }
                accepted[j] <- 0
      }
      deltalog=rbind(deltalog, delta)
    }}
###############

  plot(deltalog[,1], type="l")
  itns[which(itns[,3]==max(itns[,3])),]

  plot(itns[,1], itns[,3])
  plot(itns[,2], type="l")
  hist(itns[,1])
  mean(itns[100000:300000,2])
  mean(itns[100000:300000,1])

  plot(itns[8000:10000,2], type="l")

  plot(density(itns[100000:300000,3]),type="l",col="darkgreen",lwd=4)

  #updateparam(nparam,param,likhood,prior_r1,delta,my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length,  c(.5,.5))

  #treeedge=my.tree$edge
  #treeNnode=my.tree$Nnode
 # data=my.data
 # edgelength=my.tree$edge.length
 # base_vector=c(.5,.5)
 # param


# Function for performing single-update MH algorithm:
updateparam <- function(nparam,param, likhood, prior_r1,
                        delta,    treeedge, treeNnode, data, edgelength, base_vector, accepted){
  for (i in 1:nparam) {
    # Keep a record of the current parameter value being updated
    # Propose a new value using a random walk from normal proposal
    oldparam <- param[i]
    param[i] <- rnorm(1, param[i], delta[i])
    # If proposed value is in [0,1] calculate acceptable probability
    if (param[i] >= 0) {
      M=matrix(c(-param[1],param[2],param[1],-param[2]),2,2)
      newlikhood <- likelihood_root(treeedge, treeNnode, data, edgelength, M, base_vector)
        #calclikhood(treeedge, treeNnode, data, edgelength, base_vector,  param)
        #calclikhood(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length,  c(.5,.5), param)



      # In acceptance probability include likelihood and prior contributions
      # (proposal contributions cancel since it is a symmetric distribution)
      num <- newlikhood
      den <- likhood
     # num <- newlikhood + log(dexp(param[i],prior_r1))
      #den <- likhood +log(dexp(oldparam,prior_r1))
      A <- min(1,exp(num-den))
    }
    else { A <- 0 }
    # Simulate a random number in [0,1] and accept move with probability A;
    # else reject move and return parameter value to previous value
    u <- runif(1)
    if (u <= A) { likhood <- newlikhood; accepted[i]=accepted[i]+1 }
    else { param[i] <- oldparam }
  }
  # Set the values to be outputted from the function to be the
  # parameter values and log-likelihood value. Then output the values.
  output <- c(param, likhood, accepted)
  output
}
