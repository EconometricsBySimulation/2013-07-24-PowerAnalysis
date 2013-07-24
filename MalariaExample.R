require("AER")   # We will use the ivreg function later

# We define a function that both simulates a population and estimates 
# results and returns and returns the results.
MalariaIV <- function(
  nsim = 100,
  npop = 1000,    # Define the sampling population
  treat0 = 1/3,   # Define the proportion which is control
  treat1 = 1/3,   # Define the proportion which is treated with free bednets
  treat2 = 1/3,   # Define the proportion which is treated with 50% cost
  t0comp = .05,   # Bednet useage among control group
  t1comp = .85,   # Bednet useage if recieving free net
  t2comp = .5 ,   # Bednet useage if recieving 50% cost
  malariaRT = .85,# Rate of getting malaria without bednet
  netRT = .45,    # Rate of getting malaria with bednet
  Tdetect = 1,    # Likihood of detecting malaria if present
  Fdetect = 0,    # Likihood of detecting malaria if not present
  alpha = .05     # The alpha level that a p-value must be below
) {
  # Define how the function works
  # First, define a vector to store results
  pvalues <- NULL
  
  #
  for (i in 1:nsim) { # Repeat the simulation a number of times
    # Generate the population by technology useage
    simdata <- data.frame(
      control=c(rep(1,npop*treat0),rep(0,npop*treat1), rep(0,npop*treat2)),
      treat1 =c(rep(0,npop*treat0),rep(1,npop*treat1), rep(0,npop*treat2)),
      treat2 =c(rep(0,npop*treat0),rep(0,npop*treat1), rep(1,npop*treat2)))
    # Calculate the actual population generated (should be 999 in this case)
    npeople <- nrow(simdata)
    # Generate the rate of bednet usage.
    simdata$bednet <-
      simdata$control*rbinom(npeople,1,t0comp) + # Bednet usage control
      simdata$treat1 *rbinom(npeople,1,t1comp) + # Free
      simdata$treat2 *rbinom(npeople,1,t2comp)   # 50% cost
    # Now let's generate the rate of malaria
    simdata$malaria <-
      (simdata$bednet==0)*rbinom(npeople,1,malariaRT)+
      (simdata$bednet==1)*rbinom(npeople,1,netRT)
    # Finally generate the rate of malaria detection as a function
    # of true parasite levels.
    simdata$mdetect <- 
      (simdata$malaria==0)*rbinom(npeople,1,Fdetect)+ # False detection
      (simdata$malaria==1)*rbinom(npeople,1,Tdetect)  # True detection
    
    # Time to do our simple estimation of the effect of treatment on bednet use
    lmcoef <- summary(lm(bednet~treat1+treat2,data=simdata))$coefficients
    # Now let's try the 2SLS to estimate the effect of bednets on contraction
    # of malaria.
    
    ivregest <- ivreg(mdetect~bednet | treat1+treat2, data=simdata)
    ivregcoef <- summary(ivregest)$coefficients

    # Save the rejection of the null rates
    pvalues <- rbind(pvalues,c(treat1=lmcoef[2,4]<alpha,
                               treat2=lmcoef[3,4]<alpha,
                               iv=ivregcoef[2,4]<alpha))
    
  }
  # Calculate the mean rejection rate for each coefficient.
  apply(pvalues,2,mean)
}

MalariaIV() 
# Running it once we can see that we get a single set of results
# where we easily reject the null.  However, we want to know
# what happens when bednets are not so effective or malaria is harder
# to detect.  We can modify the parameters fed into the model to test
# these questions manually or we could use the SimpleSim function from
# a previous post.

SimpleSim(fun=MalariaIV, 
          npop=c(100,1000),
          t0comp=c(.05,.25),
          t1comp=c(.5,.75),
          t2comp=c(.25,.5),
          malariaRT=c(.85,.5),
          netRT=c(.75,.45),
          Tdetect=c(1,.7),
          Fdetect=c(0,.3),
          alpha=.05,
          nsim=10)
# This could take a little while to run since there are 256 combinations
# to try and each of them will be run 10 times.

# The above command gives back lots of data but it is not always very easy
# to understand in a matrix form.  It is often easier to just vary one
# paramter at a time.

sample.size <- SimpleSim(fun=MalariaIV, 
          npop=c(500,1000*1:10),
          t0comp=c(.25),
          t1comp=c(.75),
          t2comp=c(.375),
          malariaRT=c(.5),
          netRT=c(.35),
          Tdetect=c(.8),
          Fdetect=c(.2),
          alpha=.05,
          nsim=200)

# Looking at just sample size
require(ggplot2)

# Save the results to single long data format to be useable by ggplot2
results <- with(sample.size, data.frame(reject = as.numeric(
  c(iv,treat2,treat1)),
  id=rep(c("iv","treat2","treat1"), each=length(npop)),
  npop))

p <- ggplot(results, aes(npop, reject))
# Normally a 80% detection rate is the minimum rejection rate needed 
# to justify a study. 
p + geom_point(aes(colour =id)) + 
  geom_line(aes(group=id)) + 
  geom_hline(yintercept = .8) 

# Thus we want to ensure our study has at least 5000 participants.

# I think there might be some additional considerations for the
# number of participants required to ensure that there is no
# false rejection of the null.  However, I am no expert on Power 
# Analysis so this is what I know.
