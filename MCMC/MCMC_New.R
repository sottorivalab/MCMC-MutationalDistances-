#!/usr/bin/Rscript




data <- read.table("Distribution.txt")  # Importing the data file



repeats <- 1000 # Length of the Markov Chain

samplesize <- 1469 # Number of distances sampled


# Defining some initial values for the parameters of the distributions

N <- nrow(data)  # Length of the distribution (Maximum number of mutations)
mu <- 10        # Mutation rate, takes values between 0 and infinity
beta <- 0.1      # Survival rate of cells, takes values between 0 and 1

p <- 1         # Proliferation rate of cells, set to 1 (Time is measured in generations)



r1 <- 10      # Defining the cut off of the infinte sums of the probability density functions, It depends on the details of the data how many terms we need. Higher terms       contribute to larger distances.
k1 <- 30



##############



# Storing the MCMC parameters

list_mu <- c(mu)  # List of mutation rates
list_beta <- c(beta) # List of cell survival rates
list_p <- c(p)  # List of cell proliferation rates (kept as 1)

list_Li <- c(0) # List of proposed likelihoods
list_Li1 <- c(0) # List of accepted likelihoods



############




# Calculating the normalization factor for the probability density

# We use a recursive implementation of the probability density to avoid numerical errors due to the large Binomial coefficients and the small exponential functions


norm <- 0

for(j in 1:N)
{ 
for(r in 1:r1)
{
for(k in r:k1)
{ prod <- 1
for(l in 1:j)
{ prod <- prod * (mu*k/l)*exp(- mu*k/j) }
norm <- norm +   (exp( - exp( - p*beta*(r+1))/(beta * p)) - exp( - exp( - p*beta*(r))/(beta * p)))/(1-exp( - exp( - p*beta)/(beta * p)))    *choose(k-1,r-1)*beta^(r)*(1-beta)^(k-r)*prod
}
}
}

scaling <- norm

#################


# Calculating the likelihood of the normalized probability density function with initial parameters given the data

size <- length(data[,1])

prob <- 0 

for(j in 1:size)
{ 
if (data[j,2] >0 )
{ test <- 0
for(r in 1:r1)
{
for(k in r:k1)
{ prod <- 1
for(l in 1:j)
{ prod <- prod * (mu*k/l)*exp(- mu*k/j) }
test <- test +  (exp( - exp( - p*beta*(r+1))/(beta * p)) - exp( - exp( - p*beta*(r))/(beta * p)))/(1-exp( - exp( - p*beta)/(beta * p)))   *choose(k-1,r-1)*beta^(r)*(1-beta)^(k-r)*prod/scaling
}
}
prob <- prob + data[j,2]*log(test)
}
}

################



# Running the MCMC algorithm

for(w in 1:repeats)
{



# Creating a random set of model parameters, in "proximity" of the old parameters

v1 <- mu + rnorm(1,0,0.15)
v2 <- beta +runif(1,-0.06,0.06)
  #  v3 <- p  +runif(1,-0.06,0.06)

##############


# Making sure, the new parameters are in a meaningfull range

if (v1 > 0)
{mu1 <- v1}
if(v2 > 0.05 & v2 <= 1)
{beta1 <- v2}
 # if(v3 > 0.05 & v3 < 1)
p1 <- p # {p1 <- v3}

############



# The new probability density (new sets of parameters) needs to be normalized again

norm <- 0

for(j in 1:N)
{ 
for(r in 1:r1)
{
for(k in r:k1)
{ prod <- 1
for(l in 1:j)
{ prod <- prod * (mu1*k/l)*exp(- mu1*k/j) }
norm <- norm + (exp( - exp( - p1*beta1*(r+1))/(beta1 * p1)) - exp( - exp( - p1*beta1*(r))/(beta1 * p1)))/(1-exp( - exp( - p1*beta1)/(beta1 * p1)))  *choose(k-1,r-1)*beta1^(r)*(1-beta1)^(k-r)*prod
}
}
}

scaling <- norm 

#############


# Calcualting the likelihood of the new distribution given the data again

size <- length(data[,1])

prob1 <- 0 

for(j in 1:size)
{ 
if (data[j,2] >0 )
{ test <- 0
for(r in 1:r1)
{
for(k in r:k1)
{ prod <- 1
for(l in 1:j)
{ prod <- prod * (mu1*k/l)*exp(- mu1*k/j) }
test <- test +  (exp( - exp( - p1*beta1*(r+1))/(beta1 * p1)) - exp( - exp( - p1*beta1*(r))/(beta1 * p1)))/(1-exp( - exp( - p1*beta1)/(beta1 * p1)))  *choose(k-1,r-1)*beta1^(r)*(1-beta1)^(k-r)*prod/scaling
}
}
if (test > 0) # Avoiding a likelihood of 0
{prob1 <- prob1 + data[j,2]*log(test)}
else
{}
}
}

if(test > 0)
{ratio <- exp(samplesize*(prob1 - prob)) # Comparing the likelihood of the new set of parameters (prob1) and the old set of parameters

############
    
list_Li <- append(list_Li,prob1)


if ( ratio > runif(1,0,1))
{ 
mu <- mu1
beta <- beta1
p <- p1
prob <- prob1
list_mu <- append(list_mu,mu)
list_beta <- append(list_beta,beta)
list_p <- append(list_p,p)

list_Li1 <- append(list_Li1,prob1)



}
}
else
{}

###########

print(w)

}

###############


write.table(list_mu, "list_mu3.txt", sep="\t")
write.table(list_beta, "list_beta3.txt", sep="\t")
write.table(list_Li, "list_Li.txt", sep="\t")

write.table(list_Li1, "list_Li1.txt", sep="\t")


