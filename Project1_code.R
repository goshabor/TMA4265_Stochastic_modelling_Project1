#Task 1c
beta <- 0.01
gamma <- 0.10
alpha <- 0.005

P <- matrix(c(1-beta, beta, 0,
              0, 1-gamma, gamma,
              alpha, 0, 1-alpha), nrow=3, byrow=TRUE)

simulate_markov_chain <- function(P, n, initial_state=1) {
  current_state <- initial_state
  states <- numeric(n)
  states[1] <- current_state
  
  for (i in 2:n) {
    current_state <- sample(1:3, 1, prob=P[current_state,])
    states[i] <- current_state
  }
  
  return(states)
}

get_estimates <- function(simulation) {
  last_10_years <- simulation[(365*10 + 1):length(simulation)]
  
  S_days <- sum(last_10_years == 1)
  I_days <- sum(last_10_years == 2)
  R_days <- sum(last_10_years == 3)
  
  total_days <- length(last_10_years)
  
  c(S=S_days/total_days, I=I_days/total_days, R=R_days/total_days)
}

n_runs <- 30
n_days <- 20 * 365
estimates <- replicate(n_runs, get_estimates(simulate_markov_chain(P, n_days)))

# Compute 95% CI for each state
ci_S <- quantile(estimates[1,], c(0.025, 0.975))
ci_I <- quantile(estimates[2,], c(0.025, 0.975))
ci_R <- quantile(estimates[3,], c(0.025, 0.975))

sprintf("Confidence interval for S: %f", ci_S * 365)
sprintf("Confidence interval for I: %f", ci_I * 365)
sprintf("Confidence interval for R: %f", ci_R * 365)

#Task 1e)
# Model the spread of infection over time with no vaccinations

# Function to simulate the measles outbreak
infections = function(N, days, Y0){
  
  # Define transmission rates
  alpha = 0.005
  gamma = 0.1
  
  # Initialize time, susceptible, infected, and recovered vectors
  n = seq(1, days)
  S = rep(0, days)
  I = rep(0, days)
  R = rep(0, days)
  
  # Set initial conditions
  S[1] = Y0[1]
  I[1] = Y0[2]
  R[1] = Y0[3]
  
  # Iterate through each day and update counts
  for (i in 1:(days - 1)) {
    beta = 0.5 * I[i] / N
    S_to_I = rbinom(1, S[i], beta)
    I_to_R = rbinom(1, I[i], gamma)
    R_to_S = rbinom(1, R[i], alpha)
    
    S[i + 1] = S[i] - S_to_I + R_to_S
    I[i + 1] = I[i] + S_to_I - I_to_R
    R[i + 1] = R[i] + I_to_R - R_to_S
  }
  
  return(list(n, S, I, R))
}

# Run the simulation
result = infections(1000, 300, c(950, 50, 0))
# Extract the output results
n = result[[1]]
S = result[[2]]
I = result[[3]]
R = result[[4]]

# Plot the infection dynamics
plot(n, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", 
     ylim = c(0, 1000), main = "Measles Outbreak")
lines(n, I, col = "red")
lines(n, R, col = "green")
legend("topright", legend = c("S", "I", "R"), col = c("blue", "red", "green"), lty = 1)

#Task 1f)
# Analyzing the peak infection times and sizes

# Arrays to store peak infection times and sizes
I_max = c()
idx = c()

# Simulate the outbreak 1000 times
for (i in 1:1000) {
  result = infections(1000, 300, c(950, 50, 0))
  I = result[[3]]
  # Record peak infection size
  I_max = c(I_max, max(I, na.rm = TRUE))
  # Record time of peak infection
  idx = c(idx, which.max(I))
}

# Compute confidence intervals and means for peak infection sizes and times
result_I_max = t.test(I_max, na.rm = TRUE)
confidence_interval_I_max = result_I_max$conf.int
mean_I_max = mean(I_max, na.rm = TRUE)
result_idx = t.test(idx, na.rm = TRUE)
confidence_interval_idx = result_idx$conf.int
mean_idx = mean(idx, na.rm = TRUE)

# Output the results
confidence_interval_I_max
mean_I_max
confidence_interval_idx
mean_idx

#Task 1g)
# Modeling the spread of infection with different vaccination counts


# Function to simulate the measles outbreak with vaccinations
infections_v = function(N, v, days, Y0){
  # Define transmission rates
  alpha = 0.005
  gamma = 0.1
  
  # Initialize time, susceptible, infected, and recovered vectors
  n = seq(1, days)
  S = rep(0, days)
  I = rep(0, days)
  R = rep(0, days)
  
  # Set initial conditions accounting for vaccinations
  S[1] = Y0[1] - v
  I[1] = Y0[2]
  R[1] = Y0[3]
  
  # Iterate through each day and update counts
  for (i in 1:(days - 1)) {
    beta = 0.5 * I[i] / N
    S_to_I = rbinom(1, S[i], beta)
    I_to_R = rbinom(1, I[i], gamma)
    R_to_S = rbinom(1, R[i], alpha)
    S[i + 1] = S[i] - S_to_I + R_to_S
    I[i + 1] = I[i] + S_to_I - I_to_R
    R[i + 1] = R[i] + I_to_R - R_to_S
  }
  
  return(list(n, S, I, R))
}

# Run simulations for different vaccination numbers
result1 = infections_v(1000, 800, 300, c(950, 50, 0))
result2 = infections_v(1000, 600, 300, c(950, 50, 0))
result3 = infections_v(1000, 100, 300, c(950, 50, 0))
n = result1[[1]]
I_800 = result1[[3]]
I_600 = result2[[3]]
I_100 = result3[[3]]

# Plot the infection dynamics for different vaccination counts
plot(n, I_800, type = "l", col = "red", xlab = "Time", ylab = "Infected", 
     ylim = c(0, max(I_100)), main = "Measles Outbreak Vaccinated")
lines(n, I_600, col = "blue")
lines(n, I_100, col = "green")
legend("topright", legend = c("v=800", "v=600", "v=100"), col = c("red", "blue", "green"), lty = 1)

# Arrays to store peak infection sizes for different vaccination numbers
I_max1 = c()
I_max2 = c()
I_max3 = c()

# Simulate the outbreak 1000 times for each vaccination number
for (i in 1:1000) {
  result1 = infections_v(1000, 800, 300, c(950, 50, 0))
  result2 = infections_v(1000, 600, 300, c(950, 50, 0))
  result3 = infections_v(1000, 100, 300, c(950, 50, 0))
  I1 = result1[[3]]
  I2 = result2[[3]]
  I3 = result3[[3]]
  I_max1 = c(I_max1, max(I1, na.rm = TRUE))
  I_max2 = c(I_max2, max(I2, na.rm = TRUE))
  I_max3 = c(I_max3, max(I3, na.rm = TRUE))
}

# Compute confidence intervals and means for peak infection sizes for different vaccination numbers
result_I_max1 = t.test(I_max1, na.rm = TRUE)
confidence_interval_I_max1 = result_I_max1$conf.int
mean_I_max1 = mean(I_max1, na.rm = TRUE)
result_I_max2 = t.test(I_max2, na.rm = TRUE)
confidence_interval_I_max2 = result_I_max2$conf.int
mean_I_max2 = mean(I_max2, na.rm = TRUE)
result_I_max3 = t.test(I_max3, na.rm = TRUE)
confidence_interval_I_max3 = result_I_max3$conf.int
mean_I_max3 = mean(I_max3, na.rm = TRUE)

# Output the results
confidence_interval_I_max1
mean_I_max1
confidence_interval_I_max2
mean_I_max2
confidence_interval_I_max3
mean_I_max3


#Problem 2: Insurance claims

#Task 2a: Write a code to verify the calculation by simulating 1000 realisations of 
#the Poisson process

func2a <- function(lambda, time, number_realisations){
  number_claims <- numeric(number_realisations)
  results <- numeric(2)
  for(i in 1:1000){
    #Simulating the number of claims using Poisson process
    number_claims[i] = rpois(1, lambda*time)
  }
  
  #Probability estimation
  print(mean(number_claims>100))
  X_over_100 <- mean(number_claims>100) #P(X>100)
  
  return(X_over_100)
}

lambda <- 1.5
time <- 59
number_realisations <- 1000
X_over_100 <- func2a(lambda, time, number_realisations)

sprintf("Probability for X > 100 P(X>100) = %f", X_over_100)

#Make a figure that shows 10 realisations of X(t), 0≤t≤59, plotted in the same
#figure
t <- 0:59

#plotting
plot(t, 0*t, type="n", ylim=c(0, max(1.5*59, 100)), 
     xlab="time t", 
     ylab="Number of claims X(t)", 
     main="10 realisations of X(t)")

col_palette <- rainbow(10) #choosing the colours for the plot
for(i in 1:10){
  #cumulative sum represents total number of events up t time t
  X <- cumsum(rpois(60, lambda))
  lines(t, X, col=col_palette[i], lty = i, lwd=2)
}

#Adding a legend to the plot
legend("topright", legend=1:10, col=col_palette, lty=1:10, 
       title="Realization", cex=0.8)

#Task 2b)
#Write a code that uses 1000 simulations to estimate the probability 
#that the total claim amount exceeds 8 mill. kr. at March 1st (59 days)

func2b <- function(gamma, lambda, time, number_realisations){
  total_claims <- numeric(number_realisations)
  
  for(i in 1:number_realisations){
    #Simulating the number of claims using Poisson process
    X_t <- rpois(1, lambda *time)
    
    #Simulating the claim amounts using exponential distribution
    C_i <- rexp(X_t, rate = gamma)
    
    #Sum of all claims for each simulation
    total_claims[i] <- sum(C_i)
  }
  
  #Probability estimation
  prob_estimation <- mean(total_claims > 8)
  return(prob_estimation)
}

gamma <- 10
Z_over_8 <- func2b(gamma, lambda, time, number_realisations)
sprintf("Probability for Z > 8 P(Z>8) = %f", Z_over_8)

#Make a figure that shows 10 realisations of Z(t), 0≤t≤59, plotted in the same 
#figure

realisations <- 10

#Plotting
plot(t, 0*t, type="n", ylim=c(0, 10), 
     xlab="time t", 
     ylab="Number of claims Z(t)", 
     main="10 realisations of Z(t)")

col_palette <- rainbow(realisations) #Choosing the colours for the plot

for(i in 1:realisations){
  #Cumulative sum represents total number of events up t time t
  X <- cumsum(rpois(60, lambda))
  Z <- cumsum(rexp(X, gamma))
  
  lines(t, Z, col=col_palette[i], lty = i, lwd=2)
}

#Adding a legend to the plot
legend("topright", legend=1:10, col=col_palette, lty=1:10, 
       title="Realization", cex=0.8)

