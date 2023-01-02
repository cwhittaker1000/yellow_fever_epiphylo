# Load required libraries
library(odin)

# Loading SIR model
seir <- odin::odin("
                   
    ## Core equations for transitions between compartments:
    deriv(S) <- - beta * S * I / N
    deriv(E) <- beta * S * I / N - gamma * E 
    deriv(I) <- gamma * E - sigma * I
    deriv(R) <- sigma * I
    
    ## Total population size (odin will recompute this at each timestep:
    ## automatically)
    N <- S + E + I + R
    
    ## Initial states:
    initial(S) <- S_ini # will be user-defined
    initial(E) <- E_ini #will be user-defined
    initial(I) <- I_ini # will be user-defined
    initial(R) <- 0
    
    ## User defined parameters - default in parentheses:
    S_ini <- user(1000)
    E_ini <- user(1)
    I_ini <- user(1)
    beta <- user(0.2)
    gamma <- user(0.1)
    sigma <- user(0.1)
    
    incidence <- beta * S * I / N
    output(incidence) <- TRUE
                   
    ")

seir_mod <- seir$new(beta = 0.4, gamma = 0.1, sigma = 0.1, S_ini = 1000, E_ini = 0, I_ini = 1)
output <- seir_mod$run(1:150)
sir_col <- c("#8c8cd9", "#EDBF85", "#cc0044", "#83C166", "black")
matplot(output[, 1], output[, -1], xlab = "Time", ylab = "Number of individuals", type = "l", col = sir_col, lty = 1)
legend("topright", lwd = 1, col = sir_col, legend = c("S", "I", "R", "R", "incidence"), bty = "n")

# Generating data for model fitting
dur <- 125
set.seed(11)
# noisy_data <- round(c(0, diff(output[1:dur, "R"]))) # rnbinom(n = dur, mu = c(0, diff(output[1:dur, "R"])), size = 30) # low size = more overdispersion
noisy_data <- rnbinom(n = dur, mu = c(0, diff(output[1:dur, "R"])), size = 30) # low size = more overdispersion
matplot(output[1:dur, 1], c(0, diff(output[1:dur, "R"])), xlab = "Time", ylab = "Number of individuals", type = "l", col = sir_col[3], lty = 1, ylim = c(0, max(noisy_data)))
matpoints(output[1:dur, 1], noisy_data, xlab = "Time", ylab = "Number of individuals", col = sir_col[3], pch = 20)

# Generating data list
data_list <- list(N0 = 1, 
                  N2 = dur,
                  cases = round(output[1:dur, "incidence"]),
                  deaths = noisy_data,
                  EpidemicStart = 20,
                  pop = 1001)
# SI <- dexp(x = 1:data_list$N2, rate = seir_mod$contents()$gamma)
SI <- dgamma(x = 1:data_list$N2, shape = 2, scale = 1 / seir_mod$contents()$gamma)
# SI2 <- dexp(x = 1:data_list$N2, rate = seir_mod$contents()$sigma)
SI2 <- dgamma(x = 1:data_list$N2, shape = 2, scale = 1 / seir_mod$contents()$gamma)
data_list$SI <- SI
data_list$death_delay <- SI2

# try to work out what's going on with the weird infections values (negatives)
mod <- cmdstanr::cmdstan_model(stan_file = "7_Renewal_Equation_Development/renewal_equation.stan")
fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 100,
  list(
    list(mu = 3,
         kappa = 0.5,
         y = 2,
         tau = 10,
         phi = 5))
)
draws_df <- fit$draws(format = "df")
fit$summary()
hist(draws_df$kappa)
hist(draws_df$mu)
hist(draws_df$tau)
hist(draws_df$y)
hist(draws_df$phi)

# Infections Over Time
par(mfrow = c(1, 3))

infections_output <- draws_df[, grep("infections0*", colnames(draws_df))[1:dur]]
infections <- unname(colMeans(infections_output))
infections_lower <- unname(apply(infections_output, 2, quantile, 0.05))
infections_upper <- unname(apply(infections_output, 2, quantile, 0.95))

matplot(output[1:dur, 1], output[1:dur, "incidence"], xlab = "Time", ylab = "Number of infections", pch = 20, col = sir_col[3], 
        lty = 1, 
        ylim = c(0, max(c(output[1:dur, "incidence"], infections_upper[1:dur]))))
matlines(output[1:dur, 1], infections[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], infections_lower[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], infections_upper[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)

# Recoveries Over Time
recoveries_output <- draws_df[, grep("deaths0*", colnames(draws_df))[1:dur]]
recoveries <- unname(colMeans(recoveries_output))
recoveries_lower <- unname(apply(recoveries_output, 2, quantile, 0.05))
recoveries_upper <- unname(apply(recoveries_output, 2, quantile, 0.95))

matplot(output[1:dur, 1], noisy_data[1:dur], xlab = "Time", ylab = "Number of recoveries", col = sir_col[3], pch = 20,
        ylim = c(0, max(c(noisy_data[1:dur], recoveries[1:dur]))))
matlines(output[1:dur, 1], recoveries[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], recoveries_lower[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], recoveries_upper[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)

Rt <- draws_df[grep("Rt_adj0*", colnames(draws_df))]
plot(apply(Rt, 2, quantile, 0.95)[1:dur], type = "l", ylab = "Reff")
lines(apply(Rt, 2, quantile, 0.05)[1:dur], type = "l", ylab = "Reff")
lines(colMeans(Rt)[1:dur], type = "l", ylab = "Reff")

# set.seed(11)
# noisy_data <- rnbinom(n = 100, mu = output[1:100, "incidence"], size = 30) # low size = more overdispersion
# noisy_data <- round(output[1:100, "incidence"])
# matplot(output[1:100, 1], output[1:100, "incidence"], xlab = "Time", ylab = "Number of individuals", type = "l", col = sir_col[3], lty = 1, ylim = c(0, max(noisy_data)))
# matpoints(output[1:100, 1], noisy_data, xlab = "Time", ylab = "Number of individuals", col = sir_col[3], pch = 20)
# matlines(output[1:100, 1], infs, xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
# hist(rexp(n = 100000, rate = 2) + rexp(n = 1000, rate = 2))
# hist(rgamma(n = 100000, shape  = 2, scale = 1/2))
# mod <- cmdstanr::cmdstan_model(stan_file = "7_Renewal_Equation_Development/renewal_equation.stan")
# SI <- dgamma(x = 1:data_list$N2, shape = 2, scale = 1 / seir_mod$contents()$sigma)

