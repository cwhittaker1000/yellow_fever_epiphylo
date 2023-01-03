
# Simulates a negative binomial branching process whilst accounting for susceptible depletion
## NOTE: Taken from bp models (Seb Funk and Flavio Finger) - currently this function
##       simulates all infections (and their times) from a given infection, then does the same
##       for the next soonest infection. I.e. if one person is infected at t=5, they'll generate
##       e.g. 3 infections at t>5. All those 3 infections are now unavailable for infection by folks
##       infected (and hence infectious) at t=6, t=7 etc. Not clear to me whether this is a problem - seems
##       like it's negating competing risks by effectively giving precedence to the earlier infections 
##       - maybe not a problem, but probably worth emailing Seb Funk. 
## NOTE: Think the alternative would be to have at like every timestep (?) or something do a draw 
##       from a multinomial with infection by a specific individual proportion to their respective serial intervals
##       or something. Not clear how you'd take variation in R into account though. 
chain_sim_susc_nb <- function(
    mn_offspring,                  # Average number of secondary cases for each case
    disp_offspring,                # Dispersion coefficient (var/mean) of the number of the number of secondary cases
    serial,                        # The serial interval. A function that takes one parameter n, which is the number of serial intervals to randomly sample.
    tf,                            # End time
    initial_infections,            # Initial number of infections to seed with
    pop) {                         # Population size
  
  if (disp_offspring <= 1) { 
    stop("Offspring distribution 'nbinom' requires argument disp_offspring > 1.")
  }
    
  offspring_fun <- function(n, susc) {
    new_mn <- mn_offspring * susc / pop # Allow for susceptible depletion
    size <- new_mn / (disp_offspring - 1)
    
    ## Use a right truncated nbinom distribution to avoid more infections than susceptibles
    truncdist::rtrunc(
      n,
      spec = "nbinom",
      b = susc,
      mu = new_mn,
      size = size)
  }
  
  ## Dataframe containing infected individuals, initialised with first infected cases
  tdf <- data.frame(
    id = seq_len(initial_infections),        ## Id for those individuals
    ancestor = NA_integer_,                  ## Who infected them
    generation = 1L,                         ## Which generation of infections it occurred in
    time = 0,                                ## Time of initial infection
    offspring_generated = FALSE)             ## Tracks whether that individual has generated offspring
  
  # Continue simulation until no unsimulated has t <= tf OR there are no susceptibles left
  susc <- pop - initial_infections  ## Number of susceptibles remaining
  t <- 0                            ## Initial time
  while (any(tdf$time[!tdf$offspring_generated] <= tf) & susc > 0) {
    
    ## Selecting the case to generate offspring
    t <- min(tdf$time[!tdf$offspring_generated]) # Lowest unsimulated t in those individuals who have not generated offspring yet
    
    ## Getting the index of the first individual in the dataframe with t, extracting the relevant variables
    idx <- which(tdf$time == t & !tdf$offspring_generated)[1]
    id_parent <- tdf$id[idx]
    t_parent <- tdf$time[idx]
    gen_parent <- tdf$generation[idx]
    
    ## Generate number of offspring for the case generating offspring this step
    n_offspring <- offspring_fun(1, susc)
    current_max_id <- max(tdf$id) # calculating current max id so that new offspring can be added with ids > this
    
    ## Mark that individual as having generated offspring
    tdf$offspring_generated[idx] <- TRUE
    
    ## Adding newly infected individuals to the dataframe
    if (n_offspring > 0) {
      
      new_times <- serial(n_offspring) # Draw serial intervals for each newly infected individual.
                                       # See note above - unclear to me whether this way of doing it is a problem and negates potential competing risks
      if (any(new_times < 0)) {
        stop("Serial interval must be >= 0.")
      }
      
      new_df <- data.frame(
        id = current_max_id + seq_len(n_offspring),
        ancestor = id_parent,
        generation = gen_parent + 1L,
        time = new_times + t_parent,
        offspring_generated = FALSE)
      
      ## Add new cases to the main dataframe
      tdf <- dplyr::bind_rows(tdf, new_df)
    }
    
    ## Adjust number of susceptible individuals
    susc <- susc - n_offspring
  }
  
  ## Sort output and remove columns not needed
  tdf <- tdf[order(tdf$time, tdf$id), ]
  tdf$offspring_generated <- NULL
  
  return(tdf)
}

library(dplyr); library(tictoc)

serial_fn <- stats::rgamma
hist(rgamma(10000, shape = 5, rate = 1))
set.seed(1)
x <- chain_sim_susc_nb(mn_offspring = 3, disp_offspring = 2, serial = function(x) { serial_fn(x, 5, 1) }, 
                       initial_infections = 2, tf = 100, pop = 100)

tic()
iter <- 25
initials <- 5
R0 <- seq(1.5, 3.5, 0.5)
final_size_df <- array(data = NA, dim = c(iter, initials, length(R0)))
for (i in 1:iter) {
  set.seed(i)
  for (j in 1:5) {
    for (k in 1:length(R0)) {
      x <- chain_sim_susc_nb(mn_offspring = R0[k], disp_offspring = 10, serial = function(x) { serial_fn(x, 5, 1) }, 
                             initial_infections = j, tf = 100, pop = 87)
      final_size_df[i, j, k] <- max(x$id) 
    }
  }
}
toc()

apply(final_size_df[, 5, ], 2, mean)

apply(final_size_df[, 5, ], 2, min)
apply(final_size_df[, 5, ], 2, max)

apply(final_size_df, 2, mean)

mean(final_size_df$init_one)
mean(final_size_df$init_five)

final_size_df$init_five
