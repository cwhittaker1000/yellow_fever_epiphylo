# Loading required libraries
library(parallel); library(lubridate); library(data.table); library(mgcv); library(fitdistrplus)
library(foreach); library(readr); library(RCurl); library(rstan); library(zoo);
library(igraph); library(visNetwork)

# Loading in dummy linelist data to illustrate the methods
urlfile <- getURL("https://raw.githubusercontent.com/IzzyRou/multitree_illustrative/54e4ceb8f4d5a5ee6dd6ac251fa754e5e4b30f76/dummy_data.csv")
raw_linelist <- read.csv(text = urlfile, stringsAsFactors = FALSE)

# Formatting and processing data
linelist <- raw_linelist[complete.cases(raw_linelist), ] # removes any rows with any missing values
linelist$Date <- as.Date(linelist$Date, format = "%d/%m/%Y") # formats date appropriately
linelist$Month <- sapply(raw_linelist$Date, function(x) as.numeric(strsplit(x, "/")[[1]][2]), USE.NAMES = FALSE) # extracts month
linelist$Date.symptoms <- decimal_date(linelist$Date) # converting to decimal date
same_dates <- duplicated(linelist$Date.symptoms)
linelist$Date.symptoms[same_dates] <- linelist$Date.symptoms[same_dates] + 0.1 # jittering(ish) around duplicates - CHECK WITH IZZY
linelist <- linelist[order(linelist$Date.symptoms), ] # ordering by date of symptoms
linelist$importation <- ifelse(linelist$Imported == "N", 0, 1)
not_imported <- which(linelist$importation == 0) # tracking which cases were imported and which weren't
linelist$rel_date_symptoms <- (linelist$Date.symptoms - linelist$Date.symptoms[1]) * 365 # calculating date of symptoms relative to first calendar case
linelist$abs_date_symptoms <- linelist$Date.symptoms # keeping track of absolute symptom date (used for plotting below)

# Generating the data required for the model
data <- list(n = nrow(linelist),              # number of cases in the linelist
             q = length(not_imported),        # number of not imported cases in the linelist
             NI = not_imported,               # index in the linelist of the not imported cases
             t = linelist$rel_date_symptoms)  # timing of symptom onsets for all cases in the linelist
mat <- matrix(NA, nrow = data$n, ncol = data$q) # matrix where each row is a potential infector, and each column is an infectee (I think - this is reversed compared to Izzy's code)
for (i in 1:data$q) {                           # assigning to 1 all of the infectors who occurred earlier in time than the infectee
  for(j in 1:(not_imported[i] - 1)){
    mat[j, i] <- 1
  }
}
data$z <- sum(!is.na(mat)) # number of entries in the A matrix (which is a vector in Sam's code)

# Setting rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

# Loading stan model
m1 <- stan_model(file = "6_Routledge_Methods_Development/netrate_model_s1.stan")
fit2 <- sampling(m1, data = data, iter = 5000, verbose = TRUE, chains = 1)

# Generating the A matrix from the Stan output
alpha_params <- extract(fit2, 'A')[[1]]
mean_alpha_param <- colMeans(alpha_params)
A <- matrix(NA, nrow = data$n, ncol = data$q)
counter <- 1
for (i in 1:data$q) {
  for(j in 1:(not_imported[i] - 1)){
    A[j, i] <- mean_alpha_param[counter]
    counter <- counter + 1
  }
}
A[is.na(A)] <- 0
mean_epsilon <- mean(extract(fit2, "epsilon_edge")[[1]])
A <- rbind(A, rep(mean_epsilon, dim(A)[2]))

## Generating matrix of time-differences between all cases
time_diff_mat <- diag(0, nrow = nrow(linelist))
for(i in 1:nrow(time_diff_mat)){
  for(j in 1:ncol(time_diff_mat)){
    time_diff_mat[i,j] <- data$t[j] - data$t[i]
  }
}
time_diff_mat <- time_diff_mat[, linelist$importation == 0] # subsetting only the non-imported cases

# Generating the Likelihood of Each Infector-Infectee Pair From Their Alphas
## Note this was implicitly calculated in the Stan code but isn't there formally atm (based on what Sam gave me)
survival <- function(alpha, time_diff) {
  return(exp(-alpha * 0.5 * time_diff^2)) 
}
hazard <- function(alpha, time_diff) {
  return(alpha * time_diff)
}
epsilon_likelihood <- exp(-mean_epsilon) * mean_epsilon ## this is uncertain to me

# Generating likelihood matrix specifying likelihood for each infector-infectee pair based on estimated A matrix
likelihood_mat <- matrix(NA, nrow = nrow(A) - 1, ncol = ncol(A))
for(i in 1:nrow(likelihood_mat)){
  for(j in 1:ncol(likelihood_mat)){
    likelihood_mat[i,j] <- survival(A[i, j], time_diff_mat[i, j], 1) * hazard(A[i, j], time_diff_mat[i, j], 1)
  }
}
likelihood_mat <- rbind(likelihood_mat, rep(epsilon_likelihood, dim(likelihood_mat)[2]))

## Normalising the Likelihoods to Get Relative Probabilities
norm_likelihood_mat <- matrix(data = NA, nrow = nrow(likelihood_mat), ncol = ncol(likelihood_mat))
for(i in 1:ncol(likelihood_mat)){
  norm_likelihood_mat[, i] <- likelihood_mat[, i]/colSums(likelihood_mat, na.rm = TRUE)[i]
}
rowSums(norm_likelihood_mat) # individual R
mean(rowSums(norm_likelihood_mat)[-length(rowSums(norm_likelihood_mat))]) # average R

# Some basic plotting and visualisation
heatmap(A, Colv=NA, Rowv=NA, labRow = linelist$Date, labCol= linelist$Date[linelist$importation == 0])
heatmap(norm_likelihood_mat, Colv=NA, Rowv=NA, labRow = linelist$Date, labCol= linelist$Date[linelist$importation == 0])
plot(linelist$Date.symptoms, rowSums(norm_likelihood_mat)[-length(rowSums(norm_likelihood_mat))], type = 'b')

# Creating a matrix that formally includes the imported cases as columns (but with all 0s to reflect fact they were
# infected elsewhere) - note this produces identical results to the above
imported <- which(linelist$importation =="1")
tmpmat <- matrix(0, 35, 34 - length(data$NI)) # 34 = number of individual in linelist
colnames(norm_likelihood_mat) <- not_imported
colnames(tmpmat) <- imported
A_sym <- cbind(norm_likelihood_mat, tmpmat)
A_sym <- A_sym[, paste0(sort(as.numeric(colnames(A_sym))))]
rt <- rowSums(A_sym, na.rm = TRUE) # remember last element will be epsilon edges so need to remove when calculating average Rt
mean(rt[-length(rt)])

# Some basic plotting of Rts for imported vs non-imported cases
plot(linelist$Date.symptoms[-not_imported], rt[-c(not_imported, length(rt))], type = "l")
lines(linelist$Date.symptoms[not_imported], rt[not_imported], type = "l", col = "red")
abline(h=1)
hist(rt[-not_imported])
hist(rt[not_imported])

# Plotting the serial interval
mean_alpha_not_imported <- mean(A[not_imported, ], na.rm = TRUE)
mean_alpha_imported <- mean(A[-not_imported, ], na.rm = TRUE)
times <- seq(0, 100)
plot(times, survival(mean_alpha_imported, times) * hazard(mean_alpha_imported, times), 
     type = "l", col = "black", lwd = "3", ylab = "likelihood of transmission", xlab = "time after symptom onset")
lines(times, survival(mean_alpha_not_imported, times) * hazard(mean_alpha_not_imported, times), col = "red", lwd = "3" )
legend("topright", c("median imported infector", "median autoctonous infector"), col = c("black", "red"), lwd= 3)

########################################################
# Scrap Code

# heatmap(A_sym,Colv=NA, Rowv=NA)
# adjgraph<-graph_from_adjacency_matrix(A_sym[-35, ], weighted = T)
# visnet<-toVisNetworkData(adjgraph, idToLabel = TRUE)
# visnet$edges<-visnet$edges[which(visnet$edges$weight>0.3),]
# visnet$nodes$imports<- linelist$I
# 
# visnet$nodes$shape <- "dot" 
# visnet$nodes$shadow <- TRUE 
# visnet$nodes$color.background<-rep("gold", length(visnet$nodes$imports))
# visnet$nodes$color.background[which(visnet$nodes$imports == "1")] <- c("tomato")
# visnet$nodes$color.border <- "black"
# visnet$nodes$color.highlight.background <- "orange"
# visnet$nodes$color.highlight.border <- "darkred"
# 
# visnet$edges$width <- 10*(visnet$edges$weight)
# 
# visNetwork(visnet$nodes, visnet$edges, main = "network ")  %>%
#   visEdges(arrows = 'to')%>%
#   visOptions(highlightNearest = TRUE, 
#              selectedBy = "imports")
# fit1 <- optimizing(m1, data = data, iter = 5000, verbose = TRUE, init = list('A' = rep(0.002, data$z)))
