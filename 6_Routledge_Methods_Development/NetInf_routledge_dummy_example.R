# Loading required libraries
library(parallel); library(lubridate); library(data.table); library(mgcv); library(fitdistrplus)
library(foreach); library(readr); library(RCurl)

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

df <- list(n = nrow(linelist),              # number of cases in the linelist
           q = length(not_imported),        # number of not imported cases in the linelist
           NI = not_imported,               # index in the linelist of the not imported cases
           t = linelist$rel_date_symptoms)  # timing of symptom onsets for all cases in the linelist


# Defining the serial interval distribution
### t = time of symptom onset
### shift = minimum incubation period
### a = alpha (shaping) parameter
### epsilon = epsilon edge (unobserved infection source)

# Defining the likelihood function
likelihood_function <- function(t, shift, a, epsilon) {
  t_shift <- t - shift
  likelihood <- a * t_shift * exp(-0.5 * a * t_shift * t_shift)/epsilon
  likelihood[likelihood < 0] <- 0 # ask Izzy about this
  return(likelihood)
}

# Creating prior draws
n_prior <- 500
shift_prior <- runif(n_prior, 10, 15) # prior on the minimum incubation period
a_prior <- runif(n_prior, 0.001, 0.01) # prior on shaping parameter
e <- 1e-6 # prior on epsilon

# Running the algorithm for inference of the transmission tree with this data
ftol <- function(f,f2) (f - f2)/f # function to calculate the change in likelihood from adding additional edge

## Initially running for 1 prior sample
prior_sample_index <- round(runif(1, 1, n_prior), 0) # randomly sample an index to select from the priors
Nedges <- 5000 # maximum number of edges there could be - set to a very large number as algorithm already encourages sparsity (this should never be reached)
A <- matrix(0, nrow = length(df$NI), ncol = df$n)
R <- matrix(1, nrow = length(df$NI), ncol = df$n)
score <- rep(NA, Nedges)

delta_cji <- 0 
tol <- 1 # dummy initialisation of tol, which will be updated sequentially as the algorithm progresses

# Loops over sequentially adding in edges that increases the likelihood of the transmission network by the most
for(K in 1:Nedges){
  
  delta_ji <- matrix(0, nrow = length(df$NI), ncol = df$n)  # stores incremental likelihoods (I think); rows are potential infectees, columns are potential infectees
  k <- 1 # generic counter tracking iteration number
  
  # Looping over potential infectees, and for each, calculating the likelihood of them being infected by each potential, possible infector (possible = with calendar onset before)
  for(i in not_imported){  # this HAS to be the index in the dataframe specifying which cases haven't been imported
    
    scores <- likelihood_function(df$t[i] - df$t[1:(i - 1)], shift_prior[prior_sample_index], a_prior[prior_sample_index], e) # calculating likelihood of infection of infectee by each potential infector i.e. all infectors with symptom onset earlier in calendar time than the considered infectee
    scores <- scores * R[k, 1:(i - 1)] # multiplying by R - k is the row in R that corresponds to the infectee being considered, 1:(i-1) are all the folks who could potentially infect that infetee 			
                                       # Note that below, we assign to NA the index position in R corresponding to that infector-infectee pair, to prevent that edge being added again
                                       # scores * NA = NA and so it's just not counted for the purposes of the calculations below
    
    # No clue what's going on in this section - especially Line 78, unclear to me why delta_cji is being iteratively updated (i.e. delta_cji is both on LHS and RHS of the assignment)
    for(j in 1:(i - 1)){
      delta_cji <- delta_cji + sum(scores[-j], na.rm = TRUE)              # Calculating the score under the assumption that all possible edges to infectee OTHER than the one we want to consider has been added
      delta_ji[k,j] <- log(delta_cji + scores[j]) - log(delta_cji + 1)    # Then asking whether adding that edge in still produces an improvement in score (i.e. does adding this infector in still improve things if we consider infectee's potential connection to other infectors)
                                                                          # delta_ji is then the matrix of the increases in tree likelihood from adding edge between infector j and infectee k
                                                                          ## Note still don't fully understand why delta_cji is calculated cumulatively over all i's and all j's
                                                                          ## but Sam's quite confident it's right! :)  
    }
    k <- k + 1
  }
  greedy_index <- which(delta_ji == max(delta_ji, na.rm = TRUE), arr.ind = TRUE) # finding the infector-infectee pair (edge) that increases the score by the most
  A[greedy_index] <- delta_ji[greedy_index] # adding the best performing edge to our matrix of edges A, which stores the (unnormalised?) likelihood (probability?) of each infector-infectee pair?
  R[greedy_index] <- NA # this becomes NA because we've already assigned this edge, and so it can't be assigned again
  score[K] <- delta_cji # still unclear what this delta_cji is exactly
  
  # Evaluating the likelihood increase associated with adding the extra edge in; stops if there's little improvement  
  if (K > 1) {
    tol <- ftol(score[K], score[K-1])
    print(tol)
    if(tol < 1e-5) {
      break
    }
  }
}

# Evaluating Rt and plotting the output
As <- A/rowSums(A) # creating normalised matrix
As[is.nan(As)] <- 0 # assigning first row to 0s as initial infectee can't have been infected by anyone else in linelist
                    # nope don't the above is correct, as row 1 in the A matrix should be linelist person #2 i.e. the first non-imported case 
linelist$Rt <- colSums(As) # calculating number of infectees generated by each infector
mean(linelist$Rt)

Rt_plotting <- data.frame(year = linelist$abs_date_symptoms, Rt = linelist$Rt)
g <- gam(Rt ~ s(year),data = Rt_plotting)  
plot(Rt_plotting$year, Rt_plotting$Rt, type='b', pch = 16, col = 'red', xlab = 'Date', ylab = 'R')
lines(Rt_plotting$year, g$fitted.values, col='blue')
abline(h=1)

# I think this is actually how you calculate it
# reason Izzy has As[is.nan(As)] <- 0 is because the prob of who infected 2 was 0 for everyone, because
# it's so close in time to case 1. So epsilon edge preferred, but that's not mentioned anywhere, so 
# you just have a bunch of zeroes that return an NA when you do the calculation for As. 
# Below is what I think you have to do to actually calculate whether or not an epsilon edge is preferred
# and guarantees no rowSum is 0 (min value is e) - hence avoiding all of this issue
A_alt <- cbind(A, rep(e, dim(A)[1]))
As_alt <- A_alt/rowSums(A_alt) ## not sure if this is the right way to normalise the matrix in the way we're meant to be tbh
mean(colSums(As_alt)[-length(colSums(As_alt))])
     
match.cols<-function(val, n){
  colfunc <- colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick","darkred"), space = "rgb",bias=5)#colors
  col <- data.frame(val = seq(min(val), max(val),length.out = n), col = colfunc(n))
  out <- rep(NA, length(col))
  for(i in 1:length(val)){
    out[i] <- as.character(col[which.min(abs(col$val-val[i])), 'col'])
  }	
  return(out)
}
cols<-match.cols(As,1000)
cols[As == 0] <- 'white'
mat <- matrix(cols, nrow = nrow(As), ncol = ncol(As))

#par(mfrow=c(1,2))
x=1:ncol(As)
y=1:nrow(As)
plot(x[1],y[1],col=mat[1,1],pch=15,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),bty="n" , xaxt="n", yaxt="n",xlab='Infector',ylab='Infectee')

for(i in 1:nrow(As)){
  for(j in 1:ncol(As)){
    points(x[j],y[i],col=mat[i,j],pch=15)
  }
}
for(i in 1:nrow(As)){
  if(sum(As[i,])==0){
    for(j in 1:ncol(As)){
      points(x[j],y[i],col='yellow',pch=22)
    }
  }
}

text(x=x,y=-1 , cex=0.5,
     labels = linelist$Date , srt = 90, pos = 1, xpd = TRUE)
text(y=y,x=-2, cex=0.5,
     labels = linelist$Date[linelist$I==0],  pos = 2, xpd = TRUE,offset=-0.5)  


# perform distribution analysis
Rt<-colSums(As)
Rt[Rt==0]<-1e-4


fitg <- fitdist(Rt, "gamma")
fitg$aic
fitg <- fitdist(Rt, "exp")
fitg$aic
fitg <- fitdist(Rt, "lnorm")
fitg$aic

rbinom(1000,10,0.21)

fitg <- fitdist(Rt, "gamma")
mean(rgamma(10000,shape=fitg$estimate[1],rate=fitg$estimate[2]))
pgamma(1,fitg$estimate[1],fitg$estimate[2])

rate=fitg$estimate[2]
prob = pgamma(1,shape=fitg$estimate[1],rate=rate)

while(prob<=0.95){
  rate=rate+0.01
  prob = pgamma(1,shape=fitg$estimate[1],rate=rate)
}
mean(rgamma(10000,shape=fitg$estimate[1],rate=rate))

cbind(2016:2030,predict.gam(g,newdata=data.frame(x=2016:2030)))

##########################
## Analysing Rc by month
##########################

Rt_mat<-matrix(nrow=12,ncol=length(A_list))
for(i in 1:length(A_list)){
  A2=A_list[[i]][[1]]
  A2s<-A2/rowSums(A2)
  A2s[is.nan(A2s)]=0
  Rts<-colSums(A2s)
  Rt_mat[,i]<-as.matrix(aggregate(list(Rts),by=list(d$Month),FUN=mean))[,2]
}
avg<-apply(Rt_mat,1,mean)
q<-apply(Rt_mat,1,quantile,probs=c(0.05,0.95))

plot(1:12,avg,pch=16,ylim=c(0,2))
arrows(x, q[1,], x, q[2,], length=0.05, angle=90, code=3)
abline(h=mean(Rt_mat),col='blue')


### Parallel analyses
# A_list <- list() # list to store results
# 
# # registerDoParallel(40)
# 
# Nedges<-5000
# A_list <- foreach(xx = 1:n_prior) %dopar% {
#   A<-matrix(0,nrow=length(data$NI),ncol=data$n)
#   delta_cji<-0
#   R<-matrix(1,nrow=length(data$NI),ncol=data$n)
#   score<-rep(NA,Nedges)
#   score<-rep(NA,Nedges)
#   tol<-1
#   for(K in 1:Nedges){
#     delta_ji<-matrix(0,nrow=length(data$NI),ncol=data$n)
#     k=1
#     for(i in NI){
#       scores<-wc(data$t[i]-data$t[1:(i-1)],shift[xx],a[xx],e)
#       scores<-scores*R[k,1:(i-1)]			
#       for(j in 1:(i-1)){
#         delta_cji = delta_cji + sum(scores[-j],na.rm=TRUE)
#         delta_ji[k,j] = log(delta_cji + scores[j]) - log(delta_cji + 1)
#       }
#       k=k+1
#     }
#     greedy_index<-which(delta_ji == max(delta_ji,na.rm=TRUE), arr.ind = TRUE)
#     A[greedy_index]<-delta_ji[greedy_index]
#     R[greedy_index]<-NA
#     score[K]=delta_cji
#     if(K!=1) tol<-ftol(score[K],score[K-1]); print(tol)
#     if(tol<1e-5) break
#   }
#   return(list(A,score))
# }
# 
# ########################
# #posterior mean of matrix
# #########################
# 
# score<-list()
# A<-matrix(0,nrow=length(data$NI),ncol=data$n)
# for(i in 1:length(A_list)){
#   A=A+A_list[[i]][[1]]
#   score[[i]]<-A_list[[i]][[2]]
# }
# A=A/length(A_list)
# 
# A_avg=A
# #score analysis
# index<-450
# plot(1:index,score[[1]][1:index],pch=16,cex=0.5,type='l',ylim=c(min(unlist(score),na.rm=TRUE),max(unlist(score),na.rm=TRUE)))
# for(i in 1:length(A_list)){
#   lines(score[[i]][1:index],pch=16,col=i,cex=0.5)
# }
# 
# score
