# Loading required libraries
library(parallel); library(lubridate); library(data.table); library(mgcv); library(fitdistrplus)
library(foreach); library(readr); library(RCurl); library(rstan); library(zoo)

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

# Setting rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = 2)




m1 = stan_model(model_code = s1)

data<-list(n=nrow(d),q=length(NI),NI=as.numeric(NI),t=d$Date.symptoms)
mat<-matrix(NA,nrow=data$n,ncol=data$q)

for (i in 1:data$q) {
	for(j in 1:(NI[i]-1)){
		mat[j,i]=1
	}
}
data$z<-sum(!is.na(mat))


fit1 <- optimizing(m1,data=data,iter=5000,verbose=TRUE,init=list('A'=rep(1,data$z)))
A<-matrix(NA,nrow=data$n,ncol=data$q)

fit <- stan(model_code = s2,data=data,iter=5000,chains = 1,warmup=1000, verbose=TRUE, control = list(adapt_delta = 0.8))

param<-extract(fit,'A')[[1]]
param<-colMeans(param)

k=1
for (i in 1:data$q) {
	for(j in 1:(NI[i]-1)){
		A[j,i]=param[k]
		k=k+1
	}
}

heatmap(A, Colv=NA, Rowv=NA, labRow=d$Date, labCol= d$Date[d$Imported==0])
################
mat<-diag(0,nrow=nrow(d))

for(i in 1:nrow(mat)){
	for(j in 1:ncol(mat)){
		mat[i,j]<-data$t[j]-data$t[i]
	}
}
mat<-mat[,d$I==0]

H<-function(a,t) a*t
S<-function(a,t) exp(-a*0.5*t^2)

likelihood<-matrix(NA,nrow=nrow(A), ncol=ncol(A))

for(i in 1:nrow(likelihood)){
  for(j in 1:ncol(likelihood)){
    likelihood[i,j]<-S(A[i,j],mat[i,j])*H(A[i,j],mat[i,j])
  }
}
#normalising likelihoods
for(i in 1:ncol(likelihood)){
  likelihood[,i]=likelihood[,i]/colSums(likelihood,na.rm=TRUE)[i]}
  
  
par(mfrow=c(1,2))
plot(d$Date.symptoms,rowSums(likelihood,na.rm=TRUE),type='b',pch)

heatmap(likelihood, Colv=NA, Rowv=NA, labRow=d$Date, labCol= d$Date[d$I==0])
#normalize by likelihood of any infection occurring - row or col sums?? check literature ()  

imported<-which(d$I =="1")
tmpmat<-matrix(0,94,94-(length(NI)))
colnames(likelihood)<-NI
colnames(tmpmat)<-imported
A_sym<-cbind(likelihood,tmpmat)
A_sym <- A_sym[,paste0(sort(as.numeric(colnames(A_sym))))]


#should this be done before normalising?
rt<-(rowSums(A_sym,na.rm=T))
mean(rt)
plot(d$Date.symptoms[-NI], rt[-NI], type = "l")
lines(d$Date.symptoms[NI], rt[NI], type = "l", col = "red")
abline(h=1)
hist(rt[-NI])
hist(rt[NI])
length(d$[NI])



#should 0s be counted in mean? if not don't do this step
#A_sym[ is.nan(A_sym) ] <- 0
A[A==1.000000e+01 ]<-NA
hist(A[A<0.1 & A>0.0001])

hist(log(A))
hist(A[NI,])
hist(A[-NI,])
mean(A_sym, na.rm=T)


aNI<-median(A[NI,], na.rm=T)
aI<-median(A[-NI,], na.rm=T)

median(A[NI,], na.rm=T)

summary(A)

head(A)
likeNI<-mean(A[NI,], na.rm=T)
likI<-mean(A[-NI,], na.rm=T)

aI<- 0.01
t <-seq(0,100)
plot(t, S(likI,t)*H(likI,t), type = "l", col = "black", lwd = "3", ylab = "likelihood of transmission", xlab = "time after symptom onset")
lines(t, S(likeNI,t)*H(likeNI,t), col = "red", lwd = "3" )
legend("topright", c("median imported infector", "median autoctonous infector"), col = c("black", "red"), lwd= 3)



########################################################
library(igraph)
library(visNetwork)
heatmap(A_sym,Colv=NA, Rowv=NA)
adjgraph<-graph_from_adjacency_matrix(A_sym, weighted = T)
visnet<-toVisNetworkData(adjgraph, idToLabel = TRUE)
visnet$edges<-visnet$edges[which(visnet$edges$weight>0.3),]
visnet$nodes$imports<- d$I

visnet$nodes$shape <- "dot" 
visnet$nodes$shadow <- TRUE 
visnet$nodes$color.background<-rep("gold", length(visnet$nodes$imports))
visnet$nodes$color.background[which(visnet$nodes$imports == "1")] <- c("tomato")
visnet$nodes$color.border <- "black"
visnet$nodes$color.highlight.background <- "orange"
visnet$nodes$color.highlight.border <- "darkred"

visnet$edges$width <- 10*(visnet$edges$weight)

visNetwork(visnet$nodes, visnet$edges, main = "network ")  %>%
  visEdges(arrows = 'to')%>%
  visOptions(highlightNearest = TRUE, 
             selectedBy = "imports")


