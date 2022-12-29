#############################
# Spatial Netrate
#'
#' @param tmat timeseries matrix
#' @param dmat distance matrix
#' @param DataType distance data type - euclidian, friction, binary
#' @param SpatialKernel
#' @param fixed vector of fixed parameters
#' @param alpha prior mean and variance for serial interval parameter alpha, if fixed the fixed value
#' @param delta prior mean and variance for distance function parameter delta, if fixed the fixed value
#' @param epsilon prior mean and variance for epsilon edge, if fixed the fixed value
#'
#' @return
#' @export
#'
#' @examples




spatialnetrate <-
  function(tp = tmat,
           dp = dmat,
           DataType = euclidian,
           SpatialKernel = "nil",
           fixed = "epsilon",
           alpha = c(0.003, 0.01),
           delta = c(0.01, 0.01),
           epsilon = 1e-3) {
    
    library(tensorflow)
    library(reticulate)
    tfp <- reticulate::import("tensorflow_probability",convert=FALSE)
    
    ###############################
    #key small functions to define
    ##############################
    
    
    f1<-function(t,d,a,r) {
      likelihood = a*t*exp(-0.5*a*t*t) * (1/d)^r
      
      # survival (exp(-0.5*a*t*t) * (1/d)^(r-1) )/(r-1)
      # hazard (a*(r-1)*t)/(d)
      
      
      return(likelihood)
    }
    f2<-function(t,d,a,r) {
      likelihood = a*t*exp(-0.5*a*t*t) * exp(-d^2*r)
      
      # survival (exp(-0.5*a*t*t) * (1/d)^(r-1) )/(r-1)
      # hazard (a*(r-1)*t)/(d)	
      
      return(likelihood)
    }
    f3<-function(t,d,a,r) {
      likelihood = a*t*exp(-0.5*a*t*t) * exp(-abs(d)*r)
      
      # survival -exp(-0.5*a*t*t - r*d)
      # hazard 	
      
      
      return(likelihood)
    }
    match.cols<-function(val,n){
      colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"))
      col<-data.frame(val=seq(min(val),max(val),length.out=n),col=colfunc(n))
      out<-rep(NA,length(col))
      for(i in 1:length(val)){
        out[i]<-as.character(col[which.min(abs(col$val-val[i])),'col'])
      }	
      return(out)
    }
    
    
    ###################################################
    ## Data defining section
    ###################################################
    
    #set.seed(1234)
    tf$reset_default_graph()
    N = ncol(tmat) # number of infectors
    Q = nrow(tmat) # number of infectees
    
    tp = tf$placeholder(shape = shape(Q, N), dtype = tf$float32)
    dp = tf$placeholder(shape = shape(Q, N), dtype = tf$float32)
    
    ## A defining section
    if (fixed == "alpha") {
      A <- alpha
      #this may need more coding, maybe just A=alpha is okay too?
      A <- tf$constant(alpha, shape = shape(Q, N), dtype = tf$float32)
    } else{
      A <-
        tf$nn$relu(tf$Variable(
          tf$constant(alpha[1], shape = shape(Q, N), dtype = tf$float32),
          name = 'alpha'
        )) # length scale for space
      A <-
        tf$minimum(A, tf$constant(0.01, shape = shape(Q, N), dtype = tf$float32))
      A <- A * tf$sign(tp)
      Ad <-
        tfp$distributions$TruncatedNormal(
          loc = alpha[1],
          scale = alpha[2],
          high = 0.01,
          low = 0.0
        ) #0.005,0.001
      
    }
    
    ## D defining section
    if (fixed == "delta") {
      D <- delta
      
    } else{
      D <-
        tf$nn$relu(tf$Variable(
          tf$constant(delta[1], shape = shape(1L), dtype = tf$float32),
          name = 'distance_scale'
        )) # length scale for space
      #Dd<-tfp$distributions$Exponential(0.001)
      Dd <-
        tfp$distributions$TruncatedNormal(
          loc = delta[1],
          scale = delta[2],
          low = 0.0,
          high = 0.5
        ) #0.005,0.001
      
      
    }
    
    
    #E defining setcion
    if (fixed == "epsilon") {
      E<-tf$nn$relu(tf$Variable(tf$constant(epsilon,shape=shape(Q),dtype=tf$float32),name='epsilon_edge'))
      #E<-epsilon
    } else{
      E <-
        tf$nn$relu(tf$Variable(
          tf$constant(epsilon[1], shape = shape(Q), dtype = tf$float32),
          name = 'epsilon_edge'
        )) # length scale for space
      Ed <-
        tfp$distributions$TruncatedNormal(epsilon[1],
                                          scale = epsilon[2],
                                          high = 0.9,
                                          low = 0.0)
    }
    if(fixed == "nil"){
      print("no fixed parameters - identifiability could be an issue")
    }
    ###################################################
    ## Spatial kernel section
    ###################################################
    
    
    if (SpatialKernel == "exponential") {
      H <- A * tp * tf$exp(-dp * D)
      
    }else if(SpatialKernel == "Power Law"){
      H <- A * tp * ((1/dp)^D)
      
    } else if (SpatialKernel == "linear") {
      H <- A * tp * D * dp
      
      
    } else if (SpatialKernel == "gaussian"){
      H <- A * tp * (1/(2*3.14*(D*D)))*exp(-dp)/(2*D*D)
      
      
      
    } else {
      H <- A * tp
    }
    
    
    #######################################################
    ##tensorflow code
    #######################################################
    nll1 <-
      tf$reduce_sum(log(tf$reduce_sum(H, 1L) + E)) #first loop and second loop
    S <- tf$negative(0.5 * A * tp * tp)
    
    nll2 <-
      tf$reduce_sum(tf$reduce_sum(S, 1L) - E) #first loop and second loop
    log_prior <-
      tf$negative(tf$reduce_sum(Ad$log_prob(A)) + tf$reduce_sum(Dd$log_prob(D)))
    if(fixed != "epsilon"){
      log_prior <-
        tf$negative(tf$reduce_sum(Ad$log_prob(A)) + tf$reduce_sum(Dd$log_prob(D))+tf$reduce_sum(Ed$log_prob(E)))
    }else{ 
      log_prior <-
        tf$negative(tf$reduce_sum(Ad$log_prob(A)) + tf$reduce_sum(Dd$log_prob(D)))}
    #log_prior<- tf$negative(tf$reduce_sum(Ad$log_prob(A)))
    
    nll <- tf$negative(nll1 + nll2)
    log_posterior <- tf$add(nll, log_prior)
    
    F <-
      H * tf$exp(S) # transmission likelihood is hazard times survival
    
    train_step = tf$train$AdamOptimizer(1e-3)$minimize(log_posterior) # gradient step size !! IMPORTANT FREE PARAMETER !! 
    init = tf$global_variables_initializer()
    myconfig = tf$ConfigProto(log_device_placement = FALSE)
    sess = tf$Session(config = myconfig)
    steps = 10000
    sess$run(init) # initialise session
    nll_store <- 1:steps
    
    #############################################################
    #BFGS optimising
    #############################################################
    # 
    # optimizer = tf$contrib$opt$ScipyOptimizerInterface(
    #   log_posterior,
    #   method = 'L-BFGS-B',
    #   options = dict('maxiter' = 2000000L,
    #                  'disp' = TRUE)
    # )
    # sess = tf$Session()
    # init = tf$global_variables_initializer()
    # sess$run(init) # initialise session
    # 
    # optimizer$minimize(session = sess,
    #                    feed_dict = dict(tp = tmat,
    #                                     dp = dmat))
    # 
    # 
    # #sess$run(E)
    # Ds<-sess$run(D)
    # value <- sess$run(list(nll, log_prior, log_posterior),
    #                   feed_dict = dict(tp = tmat,
    #                                    dp = dmat))
    # print(paste(value[[1]], value[[2]], value[[3]]))
    # 
    
    #manual opt
    for (i in 1:steps) {
      value<-sess$run(list(train_step,nll,log_prior,log_posterior),
                      feed_dict=dict(
                        tp=tmat,
                        dp=dmat
                      ))
      nll_store[i]=value[[2]]
      if(i %% 1000 == 0){
        print(paste(value[[2]],value[[3]],value[[4]]))
        #			print(paste(sess$run(E)))
        
      }
    }
    # plot(nll_store)
    
    #############################
    
    As <- sess$run(list(A), feed_dict = dict(tp = tmat,
                                             dp = dmat))[[1]]
    meanA<-mean(As[As != 0])
    par(mfrow = c(1, 4))
    Fs <- sess$run(list(F), feed_dict = dict(tp = tmat,
                                             dp = dmat))[[1]]
    Ds<-sess$run(D)
    eps_edge<-sess$run(E)
    #eps_edge <- rep(E, nrow(Fs))
    
    Fs <- cbind(Fs, eps_edge)
    
    Fs = sweep(Fs, 1, rowSums(Fs), FUN = '/')
    Rt <- colSums(Fs, na.rm=TRUE)[1:(ncol(Fs) - 1)]
    
    #plotting space-time likelihood and kernels
    x=seq(15,165,length.out=100)
    d=seq(1,200,length.out=100)
    z=expand.grid(x,d)
    r= as.numeric(Ds)
    # par(mfrow=c(1,3))
    # 
    # if (SpatialKernel == "exponential") {
    #   
    #   y=f3(z[,1]-15,z[,2],meanA,r)
    #   
    #   plot(z,col=match.cols(y,1000),pch=15)
    #   plot(d,exp(-d*r),type='b')
    #   
    # }else if(SpatialKernel == "Power Law"){
    #   
    #   y=f1(z[,1]-15,z[,2],meanA,r)
    #   plot(z,col=match.cols(y,1000),pch=15)
    #   plot(d,(1/d)^r,type='b')
    #   
    # } else if (SpatialKernel == "linear") {
    #   H <- A * tp * D * dp
    #   
    #   
    # } else if (SpatialKernel == "gaussian"){
    #   y=f2(z[,1]-15,z[,2],meanA,r)
    #   plot(z,col=match.cols(y,1000),pch=15)
    #   plot(d,exp(-d^2*r),type='b')
    #   
    #   
    #   
    # } else {
    #   
    # }
    # 
    # 
    # plot(dd$t,Rt)
    # abline(h=mean(Rt))
    # 
    
    
    return(list(Rt,Ds,As,eps_edge))
    
  }

