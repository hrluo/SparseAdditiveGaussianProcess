##############
# SAGP Model #
##############

message('sagp.R module loaded!\n')

# sagp
# Y observations, n*1 vector.
# X observed locations, n*d matrix (only d=1 for now). 
# pinfo list of prior information, used for specifying the hyper-parameters.
# minfo list of MCMC information. number of MCMC steps, burn-in, etc.
# tinfo list of partition design that defines the domains of the components and the 
#   structure of the layers. 
#   E.g., 3 layers, 7 components in bipartite component design (1 comp in 1st layer, 2 comp in 2nd layer, 4 comp in 3rd layer) 
#     is represented by:
#     tinfo = list(data.frame(loc=c(0.5),
#                             rad=c(0.5)),
#                  data.frame(loc=c(0.25,0.75),
#                             rad=c(rep(0.25,2))),
#                  data.frame(loc=c(0.125,.375,.625,0.875),
#                             rad=c(rep(0.125,4))))
# Xnew is the locations over which we want to conduct our predictions.
# sampleEta/sampleRho are booleans, should eta/rho be sampled in MCMC algorithm.
# storePredX boolean, should the predictions at the observed locatios be stored? If X large, 
#   should be set to FALSE and only predictions at Xnew are stored.
sagp<-function(Y,X,pinfo,minfo,tinfo,Xnew=X,sampleEta=T,sampleRho=F,storePredX=T)
{
  #Seed, for reproducibility
  set.seed(minfo$seed)
  
  #Transform Y to have mean = 0 and sd = 1
  meanY <- mean(Y)
  sdY <- sd(Y)
  Y <- (Y-meanY)/sdY
  #At the end of the code, samples of Y are transformed back to original scale
  
  #overall number of observations in the sample.
  n =nrow(X)
  if(n!=length(Y))stop('sagpfun::sagp @ the length of X and Y shall match.')
  #overall number of new locations to be predicted.
  nNew = nrow(Xnew)
  #number of additive components - defined by the list with tree partition structure
  #N=length(unlist(tinfo))/2
  #number of pseudo inputs - if equal for all components
  m=pinfo$m
  #number of layers
  L=length(tinfo)
  #dimension of X
  D=ncol(X)
  
  #number of samples to store and return in the result list
  Nmcmc    = minfo$Nmcmc 
  #number of burn-in iterations
  burnIn   = minfo$burnIn 
  #total number of iterations to run
  Nsamples = Nmcmc+burnIn
  
  #Identify the locations (X and Xnew) corresponding to each component, for each layer
  # indexCompX is a list of L elements. For each layer l, indexCompX[[l]] is a vector with length n.
  #   For each element X[i] in X[1],...,X[n], indexCompX[[l]][i] is the index of the component in the l-th layer to which X[i] belongs.
  # indexCompNew is a list of L elements. For each layer l, indexCompX[[l]] is a vector with length nNew.
  #   For each element X[i] in X[1],...,X[nNew], indexCompX[[l]][i] is the index of the component in the l-th layer to which X[i] belongs.
  # layers is a data.frame that links layers and components
  #   e.g. layer=data.frame(comp=c(1,2,3),layer=c(1,2,2)) means the component 1 is in the 1st layer while the component 2,3 are in the 2nd layer.
  #   It is used throughtout the code to know which is the layer of each component: layers$layer[layers$comp==j] is the layer of component j.
  indexCompX=vector("list",L)
  indexCompNew=vector("list",L)
  layers <- NULL
  # m_l is a vector of length N=N_1+N_2+...+N_L recording how many PIs should there be in each componet, or the l-th component to be more precise.
  m_comp    <- c()
  #Initialize the layers dataframe. numCompTot is the iterator that records the loop;
  numCompTot <- 0
  for(l in 1:L) {
    #The outer loop loops through 1,2,...,L layer.
    
    #tinfo_layer is a data frame consisting of N_l*2 elements, 
    #   the number of rows N_l means the number of components in the l-th layer
    #   the number of PIs in each component m are determined before.
    #   the first column is recording the center of each component.
    #   the second column is recording the radius of each component.
    tinfo_layer       <- tinfo[[l]]
    indexCompX[[l]]   <- rep(NA, n)
    indexCompNew[[l]] <- rep(NA, nNew)
    #N_l is the number of components in the l-th layer.
    N_l <- nrow(tinfo_layer$loc)
    for(j_l in 1:N_l) {
      #The inner loop loops through different components j_l=1,2,...,nrow(tinfo_layer) in each layer/ the l-th layer.
      
      #loc_j_l is the 'center' of the j_l component domain
      #rad_j_l is the 'radius' of the j_l component domain.
      loc_j_l <- tinfo_layer$loc[j_l,] 
      rad_j_l <- tinfo_layer$rad[j_l,]
      
      ids_j_l <- ids_X_comp(X,loc=loc_j_l,rad=rad_j_l)
      cat('Initialize...Layer',l,'Component',j_l,'/',N_l,'\n');message('loc=',loc_j_l,' rad=',rad_j_l,'\n')
      #Check:are there points in this component that have not been already assigned to other components?
      #Need to check that the points are not already assigned because points could be on the border between
      #two components and could be assigned to both.
      idToReplace <- ids_j_l[is.na(indexCompX[[l]][ids_j_l])]
      #If yes: include it. Otherwise: do not increase the counter and do not include the component.
      if(length(idToReplace) > 0) {
        numCompTot <- numCompTot + 1
        
        indexCompX[[l]][idToReplace] <- numCompTot
        
        ids_new_j_l <- ids_X_comp(Xnew,loc=loc_j_l,rad=rad_j_l)
        indexCompNew[[l]][ids_new_j_l] <- numCompTot
        
        layers <- rbind(layers, data.frame(layer=l,comp=numCompTot))
        
        # Assign different m_comp (number of PIs) in this case we use the same number of PIs in each component.
        #m_comp <- c(m_comp,m-4*floor(log(numCompTot,base=sqrt(2))))
        # or set  all m_comp to be the same
        m_comp<-c(m_comp,m)
      } 
      
      
    }
  }
  #Note: numCompTot is the number of temporary real components, where I have at least one data point.
  
  #cat(m_comp)
  m_max <- max(m_comp)
  m_compActual<-m_comp
  if(m_max<=0)stop('sagpfuns::sagp @ the maximal number of PIs in components must be greater than 0.')
  
  #"Deactivate" components corresponding to areas where there are not enough data points.
  # Loop from last component (lowest layer, local) to first (higher layer, global).
  # For each component j, look how many components are nested in j and compute necessary number of PIs.
  # If number of data points is not sufficient, deactivate the most local components nested in j (lowest layer)
  # and check if there are enough data points. Iterate, removing one layer of components each time. 
  for(j in numCompTot:1) {
    layer_j <- layers$layer[layers$comp==j] #Layer of component j
    ids_layer_j <- (indexCompX[[layer_j]] %in% j) #ids of points in comp j
    nestedComp <- list() #list where store the nested components
    layers_below_j <- 0 #counter of nested layers
    for(l in layer_j:L) {
      #Which components are nested in j in layer l? stored in nestedComp_l
      compOfNestedPoints <- indexCompX[[l]][ids_layer_j] 
      compOfNestedPoints_noNA <- compOfNestedPoints[!is.na(compOfNestedPoints)]
      nestedComp_l <- sort(unique(compOfNestedPoints_noNA))
      #If there is one component at layer l: store it in list nestedComp.
      if(length(nestedComp_l)>0) {
        layers_below_j <- layers_below_j + 1
        nestedComp[[layers_below_j]] <-  nestedComp_l 
      } else {
        break
      }
    } 
    #Given all the components nested in j, how many PIs do we need?
    req_m <- length(unlist(nestedComp)) * m
    #How many data points are available in component j?
    numPointsComp <- sum(ids_layer_j)
    #If the component does not contain enough points for its PI and the PI of the nested components:
    if(numPointsComp<req_m) {
      
      #Deactivate components from the layer at lowest level (most local) to the highest (more global).
      for(l_deact in layers_below_j:1) {
        #Components to deactivate
        compToDeact <- nestedComp[[l_deact]]
        #Loop to deactivate the components. Also, drop the component from the dataframe layers.
        for(j_deact in compToDeact) {
          layer_j_deact <- layers$layer[layers$comp==j_deact]
          indexCompX[[layer_j_deact]][indexCompX[[layer_j_deact]] %in% j_deact] <- NA
          indexCompNew[[layer_j_deact]][indexCompNew[[layer_j_deact]] %in% j_deact] <- NA
          layers <- layers[layers$comp!=j_deact,]
        }
        #When a layer of components is removed, check: do we have enough data points now? 
        nestedComp <- nestedComp[1:(l_deact-1)]
        req_m <- length(unlist(nestedComp)) * m #Required data points now
        #If now there are enough data points: OK, move to check another component
        if(numPointsComp>=req_m) {
          break
        }
      }
      
    }
  }
  
  #Rename remaining components from 1 to the number of remaining components and match the number of PIs
  #assigned to the "old" components to the "new" components.
  #Assign "old" labels to comp_old and "new" labels to comp.
  layers$comp_old <- layers$comp 
  layers$comp <- 1:nrow(layers)
  #Same for m_comp
  m_comp_old <- m_comp
  m_comp <- rep(NA,nrow(layers))
  #Loop on the "new" labels of the components
  for(j in layers$comp) {
    layer_j <- layers$layer[layers$comp==j]
    indexCompX[[layer_j]][indexCompX[[layer_j]] %in% layers$comp_old[layers$comp==j]] <- j
    indexCompNew[[layer_j]][indexCompNew[[layer_j]] %in% layers$comp_old[layers$comp==j]] <- j
    m_comp[j] <- m_comp_old[layers$comp_old[layers$comp==j]]
  }
  #Remove "old" labels
  layers$comp_old <- NULL
  m_compActual<-m_comp
  
  #If a layer does not have any components left, remove layer
  for(l in L:1) {
    if(all(is.na(indexCompX[[l]]))) {
      indexCompX <- indexCompX[1:(l-1)]
      indexCompNew <- indexCompNew[1:(l-1)]
    }
  }
  
  #Rename objects with correct number of components:
  N <- nrow(layers)
  L <- length(indexCompX)
  
  
  #Define prior on rhos and eta using the correct number of components and layers
  #Same priors for all the dimensions
  pinfo$rho=list() # one rho prior per layer, list of values defined below
  pinfo$eta=list() # one eta prior per layer, list of values defined below
  #Fixed vector of rhos
  fixedRhos <- 10^(seq(pinfo$logRhoFirstLayer, pinfo$logRhoLastLayer,length = L))
  for(l in 1:L) {
    #rho: fixed values, equispaced in the log_10-scale between 10^-1 and 10^-50
    pinfo$rho[[l]] <-  fixedRhos[l]
    #eta: gamma prior with mean 1/ ((1-r) * r^(l-1)) in layer l
    #pinfo$eta[[l]] <- list(a=1, b = (1-pinfo$r) * (pinfo$r^(l-1)))
    pinfo$eta[[l]] <- list(a=50, b = 50*(1-pinfo$r) * (pinfo$r^(l-1)))
  }
  
  #Rhos
  #Initialize array to store rho samples, each row stores one MCMC sample for rho.
  rho=array(NA,dim=c(Nsamples,N,D))
  #First value for each component: sampled or set to fixed value 
  for(j in 1:N) {
    for(d in 1:D) {
      layer_j <- layers$layer[layers$comp==j]
      if(sampleRho) {
        rho[1,j,d]=rbeta(1,pinfo$rho[[layer_j]]$alpha, pinfo$rho[[layer_j]]$beta)    
      } else {
        rho[1,j,d]=pinfo$rho[[layer_j]]
      }
    }
  }
  #Semiwidth of MH proposal for rho (same for all component and for each dimension)
  semiWidthRho=matrix(minfo$semiWidthRho,nrow=N,ncol=D)
  
  #Etas
  #Initialize matrix to store eta samples, each row stores one MCMC sample for eta.
  eta=matrix(NA,nrow=Nsamples,ncol=N)
  #print(pinfo$eta)
  #First value for each component: sampled
  for(j in 1:N) {
    layer_j <- layers$layer[layers$comp==j]
    #if(sampleEta) {
    eta[1,j]=rgamma(1,pinfo$eta[[layer_j]]$a, pinfo$eta[[layer_j]]$b)  
    #}
  }
  #Semiwidth of MH proposal for eta(same for all component)
  semiWidthEta=rep(minfo$semiWidthEta,N)
  
  #lambda - reciprocal of overall Gaussian noise variance
  lambda=rep(NA,Nsamples)
  lambda[1]=rgamma(1,pinfo$a,pinfo$b)
  
  #D_mcmc is a matrix that stores the last set of residuals for each component, each column stores residuals for each component.
  D_mcmc=matrix(0,nrow=n,ncol=N)
  for(j in 1:N){
    D_mcmc[,j]=Y
  }
  
  #xBar/fBar are lists storing the samples of pseudo-inputs/targets (for *ALL* samples, also burn-in iterations).
  xBar=vector("list",Nsamples)
  fBar=vector("list",Nsamples)
  for(i in 1:Nsamples){
    #Initialize the xBar to be zero vectors.
    #Each term in xBar is a matrix representing the xBars. Each column of the matrix is a vector for xBar in an additive component.
    xBar[[i]]=array(0,dim=c(m_max,N,D))
    #Initializes the fBar to be zeros.
    #Each term in xBar is a matrix representing the fBars. Each column of the matrix is a vector for fBar in an additive component.
    fBar[[i]]=matrix(0,nrow=m_max,ncol=N)
    
  }
  complete_vec<-function(vec,len){
    ret<-c(vec,rep(NA,len-length(vec)))
    return(ret)
  }
  complete_mat<-function(mat,len){
    ret<-rbind(mat,matrix(NA,nrow=(len-nrow(mat)),ncol=ncol(mat)))
    return(ret)
  }
  #Sample first xBar/fBar we must complete the vector to store them in a mtrix form if m_comp differ.
  for(j in 1:N) {
    idSampleX <- sample.int(n,m_comp[j])
    xBar[[1]][,j,] <- complete_mat(matrix(X[idSampleX,],ncol=D),m_max)
    fBar[[1]][,j] <- complete_vec(mvtnorm::rmvnorm(1,mean=rep(0,m_comp[j]),sigma=diag(1,m_comp[j]) ),m_max)
  }
  
  #Matrices storing the means of each component in last MCMC iteration (at observed and new locations). 
  #Cannot initialize as non-0.
  mu.component.temp=matrix(0,nrow=n,ncol=N)
  mu.new.component.temp=matrix(0,nrow=nNew,ncol=N)
  
  #Lists to store samples of the means of each component at each MCMC iteration - only after burnin to reduce memory allocation.
  mu.component=vector("list",Nmcmc)
  mu.new.component=vector("list",Nmcmc)
     
  #Lambda_diag is a list with N components, Lambda_diag[[j]] is the diagonal of Lambda_j (diagonal matrices)
  Lambda_diag=vector("list",N)
  for(j in 1:N) {
    layer_j <- layers$layer[layers$comp==j]
    #The number of elements of Lambda_diag[[j]] is the number of data points X belonging to component j
    Lambda_diag[[j]] <- rep(1, sum(indexCompX[[layer_j]] %in% j))
  }
  
  #Vector of in-sample and out-sample means (only last iteration)
  mu.temp = matrix(0,nrow=n,ncol=1)
  mu.new.temp = matrix(0,nrow=nNew,ncol=1)
  
  #Stored in-sample and out-sample means - only stored after burn-in to reduce memory allocation
  # Means at observed locations are stored only if storePredX==T
  if(storePredX==T) {
    mu=matrix(NA,nrow=n,ncol=Nmcmc)
    mu[,1]=mu.temp
  } else {
    mu <- NULL
  }
  mu.new=matrix(NA,nrow=nNew,ncol=Nmcmc)
  mu.new[,1]=mu.new.temp
  
  #Stored in-sample and out-sample predictions - only stored after burn-in to reduce memory allocation
  # Prediction at observed locations are stored only if storePredX==T
  # No need to initialize them
  if(storePredX==T) {
    Yhat=matrix(NA,nrow=n,ncol=Nmcmc)
  } else {
    Yhat <- NULL
  }
  Yhat.new=matrix(NA,nrow=nNew,ncol=Nmcmc)
  
  #Acceptance rates of MH samples
  acceptCountRhos <- matrix(0,nrow=N,ncol=D)
  sampledRhos     <- matrix(0,nrow=N,ncol=D)
  acceptCountEtas <- rep(0,N)
  sampledEtas     <- rep(0,N)
  
  #Say some stuff
  cat("Sparse Additive Gaussian Process (SAGP) Model\n")
  cat("Number of observations: ",n,"\n")
  cat("Number of predictive locations: ",nrow(Xnew),"\n")
  
  cat("Number of layers: ",L,"\n")
  cat("Number of additive components: ",N,"\n")
  cat("Number of inducing variables measure at pseudo-inputs per component: ",m,"\n")
  cat("\n")
  cat("Running MCMC for ",Nsamples," iterations (burn-in:",burnIn,", saved:", Nmcmc, ")\n")
  
  pb <- txtProgressBar(min=2,max=Nsamples,style=3)
  # Main MCMC loop
  for(i in 2:Nsamples) {
    
    setTxtProgressBar(pb, i)
    
    #Initialize i-th step to (i-1)-th step for all the sampled parameters - so we don't need to worry about i or i-1
    xBar[[i]] <- xBar[[i-1]]
    fBar[[i]] <- fBar[[i-1]] 
    rho[i,,]   <- rho[i-1,,]
    eta[i,]   <- eta[i-1,]
    lambda[i] <- lambda[i-1]
    
    # BOTTOM-UP
    #Arrange the pseudo-inputs for the current iterations of the MCMC loop
    # Pseudo-inputs are assigned from the last (finest) to the first (roughest) layer, to make sure 
    # that there are enough pseudo-inputs in each component.
    # Pseudo-inputs are assigned to each component *WITHOUT* repetition.
    Xavail <- rep(T, n) #Available X (all at beginning)
    listIndexEligX <- list() #Store eligible X for each component
    #From last to first component:
    for(j in N:1) { 
      layer_j <- layers$layer[layers$comp==j]
      #Eligible X: which are X that happen to fall in layer_j layer.
      #Available X: which are X that are not used in all previous components.
      listIndexEligX[[j]] <- which(indexCompX[[layer_j]]==j)
      #Eligible AND available X: 
      indexesEligAndAvailX <-  which(indexCompX[[layer_j]]==j & Xavail) 
      #Sample from eligible AND available X:
      indexesSampledX <- sample(indexesEligAndAvailX, m_comp[j])
      #Sampled X are not available anymore
      Xavail[indexesSampledX] <- F
      #Store sampled xBar
      #xBar[[i]][,j]<-X[indexesSampledX]
      xBar[[i]][,j,]<-complete_mat(matrix(X[indexesSampledX,],ncol=D),m_max)
    }  
    
    # TOP-DOWN (commented)
    # #For each component, Gibbs sampler
    # Xavail <- rep(T, n) #Available X (all at beginning) for xBar sampling

    for(j in 1:N) {
      
      # Component j corresponds to layer layer_j
      layer_j <- layers$layer[layers$comp==j]
      
      ########################################
      #          Sampling xBar(START)        #
      
      # TOP-DOWN (commented)
      # #Alternative sampling scheme
      # #Arrange the pseudo-inputs for the current iterations of the MCMC loop
      # # Pseudo-inputs are assigned to each component *WITHOUT* repetition.
      # 
      # listIndexEligX <- list() #Store eligible X for each component
      # #Eligible X: which are X that happen to fall in layer_j layer.
      # #Available X: which are X that are not used in all previous components.
      # listIndexEligX[[j]] <- which(indexCompX[[layer_j]]==j)
      # #Eligible AND available X: 
      # indexesEligAndAvailX <-  which(indexCompX[[layer_j]]==j & Xavail) 
      # #Sample from eligible AND available X:
      # if(length(indexesEligAndAvailX)<m_comp[j]){
      #   #If there is no sufficient number of observations in certain component to conduct the sampling, skip whole component.
      #   message('sagpfuns::sagp There is insufficent observations in ')
      #   message(paste0('  OUTPUT: This happens for comp ',j,' at layer',layer_j,'!'))
      #   m_compActual[j]<-0
      #   next
      # }else{
      #   indexesSampledX <- sample(indexesEligAndAvailX, m_comp[j])
      # }
      # #Sampled X are not available anymore
      # Xavail[indexesSampledX] <- F
      # #Store sampled xBar
      # xBar[[i]][,j]<-complete_vec(X[indexesSampledX],m_max)
    
      #          Sampling xBar(END)          #
      ########################################
      
      # Indexes of the X corresponding to component j
      indexesEligX <- listIndexEligX[[j]]
      #This xBar_noNA fetches the xBar, whose empty entries are NA, into a shorter vector by dropping the NAs.
      #Recall that NAs comes from the complete_vec function that allows us to record the xBar in a matrix form.
      xBar_noNA<-matrix(na.omit(xBar[[i]][,j,]),ncol=D)
      
      ########################################
      #          Sampling fBar(START)        #
      # Gibbs step
      
      #Compute residuals of component j given all the other components.
      if(N==1){
        D_mcmc[,j] = Y
      }else{
        idx=(1:N)[-j]
        D_mcmc[,j] = Y-apply(mu.component.temp[,idx,drop=FALSE],1,sum)
      }
      
      
      #Training data of the j-th component: only on data points (X,Y) within that component
      X_train <- matrix(X[indexesEligX,],ncol=D)
      Y_train <- matrix(Y[indexesEligX],ncol=1)
      D_mcmc_train <- D_mcmc[indexesEligX,j]
      
      #Compute covariance matrices
    
      Ktrainm <- makeKernel(X1=X_train,
                            X2=xBar_noNA,
                            kinfo=list(eta=eta[i,j],rho=rho[i,j,]))
      Kmm <- makeKernel(X1=xBar_noNA,
                        X2=xBar_noNA,
                        kinfo=list(eta=eta[i,j],rho=rho[i,j,]))
      Knewm <- makeKernel(X1=Xnew,
                          X2=xBar_noNA,
                          kinfo=list(eta=eta[i,j],rho=rho[i,j,]))
      
      #Inverse of Kmm
      # Trick in cov/apxSparse.m of matlab code
      # Any sparse approximation to the posterior Gaussian process is equivalent to
      # using the approximate covariance function:
      #   Kt = Q + s*diag(g); g = diag(K-Q); Q = Ku'*inv(Kuu+diag(snu2))*Ku;
      # where Ku and Kuu are covariances w.r.t. to inducing inputs xu and
      # snu2 is a vector with the noise variance of the inducing inputs.
      # As a default, we fix the standard deviation of the inducing inputs
      # snu = sfu/1e3 * ones(nu,1) to be a one per mil of the signal standard
      # deviation sfu^2 = trace(Kuu)/nu of the inducing inputs. Alternatively, the
      # noise of the inducing inputs can be treated as a hyperparameter by means
      # of the field hyp.snu that can contain a scalar value log(snu) or a vector
      # of inducing point specific log noise standard deviations.
      sfu2   <- sum( diag(Kmm)/m_comp[j] )
      sfu    <- sqrt(sfu2)
      snu    <- sfu/1e3
      snu2   <- snu^2
      Kmm    <- Kmm +  snu2*diag(m_comp[j])
      KmmInv <- chol2inv(chol(Kmm,pivot=F))
      
      #Compute matrix Lambda_j
      # Definition: Lambda_j =  diag(diag(1/eta[i,j] * I_train - Ktrainm %*% KmmInv %*% t(Ktrainm)))
      # Computation of Ktrainm %*% KmmInv %*% t(Ktrainm) is waste of computations: we only need the diagonal of this. 
      # Element (iDiag,iDiag) of Ktrainm %*% KmmInv %*% t(Ktrainm) is Ktrainm[iDiag,] %*% KmmInv %*% Ktrainm[iDiag,]
      # Compute only the diagonal elements with a for loop
      nTrain <- length(Y_train) # number of training data points for component j
      diag_Knm_KmmInv_Kmn <- rep(NA,nTrain)
      for(iDiag in 1:nTrain) { 
        diag_Knm_KmmInv_Kmn[iDiag] <- matrix(Ktrainm[iDiag,],ncol=m_comp[j]) %*% KmmInv %*% matrix(Ktrainm[iDiag,],nrow=m_comp[j])
      }
      Lambda_diag[[j]] <- rep(1/eta[i,j],nTrain) - diag_Knm_KmmInv_Kmn

      #Perturbation of diagonal of Lambda_j with noise (precision lambda, variance 1/lambda)
      LambdaPerturbed_diag   = Lambda_diag[[j]] + (1/lambda[i])
      #Diagonal of (Lambda_j + noise)^-1 
      LambdaPerturbedInv_diag = (1/LambdaPerturbed_diag)
      
      #Compute Q
      # Definition: Q = Kmm + t(Ktrainm) %*% LambdaPerturbedInv %*% Ktrainm 
      # Take advantage of the fact that LambdaPerturbedInv is diagonal matrix => 
      # LambdaPerturbedInv %*% Ktrainm = 
      #    [ diag(diaLambdaPerturbedInv) * Ktrainm[,1] | diag(diaLambdaPerturbedInv) * Ktrainm[,m] | ... | diag(diaLambdaPerturbedInv) * Ktrainm[,m] ],
      # where diag(diaLambdaPerturbedInv) %*% Ktrainm[,k] is the elementwise product of the k-t column of Ktrainm times  diag(diaLambdaPerturbedInv).
      # This is computed in R with Ktrainm * LambdaPerturbedInv_diag 
      LambdaPerturbedInv_Ktrainm <- Ktrainm * LambdaPerturbedInv_diag
      Q <- Kmm + t(Ktrainm) %*% LambdaPerturbedInv_Ktrainm
      
      Qinv           = chol2inv(chol(Q,pivot=F))
      #M_j is just auxillary matrix to reduce computation.
      M_j            = Kmm %*% Qinv %*% t(Ktrainm)
      
      #Compute Mean_j       
      # Definition: Mean_j  = M_j%*% LambdaPerturbedInv %*% D_mcmc_train
      # As for the definition of Q, take advantage of diagonal structure of LambdaPerturbedInv.
      # M_j %*% LambdaPerturbedInv = elementwise product of rows of M_j times LambdaPerturbedInv_diag = t(t(M_j) * LambdaPerturbedInv_diag)
      Mean_j         = t(t(M_j) * LambdaPerturbedInv_diag) %*% D_mcmc_train
      
      M_j            = Kmm %*% Qinv %*% t(Kmm)
      #Symmetrizer to avoid numerical issue.
      Var_j          = (1/2) * (M_j+t(M_j))
      
      #Sample fBar
      fBar[[i]][,j]  = complete_vec(mvtnorm::rmvnorm(1,mean=Mean_j,sigma=Var_j),m_max)
      #This fBar_noNA fetches the fBar, whose empty entries are NA, into a shorter vector by dropping the NAs.
      #Recall that NAs comes from the complete_vec function that allows us to record the fBar in a matrix form.
      fBar_noNA<-fBar[[i]][,j][!is.na(fBar[[i]][,j])]
      
      #Define prediction of component j with STEP-TAPERING: mean is 0 outside the domain of the component j.
      # Observed locations: avoid computing Knm - mean !=0 only on Xtrain
      mu.component.temp[,j] <- 0
      mu.component.temp[indexesEligX,j]  = (Ktrainm %*% KmmInv %*% fBar_noNA) 
      # New locations: force mean to be 0 outside domain
      mu.new.component.temp[,j]          = (Knewm %*% KmmInv %*% fBar_noNA) * (indexCompNew[[layer_j]] %in% j)
      
      # 
      #          Sampling fBar(END)          #
      ########################################

      
      
      ########################################
      #          Sampling rho(START)         #
      # MH step
      if(sampleRho) {
        #This commented part is not updated to multidimensional case
        # alpha   = pinfo$rho[[layer_j]]$alpha
        # beta    = pinfo$rho[[layer_j]]$beta
        # current_rho = rho[i,j]
        # 
        # proposal_rho  <-  runif(1, min = (current_rho-semiWidthRho[j]), max = (current_rho+semiWidthRho[j]))
        # 
        # #If acceptable proposal: compute acceptance prob. Otherwise, acceptance prob = 0.
        # if(proposal_rho>0 && proposal_rho<1) {
        # 
        #   result_logLikProp <- logLikelihood(rho_j=proposal_rho,
        #                                      Y=Y_train,X=X_train,xBar_j=xBar_noNA,fBar_j=fBar_noNA,eta_j=eta[i,j],lambda=lambda[i],
        #                                      mu.component.temp = matrix(mu.component.temp[indexesEligX,],ncol=N),j_index=j)
        #   result_logLikCurr <- logLikelihood(rho_j=current_rho,
        #                                      Y=Y_train,X=X_train,xBar_j=xBar_noNA,fBar_j=fBar_noNA,eta_j=eta[i,j],lambda=lambda[i],
        #                                      mu.component.temp = matrix(mu.component.temp[indexesEligX,],ncol=N),j_index=j)
        # 
        #   logPriorProp <- dbeta(proposal_rho, alpha, beta, log = T)
        #   logPriorCurr <- dbeta(current_rho, alpha, beta,log = T)
        # 
        #   logPriorProp_fBar_j <- mvtnorm::dmvnorm(fBar_noNA,mean=rep(0,m_comp[j]),sigma=result_logLikProp$Kmm, log = T)
        #   logPriorCurr_fBar_j <- mvtnorm::dmvnorm(fBar_noNA,mean=rep(0,m_comp[j]),sigma=result_logLikCurr$Kmm, log = T)
        # 
        #   if(logPriorProp_fBar_j==-Inf && logPriorCurr_fBar_j==-Inf) {
        #     logPriorProp_fBar_j <- 0
        #     logPriorCurr_fBar_j <- 0
        #   }
        # 
        #   acceptanceProb <- min(1,exp(result_logLikProp$logLik + logPriorProp + logPriorProp_fBar_j - result_logLikCurr$logLik - logPriorCurr - logPriorCurr_fBar_j))
        # } else {
        #   acceptanceProb <- 0
        # }
        # 
        # #Number of rhos sampled
        # sampledRhos[j]<-sampledRhos[j]+1
        # 
        # #Accept or reject with prob acceptanceProb
        # if(runif(1) < acceptanceProb){
        #   
        #   #Accept! Count of accepted value + 1 
        #   acceptCountRhos[j]<-acceptCountRhos[j]+1
        #   rho[i,j]<-proposal_rho
        # 
        #   #If I update rho_j, I also need to update the mean 
        #   Knewm<-makeKernel(X1=Xnew,
        #                     X2=xBar_noNA,
        #                     kinfo=list(eta=eta[i,j],rho=rho[i,j]))
        #   
        #   #Define prediction of component j with STEP-TAPERING: mean is 0 outside the domain of the component j.
        #   mu.component.temp[,j] = 0     
        #   mu.component.temp[indexesEligX,j]  = result_logLikProp$mu.component.temp.new[,j]
        #   mu.new.component.temp[,j] = (Knewm %*% result_logLikProp$KmmInv %*% fBar_noNA) * (indexCompNew[[layer_j]] %in% j)
        # 
        # }else{
        #   #Reject! 
        #   rho[i,j]<- current_rho
        # }
        # 
        # #Adapt semiWidthRho[j]:
        # # - every burnIn/20 iterations.
        # # - leave last 5% of burnin iterations without adapting semiwidth of proposal.
        # if(  i %% (burnIn/20) == 0 && i < (burnIn*19/20+1)  )
        # {
        #   rate_rho = acceptCountRhos[j]/sampledRhos[j]
        #   if(rate_rho>.49 || rate_rho<.39) {
        #     semiWidthRho[j]=semiWidthRho[j]*rate_rho/.44
        #   }
        #   acceptCountRhos[j] <- 0
        #   sampledRhos[j] <- 0
        # }
      
      }else{
        
        #If rho is not sampled, set it to fixed value
        rho[i,j,] <- rep(pinfo$rho[[layer_j]], D)
        
      } 
      
      #           Sampling rho(END)          #
      ########################################
      
        
      ########################################
      #          Sampling eta(START)         #
      # MH step
      if(sampleEta) {
        a_eta   = pinfo$eta[[layer_j]]$a
        b_eta   = pinfo$eta[[layer_j]]$b
        current_eta = eta[i,j]

        proposal_eta  <-  runif(1, min = (current_eta-semiWidthEta[j]), max =  (current_eta+semiWidthEta[j]))

        #If acceptable proposal: compute acceptance prob. Otherwise, acceptance prob = 0.
        if(proposal_eta>0) {
          
          result_logLikProp <- logLikelihood(eta_j = proposal_eta,
                                             Y=Y_train,X=X_train,xBar_j=xBar_noNA,fBar_j=fBar_noNA,rho_j=rho[i,j,],lambda=lambda[i],
                                             mu.component.temp = matrix(mu.component.temp[indexesEligX,],ncol=N),j_index=j)
          result_logLikCurr <- logLikelihood(eta_j = current_eta,
                                             Y=Y_train,X=X_train,xBar_j=xBar_noNA,fBar_j=fBar_noNA,rho_j=rho[i,j,],lambda=lambda[i],
                                             mu.component.temp = matrix(mu.component.temp[indexesEligX,],ncol=N),j_index=j)

          logPriorProp <- dgamma(proposal_eta, a_eta, b_eta,log = T)
          logPriorCurr <- dgamma(current_eta , a_eta, b_eta,log = T)

          logPriorProp_fBar_j <- mvtnorm::dmvnorm(fBar_noNA,mean=rep(0,m_comp[j]),sigma=result_logLikProp$Kmm, log = T)
          logPriorCurr_fBar_j <- mvtnorm::dmvnorm(fBar_noNA,mean=rep(0,m_comp[j]),sigma=result_logLikCurr$Kmm, log = T)

          if(logPriorProp_fBar_j==-Inf && logPriorCurr_fBar_j==-Inf) {
            logPriorProp_fBar_j <- 0
            logPriorCurr_fBar_j <- 0
          }

          acceptanceProb <- min(1,exp(result_logLikProp$logLik + logPriorProp + logPriorProp_fBar_j - result_logLikCurr$logLik - logPriorCurr - logPriorCurr_fBar_j))
        } else {
          acceptanceProb <- 0
        }
        print(acceptanceProb)
        #Number of etas sampled
        sampledEtas[j]<-sampledEtas[j]+1
        
        #Accept or reject with prob acceptanceProb
        if(runif(1) < acceptanceProb){

          #Accept! Count of accepted value + 1 
          acceptCountEtas[j]<-acceptCountEtas[j]+1
          eta[i,j]<- proposal_eta

          #If I update eta_j, I need to update the mean
          Knewm<-makeKernel(X1=Xnew,
                            X2=xBar_noNA,
                            kinfo=list(eta=eta[i,j],rho=rho[i,j,]))

          #Define prediction of component j with STEP-TAPERING: mean is 0 outside the domain of the component j.
          mu.component.temp[,j] = 0     
          mu.component.temp[indexesEligX,j]  = result_logLikProp$mu.component.temp.new[,j]
          mu.new.component.temp[,j] = (Knewm %*% result_logLikProp$KmmInv %*% fBar_noNA) * (indexCompNew[[layer_j]] %in% j)
          #print(Knewm)
        } else {
          #Reject!
          eta[i,j]<- current_eta
        }
        
        #Adapt semiWidthEta[j]:
        # - every burnIn/20 iterations.
        # - leave last 5% of burnin iterations without adapting semiwidth of proposal.
        if( i %% (burnIn/20) == 0 && i < (burnIn*19/20+1) ) {
          rate_eta = acceptCountEtas[j]/sampledEtas[j]
          if(rate_eta>.49 || rate_eta<.39) {
            semiWidthEta[j]=semiWidthEta[j]*rate_eta/.44
          }
          acceptCountEtas[j] <- 0
          sampledEtas[j] <- 0
        }
        
      } else {
        #Set eta to fixed value
        eta[i,j] <- 1
      } 
      
      #           Sampling eta(END)          #
      ########################################  
    }
    
    # Save in-sample prediction of this iteration
    mu.temp      <- apply(mu.component.temp    ,1,sum)
    mu.new.temp  <- apply(mu.new.component.temp,1,sum)
    
    ########################################  
    #        Sampling lambda(START)        #
    
    # Draw lambda, reciprocal of error variance, using gamma normal conjugacy.
    # D_mcmc_overall is the residual in this step of MCMC loop.
    D_mcmc_overall     <- Y-mu.temp
    lambda[i]          <- rgamma(1,pinfo$a+n/2,pinfo$b+0.5*t(D_mcmc_overall)%*%D_mcmc_overall)
    
    #         Sampling lambda(END)         #
    ########################################  
    
    #Store predictions only after burnin to save space
    if(i>burnIn) {
      
      #Store predictions at observed locations only if requested (it requires *A LOT* of space if n large)
      if(storePredX==T) {
        #Mean in observed location
        mu[,i-burnIn] <- mu.temp*sdY + meanY
        #Prediction in observed location
        Yhat[,i-burnIn] <- (mu.temp + rnorm(n,sd=sqrt(1/lambda[i])))*sdY + meanY
        #Mean of each component in observed location
        mu.component[[i-burnIn]] <- mu.component.temp*sdY
      } else {
        mu.component[[i-burnIn]] <- NULL
      }
      
      #Predictions of mean in new locations
      mu.new[,i-burnIn] <- mu.new.temp*sdY + meanY
      #Prediction in observed location
      Yhat.new[,i-burnIn] <- (mu.new.temp + rnorm(nNew,sd=sqrt(1/lambda[i])))*sdY + meanY
      #Prediction of mean of each component in new locations
      mu.new.component[[i-burnIn]]<- mu.new.component.temp*sdY
    }
    
  }
  close(pb)
  
  return(list(#data provided as input:
              Y=Y*sdY + meanY, 
              X=X, 
              Xnew=Xnew, 
              #parameters provided as input:
              pinfo=pinfo, minfo=minfo, tinfo=tinfo,
              sampleRho=sampleRho, sampleEta=sampleEta,storePredX=storePredX,
              #MCMC samples:
              mu=mu,                                #mu has Nmcmc elements, only stored the samples after burn-in
              mu.component=mu.component,            #mu.component has Nmcmc elements, only stored the samples after burn-in
              Yhat=Yhat,                            #Yhat has Nmcmc elements, only stored the samples after burn-in
              mu.new=mu.new,                        #mu.new has Nmcmc elements, only stored the samples after burn-in
              mu.new.component=mu.new.component,    #mu.new.component has Nmcmc elements, only stored the samples after burn-in
              Yhat.new=Yhat.new,                    #Yhat.new has Nmcmc elements, only stored the samples after burn-in
              xBar=xBar[(burnIn+1):(burnIn+Nmcmc)],
              fBar=fBar[(burnIn+1):(burnIn+Nmcmc)],          
              draw.lambda=lambda[(burnIn+1):(burnIn+Nmcmc)],
              draw.rho=rho[(burnIn+1):(burnIn+Nmcmc),,],
              draw.eta=eta[(burnIn+1):(burnIn+Nmcmc),],
              #Info of MH samples
              acceptRateRhos=acceptCountRhos/sampledRhos,
              acceptRateEtas=acceptCountEtas/sampledEtas,
              semiWidthRho=semiWidthRho,
              semiWidthEta=semiWidthEta,
              #Data frame for correspondence between layer and components--e.g., layer=data.frame(comp=c(1,2,3),layer=c(1,2,2))
              layers=layers,
              m_comp=m_comp,m_compActual=m_compActual)
         )
}

# logLikelihood
# log-Likelihood of the model, used in MH steps to sample rhos and etas.
logLikelihood <- function(Y,X,xBar_j,fBar_j,rho_j,eta_j,lambda,mu.component.temp,j_index){
  
  Knm <- makeKernel(X1=X,
                    X2=xBar_j,
                    kinfo=list(eta=eta_j,rho=rho_j))
  Kmm <- makeKernel(X1=xBar_j,
                    X2=xBar_j,
                    kinfo=list(eta=eta_j,rho=rho_j))
  
  #Inversion of Kmm: same trick as before.
  sfu2   <- sum( diag(Kmm)/(dim(Kmm)[1]) )
  sfu    <- sqrt(sfu2)
  snu    <- sfu/1e3
  snu2   <- snu^2
  Kmm    <- Kmm +  snu2*diag((dim(Kmm)[1]))
  KmmInv <- chol2inv(chol(Kmm))
  
  #Update the j-th component of the mean
  mu.component.temp[,j_index] = Knm %*% KmmInv %*% fBar_j
  
  #Computation of Lambda_j: same trick as before
  diag_Knm_KmmInv_Kmn <- rep(NA,nrow(X))
  for(iDiag in 1:nrow(X)) { 
    diag_Knm_KmmInv_Kmn[iDiag] <- matrix(Knm[iDiag,],ncol=(dim(Kmm)[1])) %*% KmmInv %*% matrix(Knm[iDiag,],nrow=(dim(Kmm)[1]))
  }
  Lambda_j_diag <- rep(1/eta_j,nrow(X)) - diag_Knm_KmmInv_Kmn
  
  #log-likelihood 
  # We should compute this:
  # logLikelihood <- mvtnorm::dmvnorm(x=as.vector(Y), 
  #                                   mean=rowSums(mu.component.temp), 
  #                                   sigma=diag(Lambda_j_diag + 1/lambda ), log = T)  
  # However, the covariance matrix is diagonal and mvtnorm::dmvnorm does not take advantage of the diagonal structure to compute its inverse.
  # The computational time was prohibitive for n~10,000. Simply replace the evaluation of the log-likelihood with the pdf of 
  # multivariate normal
  
  #Vector of (Y[i] - X[i])^2
  vector_YminMeaSq <- (as.vector(Y) - rowSums(mu.component.temp))^2
  #Vector of sigma^2[i]
  vectorVariances <- (Lambda_j_diag + 1/lambda )
  #log-likelihood
  logLikelihood <- - .5 * sum(log(vectorVariances)) - 0.5 * nrow(X) * log(2 * pi) + #determinant 
    - 0.5 * sum(vector_YminMeaSq/vectorVariances) #part in exponent  
  
  return(list(logLik=logLikelihood,
              mu.component.temp.new=mu.component.temp,
              Kmm=Kmm,
              KmmInv=KmmInv))
  
}


#generateTree
# generate a partition tree in the domain, given:
# - domain: list of length D, each element is a vector c(Xmin,Xmax) in that dimension
# - the number of layers L
# - base: vector of length D, each element tells how many times the tree should partition that dimension
generateTree <- function(domain, L, base=rep(2,length(domain))) {
  
  #Dimension
  D <- length(domain)
  
  tinfo <- list()
  #For each layer:
  for(l in 1:L) {
    #Compute the coordinates of the centers of the components in each dimensions
    # and store them in the list list_loc_l 
    list_loc_l <- list()
    for(d in 1:D) {
      loc_temp <- seq(from = domain[[d]][1], to = domain[[d]][2], length.out = (2*base[[d]]^(l-1) + 1) )
      #Take the even components: they are the locations
      loc_d <-loc_temp[seq(from=2,to=length(loc_temp),by=2)]
      list_loc_l[[d]] <- loc_d
    }
    #The centers of the components are the D-dimensional grid based on the points in list_loc_l
    tinfo[[l]] <- list(loc = as.matrix(expand.grid(list_loc_l)))
    
    tinfo[[l]]$rad <- matrix(NA, nrow=nrow(tinfo[[l]]$loc), ncol=D)
    for(d in 1:D) {
      tinfo[[l]]$rad[,d] <- (domain[[d]][2]-domain[[d]][1])/(2*base[[d]]^(l-1))
    }
    
  }  
  return(tinfo)
}


#ids_X_comp
#Extract the ids of the X that belong to the component with 
#center 'loc' and radius 'rad'
ids_X_comp <- function(X,loc,rad) {
  D=ncol(X)
  ids_output <- which(rowSums(sapply(1:D, 
                                     function(id_col) {
                                       abs(X[,id_col, drop = FALSE] - loc[id_col]) <= rad[id_col]
                                       }
                                     )
                              ) == D)
  return(ids_output)
}


ids_X_comp_old <- function(X,loc,rad) {
  absDiff_X_loc <- matrix(apply(X, FUN=function(x){abs(x-loc)}, 1),ncol=nrow(X))
  boolean_by_dim <- matrix(apply(absDiff_X_loc,FUN=function(x){x<=rad},2),ncol=nrow(X))
  boolean_X_comp <- apply(boolean_by_dim, FUN=all, 2)
}

