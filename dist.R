########################################

dst <- function(Y, infd=F){
	#calculates distances between nodes based on path length
	#returns d=g if nodes are not connected

	g <- dim(Y)[1]
	Dst=Yr <- Y
	Dst <- Y*(Y==1) + g*(Y==0)
	for(r in 2:(g-1)) {
		Yr <- Yr%*%Y
		Dst <- Dst+(r-g)*( Yr>0 & Dst==g )
    }

	if(infd==T){
		for(i in 1:g){ 
			for(j in 1:g){ 
				if( Dst[i,j]==g ){
					Dst[i,j] <- Inf} 
			}
		}
	}
	diag(Dst) <- 0

	Dst  
}


#######################################


proc.crr <- function(Z,Z0,k=2){
	#Procrustes transform; gives rotation,reflection,trranslation 
	#of Z closest to Z0

	for(i in 1:k) { 
		Z[,i] <- Z[,i]-mean(Z[,i]) + mean(Z0[,i]) 
	} #translation

	A <- t(Z)%*%(  Z0%*%t(Z0)  )%*%Z
	eA <- eigen(A,symmetric=T)
	Ahalf <- eA$vec[,1:k]%*%diag(sqrt(eA$val[1:k]))%*%t(eA$vec[,1:k])

	t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z)) 
}



########################################

lpY <- function(Y,Z=0,alpha){
	#log probability of the graph 

	#Y is sociomatrix
	#Z are latent positions
	#alpha is intercept


	lpZ <- lpz.dist(Z)    #the function of Z that is in the linear predictor

	lpg <- alpha+lpZ

	diag(lpg) <- 0
	sum( Y*lpg - log( 1+exp(lpg) ) ) + dim(Y)[1]*log(2) 
}

#######################################

mlpY <- function(avZ,Y){
	#minus log prob of graph, with Z in vector form, 
	#to be used in by optim or nlm

	alpha <- avZ[1]
	Z <- matrix( avZ[-1],nrow=n,ncol=k)
	lpZ <- lpz.dist(Z)    #the function of Z that is in the linear predictor
	lpg <- alpha+lpZ
	diag(lpg) <- 0
	-(  sum( Y*lpg - log( 1+exp(lpg) ) ) + dim(Y)[1]*log(2)  )
}


########################################

lpz.dist <- function(Z){
	##gives the negative distance between nodes
	ZtZ <- Z%*%t(Z)

	mg <- as.matrix(diag(ZtZ))%*%rep(1,length(Z[,1])) #distances
	mg <- mg+t(mg)
	d <- sqrt((mg-2*ZtZ))
	-d             
}

########
Z.up <- function(Y,Z,alpha,zdelta,mu.z=0,sd.z=10){
	##update Z
	Znew <-Z+matrix(rnorm(k*n,0,zdelta),nrow=n,ncol=k)

	lnew <- lpY(Y,Znew,alpha)
	lold <- lpY(Y,Z,alpha)

	hr <- lnew-lold+sum( dnorm(Znew,mu.z,sd.z,log=T) )-sum( dnorm(Z,mu.z,sd.z,log=T) )
	if( runif(1)>exp(hr) ){
		Znew <- Z
        lnew <- lold 
    }
	list(Z=Znew,lik=lnew)
}

########################################

alpha.up <- function(Y,Z,alpha,adelta,a.a=1,a.b=1){
	##update alpha
	alphanew <- abs(alpha+runif(1,-adelta,adelta) )

	lnew <- lpY(Y,Z,alphanew)
	lold <- lpY(Y,Z,alpha)

	hr <- exp( lnew-lold )*( dgamma(alphanew,a.a,scale=a.b)/dgamma(alpha,a.a,scale=a.b) ) 

	if(runif(1)>hr){ 
		alphanew <- alpha
        lnew <- lold 
    }

	list(alpha=alphanew,lik=lnew)  
}

########################################

