##set random seed
set.seed(1)

##some functions
source("dist_annot.R")           

##library mva provides cdmscale, used to get starting values
library(MVA)         

##read in Padgett and Ansell's (1993) Florentine family marriage data
Y <- as.matrix(read.table("ff.m.dat"))

##take out unconnected actors
zro <- (1:dim(Y)[1])[ Y%*%rep(1,dim(Y)[1])==0 ] 
if(length(zro)>0) {
	Y <- Y[-zro,-zro]
}


n <- dim(Y)[1]    #number of nodes
k <- 2            #dimension of latent space


## starting values of Z and alpha
D <- dst(Y)
Z <- cmdscale(D, 2)
Z[,1] <- Z[,1]-mean(Z[,1])
Z[,2] <- Z[,2]-mean(Z[,2])

alpha <- 2

## now find a good starting point for the MCMC
avZ <- c(alpha,c(Z))
#simulated annealing 
avZ <- optim(avZ,mlpY,Y=Y,method="SANN")$par
#BFGS 
avZ <- optim(avZ,mlpY,Y=Y,method="BFGS")$par

avZ <- avZ*2/(avZ[1])  #MLE too extreme, keep shape but rescale
alpha <- avZ[1]
Z <- Z.mle <- matrix(avZ[-1],nrow=n,ncol=k)


##see if this Z represents Y
(-1)^(1-Y) * (alpha+lpz.dist(Z) ) >=0

## MCMC params

adelta <- .5          #params for
zdelta <- .2          #proposal distribution

nscan <- 10^6         #number of scans   
odens <- 10^3         #save output every odens step

aca=acz <- 0          #keep track of acceptance rates
Alpha <- alpha        #keep track of alpha
Lik <- lpY(Y,Z,alpha) #keep track of alpha and likelihood 

Z.post <- list()    #keep track of positions
for(i in 1:k){ 
	Z.post[[i]] <- t(Z[,i]) 
}


## MCMC

for(ns in 1:nscan){

	tmp <- Z.up(Y,Z,alpha,zdelta)                     #update z's
	if(tmp$Z!=Z){ 
    	acz <- acz+1/odens
        Z <- proc.crr(tmp$Z,Z.mle)  
    }

	tmp <- alpha.up(Y,Z,alpha,adelta,a.a=2,a.b=1)  #update alpha
	if(tmp$alpha!=alpha) {
		aca <- aca+1/odens
        alpha <- tmp$alpha 
        lik <- tmp$lik
    }

	if( ns%%odens==0 ){                            #output
    	Alpha <- c(Alpha,alpha)
        cat(ns,aca,acz,alpha,lik,"\n")
        Lik <- c(Lik,lik)
        acz=aca <- 0 
        for(i in 1:k){ 
        	Z.post[[i]] <- rbind(Z.post[[i]],t(Z[,i])) 
        } 
    }
}

Post <- list(Z=Z.post,Alpha=Alpha,Lik=Lik)
dput(Post,"ff.m.post") #output to a file

##plot results, for k=2
Zp <- Post$Z
Lik <- Post$Lik
Alpha <- Post$Alpha


##posterior mean 
Z.pm <- Z.mle
for(i in 1:k){  
	Z.pm[,i] <- rep(1/dim(Zp[[i]])[1],dim(Zp[[i]])[1])%*%Zp[[i]] 
}
Z.mle <- Z.mle*sum(diag((Z.pm%*%t(Z.mle))))/sum(diag(( t(Z.mle)%*%Z.mle)))
Z.mle <- proc.crr(Z.mle,Z.pm)

##plot likelihood and alpha
par(mfrow=c(1,2))
plot(Lik,type="l")
plot(Alpha,type="l")

##plot positions and confidence sets
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
par(mgp=c(2,1,0))

##set up color scheme based on a scaled Z.mle


r <- atan2(Z.mle[,2],Z.mle[,1])
r <- r+abs(min(r))
r <- r/max(r)
g <- 1-r
b <- (Z.mle[,2]^2+Z.mle[,1]^2)
b <- b/max(b)

plot(Z.mle[,1],Z.mle[,2],xaxt="n",yaxt="n",xlab="",ylab="",type="n", 
                          xlim=range(Zp[[1]]),ylim=range(Zp[[2]]))

for(i in 1:nrow(Y)){
	for(j in 1:nrow(Y)){ 
		lines(  Z.mle[c(i,j),1] , Z.mle[c(i,j),2] ,lty=Y[i,j]) 
	}
}

text(Z.mle[,1]*1.1,Z.mle[,2]*1.1,labels(Y)[[1]])

for(i in 1:n) {
	points( Zp[[1]][,i],Zp[[2]][,i],
			pch=46,cex=5, col=rgb(r[i],g[i],b[i]) ) 
}



