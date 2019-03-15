pocrm.sim<-function(r,alpha,prior.o,x0,stop,n,theta,nsim,tox.range){
  
sim <- sim1 <- apred <- lik <- pord <- ord <- ahat <- rpred <- next.lev <- n1 <- N <- NULL
  
d<-ncol(alpha)
s<-nrow(alpha)

if(nsim==1){
	
crm<-function(obs,alpha,prior.o,theta){

sim<<-table(obs$level,obs$tox)
ifelse(dim(sim)[2]==1,sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],1-sim[,1]),sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],sim[,2]))
names(sim1)<<-c('level','nontox','tox')

apred<<-rep(0,s)
	lik<<-rep(0,s)
	for(k in 1:s){			
		ll<-function(a){
			la<-0 #value of the log-likelihood
			for(i in sim1$level){
			index<-match(i,sim1$level)
			la<-la+sim1$tox[index]*a*log(alpha[k,][i])+sim1$nontox[index]*log((1-alpha[k,][i]**a))
			}
			la
			}
		apred[k]<<-optimize(f=ll,interval=c(0,100),maximum=T)$maximum
		lik[k]<<-ll(apred[k])
		}
pord<<-(exp(lik)*prior.o)/sum(exp(lik)*prior.o)
ord<<-which.is.max(pord)
ahat<<-apred[ord]
rpred<<-alpha[ord,]**ahat
next.lev<<-which.is.max(-(abs(rpred-theta)))
next.lev
    }
###'crm' ENDS HERE


###LOAD FUNCTION 'twostgcrm'

twostgcrm<-function(r,x0,stop,n,theta){

	n1<<-n+1
	obs<-data.frame(cbind(1:n1,rep(0,n1),rep(0,n1),rep(0,n1),rep(0,n1)))
	names(obs)<-c('patient','level','tox','a','order')

#######Beginning at dose level 1
###1st Stage:  Up and down scheme with cohort size 1

i<-1
#x0<-lapply(zones,ff)
##'initial.scheme' is a vector indicating the Stage I escalation scheme
initial.scheme<-x0
initial.scheme<-c(initial.scheme,rep(ncol(alpha),n1-length(initial.scheme)))

while(i < n1){
obs$order[i]<-99
obs$level[1]<-initial.scheme[1]
p<-runif(1)
#number of tox in 1st patient
index<-p<=r[obs$level[i]] ##determines any toxicities
obs$tox[i]<-obs$tox[i]+as.numeric(index)
if(any(obs$tox[1:i]==1) & any(obs$tox[1:i]==0)){
	q<-2
	break
}
if(all(obs$tox[1:i]==1)){
	i<-i+1
	obs$level[i]<-initial.scheme[1]
}
if(all(obs$tox[1:i]==0)){
	i<-i+1
	obs$level[i]<-initial.scheme[i]
}

if(length(obs$level[obs$level==d])==stop+1){
	MTD<-d
	break
}
}
	

##2nd stage
N<<-table(obs$level>0)[2]+1
if(any(obs$tox>0)){
		  level<-crm(obs[1:(N-1),],alpha,prior.o,theta)
		  obs$a[N-1]<-ahat
	        obs$order[N-1]<-ord

for(j in N:n1){
##assigment for remaining patients
	obs$level[j]<-level
	if(obs$level[n1]>0){
			MTD<-obs$level[n1]
		break
		}

	if(length(obs$level[obs$level==level])==stop+1){
		MTD<-level
		break
		}

	index<-runif(1)<=r[obs$level[j]]
	if(index){obs$tox[j]<-1}
	level<-crm(obs[1:j,],alpha,prior.o,theta)
	obs$a[j]<-ahat
	obs$order[j]<-ord
	
##crm dose recommendation for Nth patient
	}
	} else
	MTD<-d
	out<-list(trial=obs[obs$level>0,],MTD.selection=MTD)
	}
###'twostgcrm' ENDS HERE
}

if(nsim>1){

###Load the function 'lpocrm' 
lpocrm<-function(r,alpha,prior.o,x0,stop,n,theta){
 
# if a single ordering is inputed as a vector, convert it to a matrix
	if(is.vector(alpha)) alpha=t(as.matrix(alpha));
	
	nord.tox = nrow(alpha);
	mprior.tox = prior.o;  # prior for each toxicity ordering

bcrml<-function(a,p1,y,n){
	lik=0
	for(j in 1:length(p1)){
		lik=lik+y[j]*a*log(p1[j])+(n[j]-y[j])*log((1-p1[j]**a));
		}
	return(lik);
    }

### run a trial 	
    ncomb = ncol(alpha);   #number of combos
    y=npts=ptox.hat=comb.select=numeric(ncomb);  
    comb.curr = x0[1];  # current dose level	 
    stoprule=0; #indicate if trial stops early
    i=1

stage1<-c(x0,rep(ncol(alpha),n-length(x0)))

##Stage 1
while(i <= n){
		y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
		npts[comb.curr] = npts[comb.curr] + 1;
		
		if(sum(y)==sum(npts)){
			comb.curr<-ifelse(comb.curr==1,comb.curr,comb.curr-1)
		} else if(sum(y)==0){
			comb.curr<-ifelse(comb.curr==ncomb,comb.curr,stage1[i+1])
		} else {
			break
		}
		if(any(npts>stop)){
			stoprule<-0
			break
		}
i=i+1
}

#Stage 2
while(sum(npts) <= n)
    {
		if(sum(y)==0){
			stop=0
			break
		} else{
		like.tox= est.tox=rep(0, nord.tox);
		for(k in 1:nord.tox)
		{
			est.tox[k]<-optimize(f=bcrml,interval=c(0,100),p1=alpha[k,],y=y,n=npts,maximum=T)$maximum
			like.tox[k]<-optimize(f=bcrml,interval=c(0,100),p1=alpha[k,],y=y,n=npts,maximum=T)$objective
		}		

		postprob.tox = (exp(like.tox)*mprior.tox)/sum(exp(like.tox)*mprior.tox);
		# toxicity model selection, identify the model with the highest posterior prob
		if(nord.tox>1){ 
			mtox.sel = which.is.max(postprob.tox); 
		} else{
			mtox.sel = 1;
		}

		ptox.hat=alpha[mtox.sel,]**est.tox[mtox.sel]
						
		loss=abs(ptox.hat-theta)
		comb.curr=which.is.max(-loss)
		if(npts[comb.curr]==stop){
			stoprule<-0
			break
		}

		if(sum(npts)==n){
			stoprule=0
			break
		} else{
		# generate data for a new cohort of patients
		y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
		npts[comb.curr] = npts[comb.curr] + 1;
			}
		}
	}
	if(stoprule==0){
		comb.select[comb.curr]=comb.select[comb.curr]+1;
		}
	return(list(MTD.selection=comb.select,tox.data=y,patient.allocation=npts))
}
##########'lpocrm' end here
}

###Load the function 'lpocrm.sim' 
lpocrm.sim<-function(nsim){
	ncomb=length(r)
	
	comb.select<-y<-npts<-matrix(nrow=nsim,ncol=ncomb)
	trialsize<-rep(0,nsim)
	nstop=0
	
	for(i in 1:nsim){
		result<-lpocrm(r,alpha,prior.o,x0,stop,n,theta)
		comb.select[i,]=result$MTD.selection
		y[i,]=result$tox.data
		npts[i,]=result$patient.allocation
		trialsize[i]=sum(result$patient.allocation)
	}
return(list(true.prob=r,MTD.selection=round(colMeans(comb.select),2),patient.allocation=round(colMeans(npts)/mean(trialsize),2),percent.DLT=sum(colMeans(y))/mean(trialsize),mean.n=mean(trialsize),acceptable=sum(colMeans(comb.select)[which(round(abs(r-theta),2)<=tox.range)])))
	}
##########'lpocrm.sim' end here
if(nsim==1){
	twostgcrm(r,x0,stop,n,theta)
} else{
	lpocrm.sim(nsim)
}
}


