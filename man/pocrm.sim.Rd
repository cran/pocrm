\name{pocrm.sim}
\alias{pocrm.sim}

\title{
Simulating the PO-CRM
}
\description{
pocrm.sim is used to simulate Phase I trials of combined drugs according to the partial order continual reassessment method (PO-CRM).
}
\usage{
pocrm.sim(r, alpha, prior.o, x0, stop, n, theta, nsim, tox.range)
}

\arguments{
  \item{r}{
A vector of true toxicity probabilities.
}
  \item{alpha}{
A matrix of skeleton values corresponding to the possible orderings of the toxicity probabilities generated by getwm.
}
  \item{prior.o}{
A vector of prior probabilities on the orderings.
}
  \item{x0}{
A vector specifying an initial escalation scheme. It should be a sequence of combinations to which the trial would proceed prior to the observance of the first DLT. The length of x0 cannot exceed n.
}
  \item{stop}{
The total number of patients treated on any combination required to stop the trial. See details.
}
  \item{n}{
The maximum sample size. If stop>n, then the trial will exhaust a pre-determined, fixed sample size of n patients.
}
  \item{theta}{
The target DLT rate.
}
  \item{nsim}{
The number of simulations.
}
  \item{tox.range}{
A single numeric value used to define a range of "acceptable" DLT rates. The simulation results will report the percentage of simulated trials that recommended a combination within +/- tox.range of the target rate.
}
}
\details{
The method bases toxicity probability estimates on the power model (5) of Wages, Conaway and O'Quigley (2011a).

The initial escalation scheme, x0, is a vector specifying the order in which combinations should be tried and the size of the cohort at each combination in Stage 1. For instance, proceeding through the zones in Table 3 of Wages, Conaway and O'Quigley (2011) in cohorts of 2 would x0=c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2),rep(8,2)). This scheme will only be relied upon until the first DLT is observed, when the design switches to Stage 2. In Stage 1, if a DLT is observed on the first patient, the method will continue to enroll patients on the lowest combination until a non-DLT is observed.

The design will stop at the end of Stage 1 if escalation proceeds to the highest dose combination and stop patients are treated with no DLTs. In this case, the study is stopped and the highest dose combination is declared the MTD combination. In the Stage 2, if the recommendation is to assign the next patient to a combination that already has stop patients treated on the combination, the study is stopped and the recommended combination is declared the MTD combination.
}
\value{
  \item{trial}{
Data frame containing output of a single simulated trial. Output includes the combination on which each patient was treated, toxicity information for each patient, the estimated value of the model parameter and the estimated ordering of toxicity probabilities. Only output if nsim=1.
}
  \item{true.prob}{
True toxicity probabilities.
}
  \item{MTD.selection}{
The distribution of MTD estimates. If nsim=1, this is a single numeric value of the recommended MTD combination of a single simulated trial.
}
  \item{patient.allocation}{
The distribution of patient allocation at each combination.
}
  \item{percent.DLT}{
Average percentage of DLTs across all simulated trials.
}
  \item{mean.n}{
Average number of patients treated in the trial.
}
  \item{acceptable}{
Percentage selection of combinations within +/- tox.range of the target rate.
}
}
\references{
Wages, Conaway and O'Quigley (2011). Dose-finding design for multi-drug combinations. Clinical Trials 8(4): 380-389.

Wages, Conaway and O'Quigley (2011a). Continual reassessment method for partial ordering. Biometrics 67(4): 15551563.
}

\examples{

#All specifications refer to example in Wages, Conaway and O'Quigley (2011a).

#True toxicity rates of Scenario 1 in Table 4.
r<-c(0.06,0.08,0.10,0.15,0.10,0.12,0.30,0.45,0.15,0.30,0.50,0.60,0.50,0.55,0.60,0.70)

#Specifiy the possible orderings from Section 3.3.
orders<-matrix(nrow=3,ncol=16)
orders[1,]<-c(1,2,5,3,6,9,4,7,10,13,8,11,14,12,15,16)
orders[2,]<-c(1,5,2,3,6,9,13,10,7,4,8,11,14,15,12,16)
orders[3,]<-c(1,5,2,9,6,3,13,10,7,4,14,11,8,15,12,16)

#Specify the skeleton values.
skeleton<-getprior(0.05,0.30,8,16)

#Initial guesses of toxicity probabilities for each ordering.
alpha<-getwm(orders,skeleton)

#We consider all orders to be equally likely prior to the study.
prior.o<-rep(1/3,3)

#Initial escalation in Stage 1 proceeds according to the zones in Figure 1.
#Single patient cohorts are used.
x0<-c(rep(1,1),rep(2,1),rep(5,1),rep(3,1),rep(6,1),rep(9,1),rep(4,1),rep(7,1),rep(10,1),rep(13,1),
rep(8,1),rep(11,1),rep(14,1),rep(12,1),rep(15,1),rep(16,1))

#Number of patients used to define stopping rule
stop<-61

#Maximum sample size.
n<-60

#The target toxicity rate
theta<-0.30

#Number of simulations
nsim<-2

#Definition of acceptable DLT rates
tox.range<-0.05

fit<-pocrm.sim(r,alpha,prior.o,x0,stop,n,theta,nsim,tox.range)
fit
}
