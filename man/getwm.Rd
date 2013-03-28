\name{getwm}
\alias{getwm}

\title{
Obtain the working models for use in the partial order continual reassessment method (PO-CRM)
}
\description{
getwm returns a matrix of skeleton values corresponding to the possible orderings of the toxicity probabilities in a Phase I trial of combined drugs for use in the partial order continual reassessment method (PO-CRM).
}
\usage{
getwm(orders, skeleton)
}

\arguments{
  \item{orders}{
A matrix specifying the possible orderings of the combinations with regards to their toxicity probabilities.
}
  \item{skeleton}{
A vector of values for the initial guesses of toxicity probabilities. The location of these values will be adjusted to correspond to the orderings specified by orders.
}
}
\details{
Based on the matrix of possible orderings, the skeleton values are adjusted to correspond to the possible orderings of the toxicity probabilities. See Table 4 of Wages, Conaway and O'Quigley (2011).
}
\value{
    \item{alpha}{
A matrix of skeleton values corresponding to the possible orderings of the toxicity probabilities specified by orders.
}
}
\references{
Wages, Conaway and O'Quigley (2011). Dose-finding design for multi-drug combinations. Clinical Trials 8(4): 380-389.
}

\examples{

#All specifications refer to example in Wages, Conaway and O'Quigley (2011).

#Specify the possible orderings from Table 2
orders<-matrix(nrow=8,ncol=8)
orders[1,]<-c(1,2,3,4,5,6,7,8)
orders[2,]<-c(1,3,2,4,5,6,7,8)
orders[3,]<-c(1,2,3,5,4,6,7,8)
orders[4,]<-c(1,2,3,4,5,7,6,8)
orders[5,]<-c(1,3,2,5,4,6,7,8)
orders[6,]<-c(1,3,2,4,5,7,6,8)
orders[7,]<-c(1,2,3,5,4,7,6,8)
orders[8,]<-c(1,3,2,5,4,7,6,8)

#Specify the skeleton values provided in Table 4.
skeleton<-c(0.01,0.03,0.10,0.20,0.33,0.47,0.60,0.70)

getwm(orders,skeleton) 
}
