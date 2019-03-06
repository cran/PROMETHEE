\name{PROMETHEE}
\alias{PROMETHEE}

\title{Preference Ranking Organization METHod for Enrichment of Evaluations}

\description{Functions which can be used to support the Multicriteria Decision Analysis (MCDA) process 
    involving multiple criteria, by PROMETHEE (Preference Ranking Organization METHod for Enrichment of Evaluation)}

\usage{PROMETHEE(dataset, PreferenceF,PreferenceT,IndifferenceT,Weights,Min_Max,S_Gauss)}

\arguments{
  \item{dataset}{A matrix with data (alternatives by row and criteria by column)}
  \item{PreferenceF}{A matrix with preference functions (alternatives by row and criteria by column)}
  \item{PreferenceT}{A matrix with preference thresholds (alternatives by row and criteria by column)}
  \item{IndifferenceT}{A matrix with indifference thresholds (alternatives by row and criteria by column)}
  \item{Weights}{A matrix with weights (alternatives by row and criteria by column)}
  \item{Min_Max}{A matrix that specifies whether the criteratum should be maximized or minimized (alternatives by row and criteria by column)}
  \item{S_Gauss}{A matrix with S Gaussians (alternatives by row and criteria by column)}  
}


\value{
  \item{Outranking}{Outranking matrix (alternatives by row and criteria by column)}
  \item{Nonoutranking}{Non-outranking matrix (alternatives by row and criteria by column)}
  \item{UnicriterionNetFlows}{Unicriterion net flows matrix (alternatives by row and criteria by column)}
  \item{PROMETHEE1}{PROMETHEE I matrix (first column is Positive Flows and second column is Negative Flows, alternatives by row)}
  \item{PROMETHEE2}{PROMETHEE II vector (Global net Flows, alternatives by row)}
  }

\details{Thanks are extended to Salvatore Greco, Alessio Ishizaka, and Gianpiero Torrisi for helpful comments}

\author{Giuliano Resce, Menelaos Tasiou, Francesco Vidoli}


\examples{


# Dataset
dati<-as.data.frame(cbind(c(10,25,4,2),
                          c(250000,270000,320000,350000),
                          c(3,4,2,2),
                          c(20,10,15,7)))

colnames(dati) = c("Distance.to.work","Price","Bedrooms","Age")


# PreferenceF
PreF<-as.data.frame(rbind(c("Gaussian","Linear","V-shape","Level"),
                     c("Gaussian","Linear","V-shape","Level"),
                     c("Gaussian","Linear","V-shape","Level"),
                     c("Gaussian","Linear","V-shape","Level")))
colnames(PreF) = c("Distance.to.work","Price","Bedrooms","Age")


# PreferenceT
PreT<-as.data.frame(cbind(c(2,2,2,2),
                          c(50000,50000,50000,50000),
                          c(2,2,2,2),
                          c(5,5,5,5)))
colnames(PreT) = c("Distance.to.work","Price","Bedrooms","Age")


# IndifferenceT
IndT<-as.data.frame(cbind(c(1,1,1,1),
                          c(10000,10000,10000,10000),
                          c(0,0,0,0),
                          c(2,2,2,2)))
colnames(IndT) = c("Distance.to.work","Price","Bedrooms","Age")


#Weights
Weig<-as.data.frame(cbind(c(0.25,0.25,0.25,0.25),
                          c(0.25,0.25,0.25,0.25),
                          c(0.25,0.25,0.25,0.25),
                          c(0.25,0.25,0.25,0.25)))
colnames(Weig) = c("Distance.to.work","Price","Bedrooms","Age")


# Min_Max
MiMa<-as.data.frame(cbind(c("min","min","min","min"),
                          c("min","min","min","min"),
                          c("max","max","max","max"),
                          c("min","min","min","min")))
colnames(MiMa) = c("Distance.to.work","Price","Bedrooms","Age")


#S_Gauss
gauss<-as.data.frame(cbind(c(2,2,2,2),
                          c(0,0,0,0),
                          c(0,0,0,0),
                          c(0,0,0,0)))
colnames(gauss) = c("Distance.to.work","Price","Bedrooms","Age")


###############################

PF = PROMETHEE(dati, PreF,PreT,IndT,Weig,MiMa,gauss)

###############################

PCA_UNIC <- prcomp(PF$UnicriterionNetFlows,center = TRUE,scale. = TRUE)
GAIA=predict(PCA_UNIC)[,1:2]
rownames(GAIA)=rownames(dati)
print(GAIA)

}
