\name{PROMETHEE_OW}
\alias{PROMETHEE_OW}

\title{Preference Ranking Organization METHod for the Enrichment of Evaluations with Optimal Weights}

\description{Functions which can be used to support the Multicriteria Decision Analysis (MCDA) process 
    involving multiple criteria, by a specific tool allowing joining the consolidated procedure usually employed for environmental evaluation (PROMETHEE), with a more flexible weighing process inspired by DEA}

\usage{PROMETHEE_OW(dataset,PreferenceF,PreferenceT,IndifferenceT,Weights,Min_Max,S_Gauss,dir)}

\arguments{
  \item{dataset}{A matrix with data (alternatives by row and criteria by column)}
  \item{PreferenceF}{A matrix with preference functions (alternatives by row and criteria by column)}
  \item{PreferenceT}{A matrix with preference thresholds (alternatives by row and criteria by column)}
  \item{IndifferenceT}{A matrix with indifference thresholds (alternatives by row and criteria by column)}
  \item{Weights}{A matrix with weights (alternatives by row and criteria by column)}
  \item{Min_Max}{A matrix that specifies whether the criteratum should be maximized or minimized (alternatives by row and criteria by column)}
  \item{S_Gauss}{A matrix with S Gaussians (alternatives by row and criteria by column)}  
  \item{dir}{An element specifying whether the optimization should be "Optimistic" or "Pessimistic"}  
}


\value{
  \item{Outranking}{Outranking matrix (alternatives by row and criteria by column)}
  \item{Nonoutranking}{Non-outranking matrix (alternatives by row and criteria by column)}
  \item{UnicriterionNetFlows}{Unicriterion net flows matrix (alternatives by row and criteria by column)}
  \item{Res}{DEA weights and global score}
  }

\details{Caravaggio, N., Caravella, S., Ishizaka, A., & Resce, G. (2019). Beyond CO2: a multi-criteria analysis of air pollution in Europe. Journal of Cleaner Production.}

\author{Giuliano Resce, Menelaos Tasiou, Francesco Vidoli}


\examples{

data(Austria)

PreferenceF= matrix("Level",nrow(Austria),ncol(Austria) )
PreferenceT=matrix(0,nrow(Austria),ncol(Austria) )
IndifferenceT=matrix(0,nrow(Austria),ncol(Austria) )
Weights=matrix(1,nrow(Austria),ncol(Austria) )
Min_Max=matrix("min",nrow(Austria),ncol(Austria) )
S_Gauss=matrix(0,nrow(Austria),ncol(Austria) )

PF=PROMETHEE_OW(Austria, PreferenceF,PreferenceT,IndifferenceT,Weights,Min_Max,S_Gauss,
                "Optimistic")
PF$Res
PF=PROMETHEE_OW(Austria, PreferenceF,PreferenceT,IndifferenceT,Weights,Min_Max,S_Gauss,
                "Pessimistic")
PF$Res

}
