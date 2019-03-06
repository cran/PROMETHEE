

PROMETHEE_OW <- function(dataset,PreferenceF,PreferenceT,IndifferenceT,Weights,Min_Max,S_Gauss,dir)

{  
  
  
  ################################ CHECK #########################
  
  n_indic <- dim(as.matrix(dataset))[2]
  n_unit <- dim(as.matrix(dataset))[1]
  
  # Numeric check dataset
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(dataset[,i]))
    {
      stop(paste("dataset: Data set not numeric at column:",i))
    }
  } 

  # NA check dataset
  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (is.na(dataset[i,j]))
      {
        message(paste("Pay attention: NA values at column:",i,", row",j))
      }
    }
  }  
  
  # check for PreferenceF
  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
	if (!(PreferenceF[i,j] %in% c("Gaussian","Linear","V-shape","Level" )))
      {
        message(paste("Pay attention: PreferenceF has to be set equal to Gaussian,Linear,V-shape or Level"))
      }
    }
  }
  
  # Numeric check PreferenceT
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(PreferenceT[,i]))
    {
      stop(paste("PreferenceT: Data set not numeric at column:",i))
    }
  } 
  
  # Numeric check IndifferenceT
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(IndifferenceT[,i]))
    {
      stop(paste("IndifferenceT: Data set not numeric at column:",i))
    }
  } 
  
  # Numeric check Weights
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(Weights[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    }
  } 
  
  # check for Min_Max
  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (!(Min_Max[i,j] %in% c("min","max")))
      {
        message(paste("Min_Max: Pay attention: Min_Max has to be set equal to min or max (lower case)"))
      }
    }
  }
  
  # Numeric check S_Gauss
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(S_Gauss[,i]))
    {
      stop(paste("S_Gauss: Data set not numeric at column:",i))
    }
  } 
  
  
  ###############################################################


      k=ncol(dataset)
      n=nrow(dataset)
      
      #################
      ########PROMETHEE 
      
      outranking<-matrix(0,nrow=n,ncol=k)
      nonoutranking<-matrix(0,nrow=n,ncol=k)
      UnicriterionNetFlows<-matrix(0,nrow=n,ncol=k)
      
      #pairwise comparisons creating outranking and non outranking matrices 
      for(j in 1:k){
        for(i in 1:n){
         for(a in 1:n){
            if(a!=i){
              d<-dataset[i,j]-dataset[a,j]
              d2<-dataset[a,j]-dataset[i,j]
              if (PreferenceF[i,j]=='V-shape'&Min_Max [i,j]=='max'){
              if (d>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>0)){outranking[i,j]<-outranking[i,j]+(((d)/PreferenceT[i,j])*Weights[i,j])}
              }
              else if (PreferenceF[i,j]=='V-shape'&Min_Max [i,j]=='min'){
              if (d2>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>0)){outranking[i,j]<-outranking[i,j]+(((d2)/PreferenceT[i,j])*Weights[i,j])}
              }
              else  if (PreferenceF[i,j]=='Linear'&Min_Max [i,j]=='max'){
              if (d>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>IndifferenceT[i,j])){outranking[i,j]<-outranking[i,j]+(((d-IndifferenceT[i,j])/(PreferenceT[i,j]-IndifferenceT[i,j]))*Weights[i,j])}
              }
              else  if (PreferenceF[i,j]=='Linear'&Min_Max [i,j]=='min'){
              if (d2>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>IndifferenceT[i,j])){outranking[i,j]<-outranking[i,j]+(((d2-IndifferenceT[i,j])/(PreferenceT[i,j]-IndifferenceT[i,j]))*Weights[i,j])}
              }
              else if (PreferenceF[i,j]=='Level'&Min_Max [i,j]=='max'){
              if (d>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>IndifferenceT[i,j])){outranking[i,j]<-outranking[i,j]+(Weights[i,j]/2)}
              }
              else  if (PreferenceF[i,j]=='Level'&Min_Max [i,j]=='min'){
              if (d2>PreferenceT[i,j]){outranking[i,j]<-outranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>IndifferenceT[i,j])){outranking[i,j]<-outranking[i,j]+(Weights[i,j]/2)}
              }
              else if (PreferenceF[i,j]=='Gaussian'&Min_Max [i,j]=='max'){
              if (d>0){outranking[i,j]<-outranking[i,j]+((1-exp(-((d^2)/(2*S_Gauss[i,j]^2))))*Weights[i,j])}
              else if (d<=0){outranking[i,j]<-outranking[i,j]+0}
              }
              else  if (PreferenceF[i,j]=='Gaussian'&Min_Max [i,j]=='min'){
              if (d2>0){outranking[i,j]<-outranking[i,j]+((1-exp(-((d2^2)/(2*S_Gauss[i,j]^2))))*Weights[i,j])}
              else if (d2<=0){outranking[i,j]<-outranking[i,j]+0}
              }
              
              if (PreferenceF[i,j]=='V-shape'&Min_Max [i,j]=='max'){
              if (d2>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>0)){nonoutranking[i,j]<-nonoutranking[i,j]+(((d2)/PreferenceT[i,j])*Weights[i,j])}
              }
              else  if (PreferenceF[i,j]=='V-shape'&Min_Max [i,j]=='min'){
              if (d>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>0)){nonoutranking[i,j]<-nonoutranking[i,j]+(((d)/PreferenceT[i,j])*Weights[i,j])}
              }
              else   if (PreferenceF[i,j]=='Linear'&Min_Max [i,j]=='max'){
              if (d2>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>IndifferenceT[i,j])){nonoutranking[i,j]<-nonoutranking[i,j]+(((d2-IndifferenceT[i,j])/(PreferenceT[i,j]-IndifferenceT[i,j]))*Weights[i,j])}
              }
              else if (PreferenceF[i,j]=='Linear'&Min_Max [i,j]=='min'){
              if (d>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>IndifferenceT[i,j])){nonoutranking[i,j]<-nonoutranking[i,j]+(((d-IndifferenceT[i,j])/(PreferenceT[i,j]-IndifferenceT[i,j]))*Weights[i,j])}
              }
              else   if (PreferenceF[i,j]=='Level'&Min_Max [i,j]=='max'){
              if (d2>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d2<=PreferenceT[i,j])&(d2>IndifferenceT[i,j])){nonoutranking[i,j]<-nonoutranking[i,j]+(Weights[i,j]/2)}
              }
              else  if (PreferenceF[i,j]=='Level'&Min_Max [i,j]=='min'){
              if (d>PreferenceT[i,j]){nonoutranking[i,j]<-nonoutranking[i,j]+Weights[i,j]}
              else if ((d<=PreferenceT[i,j])&(d>IndifferenceT[i,j])){nonoutranking[i,j]<-nonoutranking[i,j]+(Weights[i,j]/2)}
              }  
              else if (PreferenceF[i,j]=='Gaussian'&Min_Max [i,j]=='max'){
              if (d2>0){nonoutranking[i,j]<-nonoutranking[i,j]+((1-exp(-((d2^2)/(2*S_Gauss[i,j]^2))))*Weights[i,j])}
              else if (d2<=0){nonoutranking[i,j]<-nonoutranking[i,j]+0}
              }
              else  if (PreferenceF[i,j]=='Gaussian'&Min_Max [i,j]=='min'){
              if (d>0){nonoutranking[i,j]<-nonoutranking[i,j]+((1-exp(-((d^2)/(2*S_Gauss[i,j]^2))))*Weights[i,j])}
              else if (d<=0){nonoutranking[i,j]<-nonoutranking[i,j]+0}
              }
            }
          }
        }
      }
      
      UnicriterionNetFlows=outranking-nonoutranking
      
      #################### OW ###########################
      
      n <- nrow(UnicriterionNetFlows)
      Y = UnicriterionNetFlows/(n-1) 
      Y<-as.matrix(Y)
      s <- ncol(Y)
      Risultati <- matrix(, nrow = 1, ncol = (s+1))
      
      ################ Optimistic
      
      if (dir=="Optimistic")
      {
        for(i in 1:n)
          {                
            F_ob<-c(Y[i,])
            Aeq_crs=rbind(c(rep(1,s)))
            const.crs=c("=")
            beq_crs=c(1)
            z=lp(direction = "max", F_ob, Aeq_crs, const.crs, beq_crs)
            Risultati<-rbind(Risultati,cbind(t(z$solution),z$objval))
          }
          Risultati=Risultati[-1,] 
          colnames(Risultati) <- c(colnames(dataset),"GlobalScore")
      }
      
      
      ################ Pessimistic
      
      if (dir=="Pessimistic")
      {
        for(i in 1:n)
        {                
          F_ob<-c(Y[i,])
          Aeq_crs=rbind(c(rep(1,s)))
          const.crs=c("=")
          beq_crs=c(1)
          z=lp(direction = "min", F_ob, Aeq_crs, const.crs, beq_crs)
          Risultati<-rbind(Risultati,cbind(t(z$solution),z$objval))
        }
        Risultati=Risultati[-1,] 
        colnames(Risultati) <- c(colnames(dataset),"GlobalScore")
      }

            
      r<-list(Outranking=outranking,
           Nonoutranking=nonoutranking,
           UnicriterionNetFlows=UnicriterionNetFlows,
           Res=Risultati)

      r$call<-match.call()
      class(r)<-"PROMETHEE_OR"
      return(r)
}      



 
