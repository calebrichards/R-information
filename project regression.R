Data_Cortex_Nuclear <- read.csv("~/Grad School/Genomics/Data_Cortex_Nuclear.csv",header= TRUE,stringsAsFactors = FALSE)

mice.dat <- data.frame(cbind(Data_Cortex_Nuclear$Genotype=="Control",
                                Data_Cortex_Nuclear$Treatment=="Memantine",
                                Data_Cortex_Nuclear$Behavior=="C/S",
                             Data_Cortex_Nuclear$Treatment=="Memantine"&Data_Cortex_Nuclear$Behavior=="C/S",
                                Data_Cortex_Nuclear[,1:78]))
mice.names <- rep("",dim(mice.dat)[2])
mice.names[1] <- "Group"
mice.names[2] <- "Treatment"
mice.names[3] <- "Behavior"
mice.names[4] <- "Learning"
mice.names[5] <- "mouse"



#Check for duplicate rows
copies <- cbind(rep(0,78*78),rep(0,78*78))
idxs <- 1
for(i in 2:78){
  for(j in (i+1):78){
    if(sum(Data_Cortex_Nuclear[,i]==Data_Cortex_Nuclear[,j], na.rm = TRUE)>100){
      copies[idxs,1] <- i
      copies[idxs,2] <- j
      idxs <- idxs+1
    }
  }
}
copies <- copies[copies[,1]!=0,]  


for(i in 6:dim(mice.dat)[2]){
  mice.names[i] <- paste("protein",i,sep = "")
}

mice.dat <- setNames(mice.dat,mice.names)
mice.dat <- mice.dat[mice.dat$Group==FALSE,]

A <- matrix(unlist(strsplit(mice.dat$mouse,"_")),ncol =2,byrow=TRUE)
A <- unique(A[,1])

mice.stats <- matrix(nrow = 34, ncol = 79)
for(i in 1:length(A)){
  
  if(mice.dat$Treatment[i*15]==TRUE){
    mice.stats[i,1] <-  TRUE
  } else{
    mice.stats[i,1] <-  FALSE
  }
  
  if(mice.dat$Treatment[i*15]==TRUE&mice.dat$Behavior[i*15]==TRUE){
    mice.stats[i,2] <- TRUE
  } else {
    mice.stats[i,2] <- FALSE
  }
  
  mice.stats[i,3:79] <- apply(mice.dat[substr(mice.dat$mouse,1,nchar(A[i]))==A[i],6:82],2,mean,na.rm=TRUE)
}
mice.sd <- apply(mice.dat[substr(mice.dat$mouse,1,nchar(A[i]))==A[i],6:82],2,sd,na.rm=TRUE)
#mice.max <- apply(mice.stats,MARGIN = 2,max,na.rm=TRUE)
#mice.min <- apply(mice.stats,MARGIN = 2,min,na.rm=TRUE)

#for(i in 3:dim(mice.stats)[2]){
#  mice.stats[,i] <- (mice.stats[,i]-mice.min[i])/(mice.max[i]-mice.min[i])
#}

pvals <- cbind(rep(0,77),rep(0,77),rep(0,77))
for(i in 3:dim(mice.stats)[2]){
  mice.lm <- lm(mice.stats[,i] ~ mice.stats[,1] + mice.stats[,2])
  pvals[i-2,1:3] <- summary(mice.lm)$coefficients[,4]
}

sig.2 <- (1:77)[pvals[,2]<.01]
sig.3 <- (1:77)[pvals[,3]<.01]

id <- 1
red.id <- rep(0,length(sig.3))
red.id <- sig.3[is.element(sig.3,sig.2)==FALSE]
red.id <- red.id[red.id!=0]

pvals.reduced <- cbind(rep(0,length(red.id)),rep(0,length(red.id)))
for(i in red.id){
  mice.lm <- lm(mice.stats[,i+2] ~ mice.stats[,2])
  pvals.reduced[id,1:2] <- summary(mice.lm)$coefficients[,4]
  id <- id+1
}


sig.ids <- c(red.id[pvals.reduced[,2]<.01],sig.2[is.element(sig.2,sig.3)])
sig.proteins <- names(Data_Cortex_Nuclear)[sig.ids+1]
sig.pvals <- pvals[sig.ids,]

pvals.final <- rep(0,length(sig.ids))
beta.1 <- cbind(rep(0,length(sig.ids)),rep(0,length(sig.ids)))
beta.ci <- cbind(rep(0,length(sig.ids)),rep(0,length(sig.ids)))
sigred.id <- red.id[pvals.reduced[,2]<.01]
for(i in 1:length(sigred.id)){
  mice.lm <- lm(mice.stats[,sig.ids[i]+2] ~ mice.stats[,2])
  beta.1[i,2] <- summary(mice.lm)$coefficients[2,1]
  beta.ci[i,] <- confint(mice.lm,"mice.stats[, 2]",level = 0.95)
  pvals.final[i] <- summary(mice.lm)$coefficients[2,4]
}

remainder.ids <- sig.2[is.element(sig.2,sig.3)]
betafull.ci <- cbind(rep(0,length(remainder.ids)),rep(0,length(remainder.ids)))
for(i in 1:length(remainder.ids)){
  mice.lm <- lm(mice.stats[,remainder.ids[i]+2] ~ mice.stats[,1] + mice.stats[,2])
  beta.1[i+length(sigred.id),] <- summary(mice.lm)$coefficients[2:3,1]
  beta.ci[i+length(sigred.id),] <- confint(mice.lm,"mice.stats[, 2]",level = 0.95)
  betafull.ci[i,] <- confint(mice.lm,"mice.stats[, 1]",level = 0.95)
  pvals.final[i+length(sigred.id)] <- 1-pf(summary(mice.lm)$fstat[1],summary(mice.lm)$df[1]-1,summary(mice.lm)$df[2])
}

a.bon <- .01/76
pvals.corrected <- pvals.final[pvals.final<a.bon]
sig.corrected <- sig.ids[pvals.final<a.bon]
sigprot.corrected <- names(Data_Cortex_Nuclear)[sig.corrected+1]

#write.csv(cbind(sigprot.corrected,pvals.corrected),file = "Corrected Pvalues.csv")
#write.csv(rbind(beta.ci,betafull.ci),file = "Confints.csv")
#write.csv(beta.1,file = "Beta values.csv")
#write.csv(cbind(sig.proteins,pvals.final), file = "Significant Proteins.csv")

sigpvals.reduced <- rep(0,length(sig.ids))
for(i in 1:length(sig.ids)){
  mice.lm <- lm(mice.stats[,sig.ids[i]+2] ~ mice.stats[,2])
  sigpvals.reduced[i] <- summary(mice.lm)$coefficients[2,4]
}

plot(mice.stats[,2],mice.stats[,sig.ids[1]+2],xlim = c(-1,2),xaxt="n",
     xlab = "Learning Outcome",ylab = "Protein Level", main = "Reduced Model Regression for Protein DYRK1A_N")
axis(1, at = c(0,1))
mice.lm <- lm(mice.stats[,sig.ids[1]+2] ~ mice.stats[,2])
abline(mice.lm)

plot(mice.stats[,2],mice.stats[,sig.ids[11]+2],xlim = c(-1,2),xaxt="n",
     xlab = "Learning Outcome",ylab = "Protein Level", main = "Full Model Regression for Protein SOD1_N")
axis(1, at = c(0,1))
mice.lm <- lm(mice.stats[,sig.ids[11]+2] ~ mice.stats[,2])
abline(mice.lm)
