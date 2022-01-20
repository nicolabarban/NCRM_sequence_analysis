
# ---------------------------------------------------------------------
# Program: SyntaxLab2.r
#  Author:Nicola Barban (n.barban@unibo.it)
#    Date: Jan 21 2022
# --------------------------------------------------------------------
# --------------------------------------------------------------------

#### loading traminer package. 
library(TraMineR)


#add the labels to be used in the graphics
mvad.labels=c("employment","further education", 
              "higher education", "joblessness", "school", "training")

#add abbreviate codes to be used in graphics
mvad.scode=c("EM","FE", "HE", "JL", "SC", "TR")

####################################
#create sequence object WEIGHTED
mvad.seq=seqdef(mvad,17:86,states=mvad.scode, labels=mvad.labels, weights=mvad$weight )
alphabet(mvad.seq)
stlab(mvad.seq)
help(seqdef)


####################################
################## Statistics
##first ten most common sequences
seqtab(mvad.seq)

# transversal statistics
seqstatd(mvad.seq[, 1:8])

# duration list
duration<-seqistatd(mvad.seq)

## transition matrix
mvad.trate = seqtrate(mvad.seq) 
round(mvad.trate, 2)


par(mfrow=c(1,2))
{
  #entropy of state distribution
  Entropy= seqstatd(mvad.seq)$Entropy
  plot(Entropy,main="Entropy of the state distribution ",
       col="black", xlab="Time in months", ylab="Entropy", type="l")
  #sequence turbulence
  Turbulence=seqST(mvad.seq)
  hist(Turbulence, col="grey", main="Histogram of sequence turbulence")
}
par(mfrow=c(1,1))

###plot longitudinal entropy


####################################
# Compute Optimal Matching with transition rates
####################################
submat=seqsubm(mvad.seq, method="TRATE")
dist.om1=seqdist(mvad.seq, method="OM", indel=1, sm=submat)

dist.om1[1:10,1:10]

dist.lcs=seqdist(mvad.seq, method="LCS")
dist.dhd=seqdist(mvad.seq, method="DHD")


##################################################
#Clustering
library(cluster) 
mvad.clusterward = agnes(dist.om1, diss = T, method = "ward") 
plot(mvad.clusterward, ask = F, which.plots = 2)

### Choose the number of clusters
mvad.cl5 = cutree(mvad.clusterward, k = 5)
mvad.cl5[1:10]

## Cluster 1 HE dominated
## Cluster 2 FE dominated
## Cluster 3 Employment Dominated
## Cluster 4 Long  training
## Cluster 5 Joblessness dominated
###########

seqdplot(mvad.seq, group=mvad.cl5, border=NA, space=0)
levels(mvad.cl5)=c(" HE dominated", "FE dominated", "Employment Dominated", "Long  training", "Joblessness dominated")
seqdplot(mvad.seq, group=mvad.cl5, border=NA, space=0)





###medoids
dc= disscenter(dist.om1, group= mvad.cl5, medoids.index="first")
print(mvad.seq[dc,], format = "SPS")

############ Tables
seqtab(mvad.seq[mvad.cl5==1,])
seqtab(mvad.seq[mvad.cl5==2,])



#######logistic regression joblesseness

jobless = mvad.cl5 == 5
jobless.reglog = glm(jobless ~ male + funemp + gcse5eq, family = binomial(link = logit),data = mvad)
summary(jobless.reglog)

###Add clusters to data.frame anc export in stata

library(foreign)
write.dta(mvad, file="mvad_cluster.dta")

##### Multichannel sequence analysis
data(biofam)

## Building one channel per type of event left, children or married
bf <- as.matrix(biofam[, 10:25])
children <-  bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6

## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)

## Using transition rates to compute substitution costs on each channel
mcdist <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
  method="OM", sm =list("TRATE", "TRATE", "TRATE"))

## Using a weight of 2 for children channel and specifying substitution-cost
smatrix <- list()
smatrix[[1]] <- seqsubm(child.seq, method="CONSTANT")
smatrix[[2]] <- seqsubm(marr.seq, method="CONSTANT")
smatrix[[3]] <- seqsubm(left.seq, method="TRATE")
mcdist2 <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
  method="OM", sm =smatrix, cweight=c(2,1,1)) 


####### Measuring alternative cluser solutions

mvad.cl4 <- cutree(mvad.clusterward, k = 4)
mvad.cl6 <- cutree(mvad.clusterward, k = 6)
seqdplot(mvad.seq, group=mvad.cl5, border=NA, space=0)

seqdplot(mvad.seq, group=mvad.cl4, border=NA, space=0)
seqdplot(mvad.seq, group=mvad.cl6, border=NA, space=0)

install.packages(WeightedCluster)
library(WeightedCluster)
####### Statistics for cluster solution
qual<-list()

R2<-c()
ASW<-c()

for(j in 2:10){
  qual[[j]] <- wcClusterQuality(dist.om1, cutree(dist.ward, k = j), weights=mvad$weights)
  R2[j]<-qual[[j]]$stats[7]
  ASW[j]<-qual[[j]]$stats[4]

}

for(j in 2:10){
    ASW[j]<-qual[[j]]$stats[4]
}
par(mfrow=c(1,2))
plot(1:10,ASW, type="l")

wardRange <- as.clustrange(mvad.clusterward, diss = dist.om1, ncluster = 10, weights=mvad$weights)
summary(wardRange, max.rank = 2)
plot(wardRange, stat = c("ASW", "HG", "PBC", "HC"))


