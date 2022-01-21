
# ---------------------------------------------------------------------
# Program: SyntaxLab2.r
#  Author:Nicola Barban (n.barban@unibo.it)
#    Date: Jan 21 2022
# --------------------------------------------------------------------
# --------------------------------------------------------------------

#### loading traminer package. 
install.package("TraMineR")
library(TraMineR)
data(mvad)

#add the labels to be used in the graphics
mvad.labels=c("employment","further education", 
              "higher education", "joblessness", "school", "training")

#add abbreviate codes to be used in graphics
mvad.scode=c("EM","FE", "HE", "JL", "SC", "TR")

####################################
#create sequence object WEIGHTED
library(RColorBrewer)

mvad.seq<-seqdef(mvad,
                15:86,
                states=mvad.scode,
                labels=mvad.labels, 
                weights=mvad$weight,
                cpal=brewer.pal(6,"Set1"))

alphabet(mvad.seq)
stlab(mvad.seq)
help(seqdef)


####################################
################## Statistics
##first ten most common sequences
seqtab(mvad.seq)

# transversal statistics


# duration list
mvad.duration<-seqistatd(mvad.seq)
mvad.duration

EM.duration<-mvad.duration[,1]
SC.duration<-mvad.duration[,5]
JL.duration<-mvad.duration[,4]



seqIplot(mvad.seq)

seqIplot(mvad.seq, sortv=EM.duration)
seqIplot(mvad.seq, sortv=-EM.duration)

seqIplot(mvad.seq, sortv=c(EM.duration,SC.duration,JL.duration))




## transition matrix
mvad.trate = seqtrate(mvad.seq) 
round(mvad.trate, 2)



  #entropy of state distribution
  Entropy= seqstatd(mvad.seq)$Entropy
  plot(Entropy,main="Entropy of the state distribution ",
       col="black", xlab="Time in months", ylab="Entropy", type="l")


  #sequence turbulence
  Turbulence=seqST(mvad.seq)
  hist(Turbulence, col="grey", main="Histogram of sequence turbulence")
seqIplot(mvad.seq, sortv=Turbulence)



####################################
# Compute Optimal Matching with transition rates
####################################

trcost <- seqcost(mvad.seq, method="TRATE")

dist.om1=seqdist(mvad.seq, method="OM", indel=trcost$indel, sm=trcost$sm)





dist.om1[1:10,1:10]

dist.lcs=seqdist(mvad.seq, method="LCS")
dist.dhd=seqdist(mvad.seq, method="DHD")


##################################################
#Clustering
library(cluster) 
mvad.clusterward = agnes(dist.om1, diss = T, method = "ward") 
plot(mvad.clusterward, ask = F, which.plots = 2)

### Choose the number of clusters
mvad.cl4 = cutree(mvad.clusterward, k = 4)
mvad.cl4[1:10]



###########

seqdplot(mvad.seq, group=mvad.cl4, border=NA, space=0)
seqmtplot(mvad.seq, group=mvad.cl4, border=NA, space=0)

levels(mvad.cl4)=c(" EM dominated", "HE dominated", "FE dominated", "Joblessness dominated")
seqdplot(mvad.seq, group=mvad.cl4, border=NA, space=0, main=levels(mvad.cl4))

#### Trying different cluster solution

####### Measuring alternative cluser solutions

mvad.cl5 <- cutree(mvad.clusterward, k = 5)
mvad.cl3 <- cutree(mvad.clusterward, k = 3)
seqdplot(mvad.seq, group=mvad.cl4, border=NA, space=0)
dev.new()
seqdplot(mvad.seq, group=mvad.cl3, border=NA, space=0)
dev.new()

seqdplot(mvad.seq, group=mvad.cl5, border=NA, space=0)

###medoids
dc= disscenter(dist.om1, group= mvad.cl4, medoids.index="first")
print(mvad.seq[dc,], format = "SPS")



############ Tables
seqtab(mvad.seq[mvad.cl4==1,])
seqtab(mvad.seq[mvad.cl4==2,])



#######logistic regression joblesseness

jobless = mvad.cl4 == 4
jobless.reglog = glm(jobless ~ male + funemp + gcse5eq, family = binomial(link = logit),data = mvad)
summary(jobless.reglog)

###Add clusters to data.frame anc export in stata

library(foreign)
write.dta(mvad, file="mvad_cluster.dta")

################################################
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


biofam.clusterward = agnes(mcdist, diss = T, method = "ward") 
plot(biofam.clusterward, ask = F, which.plots = 2)
biofam.cl4 <- cutree(biofam.clusterward, k = 4)

########
seqdplot(child.seq, group=biofam.cl4)
seqdplot(marr.seq, group=biofam.cl4)
seqdplot(left.seq, group=biofam.cl4)


#########

