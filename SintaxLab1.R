# ---------------------------------------------------------------------
#  Program: SyntaxLab1.r
#  Author:Nicola Barban (n.barban@unibo.it)
#  Date: Jan 20 2022
#  NCRM Introduction to Sequence Analysis for Social Sciences
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# --------------------- PART 1: Refreshing R -------------------------
# --------------------------------------------------------------------

# Arithmetic with two numbers
2+2

2-4

2*3

8/3


# In R you can make your own functions.
mySqrt <- function(x) { x^.5 }

# Now we can use it
mySqrt(9)
# ----

# In R you can make your own functions.
mySqrt <- function(x) { x^.5 }

# Now we can use it
mySqrt(9)

# ---

# Let's make a variable called Fred
Fred <- 2+2
Fred

# One variable can be assigned to another.
a <- Fred

# Now 'a' is a copy of Fred.
a

# We read the next line as "A gets square root of 'a' times five"
A <- sqrt(a)*5
A

# This one is a little tricky.  We can re-assign A.
A <- A*2
A

# Here, we create a vector of length 13; notice what ":" does
x2 <- 1:13
x2

# Vectors can also be created by the "c()" function.
# "c" stands for "concatenate".
myvec <- c(8,13,2,1,6)
myvec
# ---


# ---------------------------------------------------------------------
#  READING AND WRITING DATA ON THE DISK

# First, we need to set the Project Folder (Working Directory) for R.
# You can either set the Project Folder in the GUI or in a script.
#

# To display the Project Folder
getwd()


# To list all files in your working directory
list.files()         



# We will use "csv" files (comma separated values) for this workshop.  
#  The csv format is the most general format for data files.s
#  You can save csv files from SPSS, SAS, or Excel.  
#  And you can even open it up in a text editor.

# Now we will read our first data file.

data_DHS <- read.csv(file="http://nicolabarban.com/NCRM_sequence_analysis/DataLab1.csv", header=TRUE)

# We have created a "data.frame" called mvad.
is.data.frame(data_DHS)

# Summary is one way of finding out about any object
summary(data_DHS)



# ---------------------------------------------------------------------
#   SUBSETTING & INDEXING

# We can select a variable from mvad by writing:
id <- data_DHS$filenw

# This selects the third element of zygo by using an index.
id[3]

# We could also select the third row and first column of mvad
#   and obtain the same cell of the dataframe.
data_DHS[3,1]

# If we omit the number after the comma, we obtain all columns.
data_DHS[3, ]

# ---------------------------------------------------------------------
# Installing Packages


### sometimes we need functions that are not in the base packages.  

### First, we need to install them in the system

install.packages("foreign")
library(foreign)  ### also require(namepackage)

### foreign is a package that can be used to import dataset from other statistical packages,e.g. Stata and SPSS 
help(package="foreign")



# --------------------------------------------------------------------
# --------------------- PART 2: Setting up a sequence analysis -------
# --------------------------------------------------------------------

#### Installing traminer 
install.packages("TraMineR")
library(TraMineR)


#### mvad data
data(mvad)
help(mvad)
names(mvad)

mvad[1:3,]
t1=table(mvad$gcse5eq,mvad$catholic)
t1
prop.table(t1, 1)
prop.table(t1, 2)

#-----

data(biofam)
help(biofam)
names(biofam)
biofam$age=2002-biofam$birthyr

min(biofam$age)
max(biofam$age)

min(biofam$age[biofam$sex=="woman"])
max(biofam$age[biofam$sex=="woman"])

t2=table(biofam$p02r04, biofam$sex)
prop.table(t2,2)

################# Sequence analysis in R

library(TraMineR)
data(mvad)

#add the labels to be used in the graphics
mvad.labels=c("employment","further education", 
"higher education", "joblessness", "school", "training")

#add abbreviate codes to be used in graphics

mvad.scode=c("EM","FE", "HE", "JL", "SC", "TR")

####################################
#create sequence object
mvad.seq=seqdef(mvad,17:86,states=mvad.scode, labels=mvad.labels)
alphabet(mvad.seq)
stlab(mvad.seq)
help(seqdef)

####################################
###change sequence representation
mvad.seq[1:3, ] 

 print(mvad.seq[1:3, ], format = "SPS")
 # distinte state representation
seqdss(mvad.seq[1:3, ])

####################################
####plot distribution and freq plot

#plot distribution plot
seqdplot(mvad.seq)
#plot frequencies plot
seqfplot(mvad.seq)



#place 3 plots plus legend in a figure
par(mfrow=c(2,2))
{
#plot first ten sequences
seqiplot(mvad.seq,
    with.legend=F, 
    border=NA ,
    space=0, 
    main="index plot (first ten sequences)" )
# ten  most frequent sequences
seqfplot(mvad.seq,
    with.legend=F,
    border=NA ,
    space=0,
    pbarw=T, 
    main="Sequence frequency plot" )
#state distribution by time
seqdplot(mvad.seq,
    with.legend=F,
    border=NA ,
    space=0,
    main="State distribution plot" )
#legend as separate graphic

seqlegend(mvad.seq,
     cex=.75)
}
par(mfrow=c(1,1))



library(RColorBrewer)
#colors in the graph
display.brewer.all()
cpal(mvad.seq)=brewer.pal(6,"Greys")
seqdplot(mvad.seq)
colours()
cpal(mvad.seq)=brewer.pal(6,"Accent")


#### Bivariate graphs
levels(mvad$male)=c("Women", "Men")
seqdplot(mvad.seq,
 group=mvad$male,
 border=NA ,
    space=0)

#average time spent in each state
seqmtplot(mvad.seq, group=mvad$male)

#modal state state
seqmsplot(mvad.seq, group=mvad$male)

#index plot
seqIplot(mvad.seq)

####################################
################## Statistics
##first ten most common sequences
seqtab(mvad.seq)

# transversal statistics
seqstatd(mvad.seq[, 1:8])

# number of transitions
 seqtransn(mvad.seq[1:10,]) 

## transition matrix
 mvad.trate = seqtrate(mvad.seq) 
  round(mvad.trate, 2)


###plot for frequent subsequences and plot first 15 most frequent
mvad.seqe=seqecreate(mvad.seq)
fsubseq=seqefsub(mvad.seqe, pmin.support=0.05,)
plot(fsubseq[1:15], col="green")


####################################
# Compute Optimal Matching with transition rates
####################################
submat=seqsubm(mvad.seq, method="TRATE")
dist.om1=seqdist(mvad.seq, method="OM", indel=1, sm=submat)

dist.om1[1:10,1:10]

############################################################################################################
