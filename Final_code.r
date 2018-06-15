#################################################################
#################################################################
#	Code for Project for the Mathematical Modelling module 	#
#	Author:Gianpiero Cea                                   	#
#################################################################
#
#################################################################
# First we install all the required packages.
# This should work but if there are any problems refer
# to the individual documentations
#################################################################
print('Install all required packages')


install.packages("mlbench")
install.packages("Amelia")
install.packages("bnlearn")
install.packages("BDgraph")
install.packages("corrplot")
install.packages("caret")

#################################################################
# Importing all the required libraries
#################################################################
print('#########################################################')
print('Part 0: Importing all libraries that we need')

#Part 0: Importing all libraries that we need:
#mlbench used to import the dataset
library(mlbench)
#used to do imputation of the values
library(Amelia)
#main library to work with Bayesian Network
library(bnlearn)
#possibly used for indirected graphs
library(BDgraph)
#used to plot a correlation plot
library(corrplot)
#used to test accuracy and other metric of models
library(caret)

##################################################################
# Cleaning the data
##################################################################
print('#########################################################')
print('Part 1: Cleaning Data: imputation/omission, visualisation')

# Part 1: Cleaning Data: imputation/omission, visualisation

# save the data in dataframe df
data(BreastCancer)
df <- BreastCancer
# check the missing values, OUTPUT: 16 NA's in 'Bare.nuclei'
summary(df)

#Omitting missing values
df.omit = na.omit(df)

#Imputing missing values
print('Imputing missing values')
## Using Amelia II:
print('Imputing values using Amelia II...')
#to.nom(df) turns a df of factors into a df of numerical values
to.nom <-function(df){
  df2 <- data.matrix(df)
  df2 <- data.frame(df2)
  return(df2)
}
df.num <- to.nom(df)
#to use amelia we can't use factors
res = amelia(df.num)
df.imputed <-df
df.imputed$Bare.nuclei <- as.factor(abs(floor(res$imputations[[1]]$Bare.nuclei)))
#we can have a comparison of our data before and after imputations
print('Comparison of the data before and after Amelia imputations')
summary(df)
summary(df.imputed)

## Using structural.em to impute data, learn from missing data
print('Imputing values using structural.em')
res<-structural.em(df[,2:11],return.all=TRUE)
df.imputed2 <- res$imputed
bn.fit.missing <-res$fitted
bn.missing <- res$dag

# Data correlation
print('Visualising correlation plot')
bc = cor(na.omit(df.num)[,2:11]) #create an object of the features
corrplot.mixed(bc)

#Split the data , p%-(1-p)%
	split <- function(df,p){
	a = floor(nrow(df)*(p/100))
	df.train <- df[1:a,]
	df.test <- df[(a+1):nrow(df),]
	return(list(train =df.train,test =df.test))
	}

df.omit.train <- split(df.omit[,2:11],20)$train
df.omit.test <- split(df.omit[,2:11],20)$test

df.imputed.train <- split(df.imputed[,2:11],20)$train
df.imputed.test <- split(df.imputed[,2:11],20)$test

df.imputed2.train <- split(df.imputed2,20)$train
df.imputed2.test <- split(df.imputed2,20)$test
##################################################################
# Structure learning
#################################################################
print('#########################################################')
print('Structural learning')
# Part 2: Structure Learning

#the following functions applies several algorithms  to the data and returns
#the structure learned
struc.learn <- function(df){
	#Classifiers

	bn.naive.bayes <- naive.bayes(df,"Class")
	bn.tree.bayes <- tree.bayes(df,"Class")

	#Score Based
	bn.hc <-hc(df)
	bn.tabu <-tabu(df)

	#Constraint Based
	#here we give the option parameter test
	bn.gs <- gs(df)
	bn.iamb <-iamb(df)
	bn.fast.iamb <- fast.iamb(df)
	bn.inter.iamb <-inter.iamb(df)
	bn.mmpc <- mmpc(df)
	bn.si.hiton.pc <- si.hiton.pc(df)

	return(list(naive.bayes = bn.naive.bayes,tree.bayes = bn.tree.bayes, hc = bn.hc,tabu=bn.tabu, gs=bn.gs,iamb =bn.iamb,fast.iamb =bn.fast.iamb,mmpc=bn.mmpc,si.hiton.pc=bn.si.hiton.pc))
}

# Here we do strucutural learning for the three different datasets
# ...omit, imputed and imputed2

structures.omit <- struc.learn(df.omit.train)
#we also had the bn learned from the missing data!
structures.omit[[length(structures.omit)+1]] <- bn.missing
print('Structures learnt for omit dataset')

structures.imputed <-struc.learn(df.imputed.train)
print('Structures learnt for imputed dataset')

structures.imputed2 <- struc.learn(df.imputed2.train)
print('Structures learnt for imputed2 dataset')



# We plot all the different structures learned for structures.omit
par(mfrow=c(3,3))
for(i in 1:9){plot(structures.omit[[i]],main=names(structures.omit)[[i]])}
plot(structures.omit[[10]],main=names(structures.omit)[[10]])

# We plot all the different structures learned for structures.imputed
par(mfrow=c(3,3))
for(i in 1:9){plot(structures.imputed[[i]],main=names(structures.imputed)[[i]])}

# We plot all the different structures learned for structures.imputed2
par(mfrow=c(3,3))
for(i in 1:9){plot(structures.imputed2[[i]],main=names(structures.imputed2)[[i]])}

##Compare the hc and tabu graph. OUTPUT: they are the same
compare(structures.imputed2$hc,structures.imputed2$tabu)
##Get a string representing the model:
##OUTPUT:"[Cell.size][Class|Cell.size][Cl.thickness|Class][Cell.shape|Class][Marg.adhesion|Class][Epith.c.size|Class][Bare.nuclei|Class][Bl.cromatin|Class][Normal.nucleoli|Class][Mitoses|Class]"
modelstring(structures.imputed2$tabu)

# Compare graphically the structures
graphviz.compare(structures.omit$naive.bayes,structures.omit$tree.bayes)
##################################################################
# Parameter learning
##################################################################
print('#########################################################')
print('Part 3: Parameter Learning')
# Part 3: Parameter Learning
## Maximum Likelihood Estimate
#we pick an ordering of the variables to assign a direction to the Constraint
#based structures (gs,iamb...)
ordering = node.ordering(structures.imputed2$tree.bayes)
#lapply({function(x)pdag2dag(structures.omit$x,names(df.omit.train))},df.omit.train[4:9])
#bn.fit.omit <-lapply({function(dag)bn.fit(dag,method='mle')},structures.omit)
 gs<-pdag2dag(structures.omit$gs,ordering)
 iamb<-pdag2dag(structures.omit$iamb,ordering)
 fast.iamb<-pdag2dag(structures.omit$fast.iamb,ordering)
 mmpc<-pdag2dag(structures.omit$mmpc,ordering)
 si.hiton.pc<-pdag2dag(structures.omit$si.hiton.pc,ordering)

structures.omit$gs = gs
structures.omit$iamb = iamb
structures.omit$fast.iamb = fast.iamb
structures.omit$mmpc = mmpc
structures.omit$si.hiton.pc = si.hiton.pc

structures.imputed2$gs = gs
structures.imputed2$iamb = iamb
structures.imputed2$fast.iamb = fast.iamb
structures.imputed2$mmpc = mmpc
structures.imputed2$si.hiton.pc = si.hiton.pc

structures.imputed$gs = gs
structures.imputed$iamb = iamb
structures.imputed$fast.iamb = fast.iamb
structures.imputed$mmpc = mmpc
structures.imputed$si.hiton.pc = si.hiton.pc

#Parameter Learning for omit data
bn.fit.omit <-lapply(structures.omit,{function(dag)bn.fit(dag,method='mle',data=df.omit.train)})
bn.fit.omit.bayes <-lapply(structures.omit,{function(dag)bn.fit(dag,method='bayes',data=df.omit.train)})
#Parameter Learning for imputed data
bn.fit.imputed <-lapply(structures.imputed,{function(dag)bn.fit(dag,method='mle',data=df.imputed.train)})
bn.fit.imputed.bayes <-lapply(structures.imputed,{function(dag)bn.fit(dag,method='bayes',data=df.imputed.train)})
#Parameter Learning for imputed2 data
bn.fit.imputed2 <-lapply(structures.imputed2,{function(dag)bn.fit(dag,method='mle',data=df.imputed2.train)})
bn.fit.imputed2.bayes <-lapply(structures.imputed2,{function(dag)bn.fit(dag,method='bayes',data=df.imputed2.train)})

#Visualise as an example one of the node Conditional probability table:
bn.fit.barchart(bn.fit.imputed2$tabu$Cl.thickness)

##################################################################
# Accuracy comparison & Inference
##################################################################
print('#########################################################')
print('Part 4: Accuracy comparison & Inference')
# Part 4: Accuracy comparison & Inference
# gets the predictions on the three kinds of test sets (omit,imputed, imputed2)
predicted.omit <- lapply(bn.fit.omit,{function(fit)predict(fit, node = "Class", data = df.omit.test)})
predicted.omit.bayes <- lapply(bn.fit.omit.bayes,{function(fit)predict(fit, node = "Class", data = df.omit.test)})

predicted.imputed <- lapply(bn.fit.imputed,{function(fit)predict(fit, node = "Class", data = df.imputed.test)})
predicted.imputed.bayes <- lapply(bn.fit.imputed.bayes,{function(fit)predict(fit, node = "Class", data = df.imputed.test)})

predicted.imputed2 <- lapply(bn.fit.imputed,{function(fit)predict(fit, node = "Class", data = df.imputed.test)})
predicted.imputed2.bayes <- lapply(bn.fit.imputed.bayes,{function(fit)predict(fit, node = "Class", data = df.imputed2.test)})

#print everything
##############################################################################
print('imputed')
for(i in 1:9){
	print(names(predicted.imputed)[[i]])
	print(confusionMatrix(unlist(predicted.imputed[i]),unlist(df.imputed.test$Class)))
}
for(i in 1:9){
	print(names(predicted.imputed.bayes)[[i]])
	print(confusionMatrix(unlist(predicted.imputed.bayes[i]),unlist(df.imputed.test$Class)))
}
##############################################################################
print('omit')
for(i in 1:10){
	print(names(predicted.omit)[[i]])
	print(confusionMatrix(unlist(predicted.omit[i]),unlist(df.omit.test$Class)))
}
for(i in 1:10){
	print(names(predicted.omit.bayes)[[i]])
	print(confusionMatrix(unlist(predicted.omit.bayes[i]),unlist(df.omit.test$Class)))
}
##############################################################################
print('imputed2')
for(i in 1:9){
	print(names(predicted.imputed2)[[i]])
	print(confusionMatrix(unlist(predicted.imputed2[i]),unlist(df.imputed2.test$Class)))
}
for(i in 1:9){
	print(names(predicted.imputed2.bayes)[[i]])
	print(confusionMatrix(unlist(predicted.imputed2.bayes[i]),unlist(df.imputed2.test$Class)))
}

#Approximate inference

#example: What is the probability that the tumor is malignant if we a Cell.size
#of 4 and Bare.nuclei of 2

#we need to repeat since the function cpquery is stochastic

l<-vector("list", 2000)
for( i in 1:2000)
{
l[i] <- cpquery(bn.fit.imputed$naive.bayes,event=(Class == "malignant"),evidence=(Cell.size=="4")&(Bare.nuclei=="2"))
}

mean(as.numeric(l))
############################################################################
# UNDIRECTED GRAPH:

library(BDgraph)
#problem with bigger iter
 res<-bdgraph(df.imputed2,method ="gcgm", iter= 100, print = 10)
summary(res)

#we can also learn another structure using Marginal Pseudo Likelihhod
res2<-bdgraph.mpl(df.imputed2,iter= 1000)
summary(res2)

#################################################################
