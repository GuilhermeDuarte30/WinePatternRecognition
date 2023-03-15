#Bibliotecas
#install.packages("psych")
library(psych)
#install.packages("corrplot")
library(corrplot)
#install.packages("xlsx")  
library(xlsx)
#install.packages('gpairs')
library(gpairs)
library(cluster)
library(DMwR)

#Ler os dataset
df_total_vinhos<-read.delim(file.choose(),header=TRUE, sep=";", dec=",")
df_total_vinhos
df_total_vinhos<-df_total_vinhos[-c(1)]
df_total_vinhos

#Fixar a semente
set.seed(123)

head(df_total_vinhos)

summary(df_total_vinhos)


#Summary and covariance matrix
#alcool com MSA baixo e impedia que KMO fosse superior a 0,5 e como tal foi removido

data <- df_total_vinhos[,1:10]
summary(data)
round(var(data),2)

#scatterplot
pairs(data, pch = 19, lower.panel = NULL)

#correlation matrix
correlation <- cor(data)
round(correlation,3)

#corrplot
par(oma = c(2, 2, 2, 2)) # space around for text
corrplot.mixed(correlation, 
               order = "hclust", #order of variables
               tl.pos = "lt", #text left + top
               upper = "ellipse")

#Correlation matrix
round(correlation, 3)

#Bartlett test and KMO
#Input is the correlation matrix
cortest.bartlett(correlation)
KMO(correlation)

#Extraction and number of components
# Scale the data
dataZ <- scale(df_total_vinhos[,1:10])
# Assume the number of components (nfactors) = number of variables, i.e., D = 10,
# Always a unrotated extraction
pc10 <- principal(dataZ, nfactors=10, rotate="none", scores=TRUE)  

#pc10
#Eigenvalues - Variances of the principal components 
round(pc10$values,3)

crit_kaiser=data.frame(round(pc10$values,3))

colnames(crit_kaiser) <- c("valor próprio")

#Screeplot - Find the elbow
plot(pc10$values, type = "b", main = "Scree plot for Wine dataset",
     xlab = "Number of PC", ylab = "Eigenvalue")


# Eigenvectors - component loadings
# Each value is the correlation between original variable and component,
# It can be thought as the contribution of each original variable to the component (when squared)
# It summarizes the explained variance too.
pc10$loadings

pc10$communality

#Let's extract a 3 component solution
pc3 <- principal(dataZ, nfactors=3, rotate="none")
pc3

#Let's rotate the 3 components using varimax
pc3r <- principal(dataZ, nfactors=3, rotate="varimax")
pc3r$loadings

round(pc3r$loadings,2)
round(pc3r$communality,2)

#Let's extract a 4 component solution
pc4 <- principal(dataZ, nfactors=4, rotate="none")
pc4

#Let's rotate the 4 components using varimax
pc4r <- principal(dataZ, nfactors=4, rotate="varimax")
pc4r$loadings

round(pc4r$loadings,2)
round(pc4r$communality,2)


#Compute the scores (eigenvalues for each observation)
pc4sc <- principal(dataZ, nfactors=4, rotate="none", scores=TRUE)  
round(pc4sc$scores,3)
mean(pc4sc$scores[,1])
sd(pc4sc$scores[,1])


#Add scores to the data set as new variables
df_total_vinhos$nivel_de_docura<- pc4sc$scores[,1]
df_total_vinhos$densidade<- pc4sc$scores[,2]
df_total_vinhos$nivel_de_acidez<- pc4sc$scores[,3]
df_total_vinhos$nivel_de_protecao_microbiano<- pc4sc$scores[,4]

View(df_total_vinhos)

# Scatterplot
gpairs(df_total_vinhos[14:17])

comp_princ=df_total_vinhos[14:17]

View(comp_princ)

#############################################
# Hierarchical clustering - hclust

# Standardized * Euclidian * Ward.D2 
demodist2 <- dist(comp_princ) # compute distance
hclust_demo2 <- hclust(demodist2,method='ward.D2')
plot(hclust_demo2,hang=-1)

# Cut the dendrogram
groups.k3 <- cutree(hclust_demo2 , k=3) # cut tree into 3 clusters
rect.hclust(hclust_demo2, k=3, border="red") 
aggregate(comp_princ,list(groups.k3), mean)

# Cut the dendrogram
#groups.k4 <- cutree(hclust_demo2 , k=4) # cut tree into 4 clusters
#rect.hclust(hclust_demo2, k=4, border="red") 
#aggregate(comp_princ,list(groups.k4), mean)

# Cut the dendrogram
#groups.k5 <- cutree(hclust_demo2 , k=5) # cut tree into 5 clusters
#rect.hclust(hclust_demo2, k=5, border="red") 
#aggregate(comp_princ,list(groups.k5), mean)
#hierar=data.frame(aggregate(comp_princ,list(groups.k5), mean))

#clusters
hierar=data.frame(aggregate(comp_princ,list(groups.k3), mean))

#Silhouette

plot(silhouette(groups.k3,demodist2))
#plot(silhouette(groups.k4,demodist2))
#plot(silhouette(groups.k5,demodist2))

################################################################################
###K-means
# K-Means=3
kmeans.k3 <- kmeans(comp_princ, 3,nstart=100)
# k = 3 from hclust, nstart = initial random solutions
#attributes(kmeans.k3)  # all elements of the cluster output

#cluster means
kmeans.k3$centers
#kmeans.k3$cluster
kmeans.k3$size
kmeansvinho=table(df_total_vinhos$vinho,kmeans.k3$cluster)
kmeansqualidade=table(df_total_vinhos$quality,kmeans.k3$cluster)

kmeansvinho
kmeansqualidade

#Silhouette
plot(silhouette(kmeans.k3$cluster,demodist2)) 

# K-Means: number of clusters
wssplot <- function(xx, nc=15, seed=1234){
  wss <- (nrow(comp_princ)-1)*sum(apply(comp_princ,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(xx, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(comp_princ, nc=10) 


## Plots
plot(df_total_vinhos$nivel_de_docura, df_total_vinhos$nivel_de_acidez, type="n", xlab="nivel_de_docura", ylab="nivel_de_acidez")
text(x=df_total_vinhos$nivel_de_docura, y=df_total_vinhos$nivel_de_acidez, labels=df_total_vinhos$vinho,col=kmeans.k3$cluster+1)

#plot(df_total_vinhos$nivel_de_docura, df_total_vinhos$densidade, type="n", xlab="nivel_de_docura", ylab="densidade")
#text(x=df_total_vinhos$nivel_de_docura, y=df_total_vinhos$densidade, labels=df_total_vinhos$vinho,col=kmeans.k3$cluster+1)

#plot(df_total_vinhos$nivel_de_docura, df_total_vinhos$nivel_de_protecao_microbiano, type="n", xlab="nivel_de_docura", ylab="nivel_de_protecao_microbiano")
#text(x=df_total_vinhos$nivel_de_docura, y=df_total_vinhos$nivel_de_protecao_microbiano, labels=df_total_vinhos$vinho,col=kmeans.k3$cluster+1)

#plot(df_total_vinhos$nivel_de_acidez, df_total_vinhos$nivel_de_densidade, type="n", xlab="nivel_de_acidez", ylab="nivel_de_densidade")
#text(x=df_total_vinhos$nivel_de_acidez, y=df_total_vinhos$nivel_de_densidade, labels=df_total_vinhos$vinho,col=kmeans.k3$cluster+1)

#plot(df_total_vinhos$nivel_de_acidez, df_total_vinhos$nivel_de_protecao_microbiano, type="n", xlab="nivel_de_acidez", ylab="nivel_de_protecao_microbiano")
#text(x=df_total_vinhos$nivel_de_acidez, y=df_total_vinhos$nivel_de_protecao_microbiano, labels=df_total_vinhos$vinho,col=kmeans.k3$cluster+1)

#######################################################################################################
##PAM clustering

pam.k3 <- pam(comp_princ,3)
#pam.k3 more information using help
table(df_total_vinhos$vinho,pam.k3$clustering)
table(df_total_vinhos$quality,pam.k3$clustering)
plot(df_total_vinhos$nivel_de_docura, df_total_vinhos$nivel_de_acidez, type="n", xlab="nivel_de_docura", ylab="nivel_de_acidez")
text(x=df_total_vinhos$nivel_de_docura, y=df_total_vinhos$nivel_de_acidez, labels=df_total_vinhos$vinho,col=pam.k3$clustering+1)

# PCA and Clustering
clusplot(pam.k3, labels = 3, col.p = pam.k3$clustering)