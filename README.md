---
title: "Viper Fang Length Evolution"
author: "Matthew Holding"
date: "5/30/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

# Analyses used to test biological correlates of fang length in Viperidae

### Load necessary R packages

```{r, eval=TRUE, message=FALSE, results='hide'}
#Statistics
library(nlme)
library(ape)
library(phytools)
library(geiger)
library(rptR)
library(picante)
library(car)
library(tidyr)
library(robCompositions)
library(zCompositions)
library(compositions)
library(tidyr)

#environmental data
library(rgbif)
library(tidyverse)
library(raster)
library(sp)
library(elevatr)
library(rgdal)

#data visualization
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(factoextra)

```

### Process species' mean morphological dataset and accompanying phylogenetic tree
```{r, message=FALSE,results='hide', warning=FALSE}
#Read in species mean fang and other morphological data for each species with species and rownames for tree pruning
data <- read.table(file = "~/Desktop/FangData/FangDataFinal.txt", row.names=1,header = TRUE, sep = "\t")

data$Species <- rownames(data)

#log transform head and fang lengths
data$logHL <- log(data$HL)
data$logMFL <- log(data$MFL)

#Read name-edited version of the ultrametric Viperidae species tree from Alencar et al. 2016
tree <- read.newick(file="~/Desktop/FangData/Vipertree_named.tre")
check<-name.check(tree, data)
check
#prune tree to contain only tips for which we have fang data
subtree <- drop.tip(tree, check$tree_not_data)
#correct branch length rounding errors leading to failure of is.ultrametric()
ult_tree<-force.ultrametric(subtree,method="extend")
```

### Calculate allometric slopes for head length vs. fang length and plot (shown as Fig 1A)
```{r}
#Run phytools phyl.RMA() to get phylogenetic reduced major axis regression slope
y <- data$logMFL
names(y) <- rownames(data)
x <- data$logHL
names(x) <- rownames(data)
phyrma <- phyl.RMA(x = x, y=y, tree = ult_tree, method = "lambda")
phyrma

#Run pgls() function implemented in nlme package to get a pgls regression slope
#First, fit Pagel's lambda correlation structure with pruned species tree
pagel <-corPagel(1, phy=ult_tree, fixed = FALSE, form =~Species)

#run pgls and summarize to extract model estimates
pglsModelpa <- gls(logMFL ~ logHL , correlation = pagel,
                   data = data, method = "ML")
summary(pglsModelpa)


#plot relationship between head length and fang length and add lines with phylogenetically-corrected slopes
plot(data$logHL,data$logMFL, xlim=c(2.5,4.5), ylim=c(1,3.3), xlab="log Head Length (mm)",ylab="log Fang Length (mm)")

#write phyrma.RMA equation
phyRMA_eq = function(x){-2.992758 + ( 1.472002 *x)}
par(new=TRUE)
plot(phyRMA_eq(1:10), type='l', col = "red", lty=3, lwd=3, xlim=c(2.5,4.5), ylim=c(1,3.3),ylab="", xlab="")

#write regression equation for plotting regression line from nlme pgls
pgls_eq = function(x){-2.594213 + ( 1.353454 *x)}
par(new=TRUE)
plot(pgls_eq(1:10), type='l', col = "blue", lty=3, lwd=3,  xlim=c(2.5,4.5), ylim=c(1,3.3),ylab="", xlab="")

#Add a line with a slope=1 to emphase departure from isometry
segments(2.40,0.90, 5.40, 3.90)
```


```{r}
#Read in dataset of individual, adult snake head and fang lengths
data.ind <- read.table("~/Desktop/FangData/Adults_allometry.tsv", sep = "\t", header = TRUE)
data.ind$logMFL <- log(data.ind$MFL)
data.ind$logHL <- log(data.ind$HLmm)

### Fit models of head length vs. fang length for each species with 10 or more individuals measured (Fig 1B)

#fit relationships between log fang length and log head length for each species
fits <- nlme::lmList(logMFL ~ logHL | CORRECT_scientificName, data=data.ind, )
fits <- summary(fits); fits <- as.data.frame(fits$coefficients[,,2])
#extract slope estimate for plotting
slope <- fits$Estimate

#Make histogram of species' allometric coefficients + evolutionary allometric coefficients
hist(slope, breaks = 20, col = NULL, main = NULL, xlab = "Allometric Coefficient")
abline(v=mean(slope),col="black", lty=3,lwd=2)
n <- length(slope)
margin <- qt(0.95,df=n-1)*sd(slope)/sqrt(n)
lower<-mean(slope)-margin
upper<-mean(slope)+margin
x <- c(lower, lower, upper, upper)
y <- c(0,6.2,6.2,0)
polygon(x,y, col=NULL, border = "darkgrey", lty = 3, lwd=2)
B.pgls <- 1.35
B.RMA <- 1.47
abline(v=B.pgls,col="blue", lty=3, lwd=2)
abline(v=B.RMA,col="red", lty=3, lwd=2)
```


### Species' repeatability of relative fang length metric 

```{r}

data.ind <- read.table("~/Desktop/FangData/Adults_allometry.tsv", sep = "\t", header = TRUE)
#Repeatability by boostrapping the dataset 10,000x
rpt(FL_HL ~ (1|CORRECT_scientificName), grname = c("CORRECT_scientificName"), data = data.ind, datatype = c("Gaussian"),
              CI = 0.95, nboot = 10000, npermut = 1000, parallel = TRUE,
              ncores = 7, ratio = TRUE, adjusted = TRUE, expect = "meanobs",
              rptObj = NULL, update = FALSE)
```

### Calculate phylogenetic signal in relative fang length evolution across Viperidae and plot in phenogram

```{r}

#make a named vector of relative fang length data:
RelFangL <- data[,4] #ratio
names(RelFangL) <- data$Species
RelFangL <- RelFangL[order(factor(names(RelFangL), levels = factor(ult_tree$tip.label)))]

#calculate Pagel's lambda and Blomberg's K for relative fang length and test significance
phylosig(ult_tree, RelFangL, method = "lambda", test = T)
phylosig(ult_tree, RelFangL, method = "K", test = T)

#plot phenogram of relative fang length evolution over time
phenogram(ult_tree,RelFangL,spread.labels=TRUE,spread.cost=c(1,0))

```

### Continuous trait mapping of relative fang length to the Viperidae phylogenetic tree
```{r}
obj<-contMap(ult_tree,RelFangL)
#invert colors so red is larger values and blue is smaller
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}
plot(setMap(obj,invert=TRUE))
```


### Get localities from GBIF database for each species in the dataset 

```{r, eval=FALSE}
splist <- read_delim("~/Desktop/FangPaper_Clean/LatRecs/SpList_all.txt", delim = "\t")%>%
  unlist(use.names = F)

#get taxon keys
keys <- lapply(as.list(splist), function(x) name_suggest(x, rank = "species"))%>%
  bind_rows()%>%
  filter(canonicalName %in% splist)%>%
  dplyr::select(key)%>%
  unlist(use.names = F)

#obtain records
sp.occ <- occ_search(taxonKey = keys, limit=10000, return = 'data', hasGeospatialIssue = FALSE, fields = c('scientificName', 'decimalLatitude','decimalLongitude'), hasCoordinate = T) # remember to change the limit
bomb.out <- sp.occ[sapply(sp.occ, "length", USE.NAMES = F) != 1]%>%
  bind_rows()

df <- do.call(rbind, sp.occ)
#write.table(df, file = "~/Desktop/FangPaper_Clean/LatRecs/LatRecords.txt", quote = FALSE, sep = "\t")
```

```{r}
library(raster)
library(sp)
library(elevatr)
library(rgdal)

#Pull down 2.5 resolution WorldClim environmental data
r <- getData("worldclim",var="bio",res=2.5)
r <- r[[c(1:19)]]
names(r) <- c(1:19)

#Associate species' localities with environmental data rasters
data <- read.table(file = "~/Desktop/FangData/LatRecords.txt", header = TRUE, sep = "\t", row.names = NULL)
lats <- data$decimalLatitude
lons <- data$decimalLongitude

coords <- data.frame(x=lons,y=lats)
coords <- coords[c(1:13940,13942:76315,76317:length(coords$x)),] #remove two bad points in the ocean
points <- SpatialPoints(coords, proj4string = r@crs)

values <- raster::extract(r,points)

#Filter raster data set to exclude some non-sensical values (ocean points, funny species' names)
pcadf <- as.data.frame(values)
pcadf <- cbind(data[c(1:13940,13942:76315,76317:length(data$decimalLongitude)),], pcadf)
pcadf <- pcadf[complete.cases(pcadf), ]
good.names <- unique(pcadf$Species)
good.names <- good.names[c(1:45,73:331)]
pcadf <- pcadf[is.element(pcadf$Species,good.names),]

#conduct environmental principal components analysis (thousands of snake localities and 19 BioClim variables)
pca <- princomp(pcadf[,4:22])

pc1 <- pca$scores[,1]
pc2 <- pca$scores[,2]
pc3 <- pca$scores[,3]

#write out PCA loadings
#dir.create("~/Desktop/FangData/OutputTables")
#write.table(loadings(pca,cutoff=0),"~/Desktop/FangData/OutputTables/env_pca_loadings.txt", quote=FALSE,sep = "\t")

#aggregate data by species for downstream analysis, extracting medians for each data type
dfp <- cbind.data.frame(pcadf,pc1, pc2, pc3)
all_med <- aggregate(.~ Species, data = dfp, median)
#write.table(all_med,"~/Desktop/FangData/OutputTables/envir_values_median.txt",row.names = FALSE,quote=FALSE,sep = "\t") #Tabled pared down to the focal species for the study in Excel and added to "FangDataFinal.txt"

```

### Fit models predicting mean fang length (MFL) of 199 species full dataset

```{r}
data <- read.table(file = "~/Desktop/FangData/FangDataFinal.txt", row.names=1,header = TRUE, sep = "\t")
data$logMFL <- log(data$MFL)
data$logHL <- log(data$HL)

#Remove rows with NAs
data <- data[complete.cases(data), ]
data$Species <- rownames(data)

env_data <- cbind(data$Temp,data$Prec,data$PC1)
colnames(env_data) <- c("Temperature", "Precipitation","environPC1")
pairs(env_data)

treepath <- "~/Desktop/FangPaper_Clean/Posterior_trees_vipers/Test/"
treelist <- dir(treepath, pattern ="tree")
```

#### Note that below this point code is written to be run on a single 'pglsModel' at a time while looping through a posterior set of 200 trees from Alencar et al. 2016. This portion of code *will not run* with a single click. Instead, specify pglsModel anew each time.

``` {r, eval=FALSE}
#model list:
#pglsModel = gls(MFL ~  HL, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Prec, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + PC1, correlation = corr, data = data, method = "ML")

out.file <- ""
HL.out <- ""
temp.out <- ""
PrecOrPC.out <- ""
AIC.out <- ""
tree.out <- ""
lambda.out <- ""
resid.out <- ""
for(i in 1:length(treelist)) {
  tryCatch({
    tree <- read.newick(file = paste(c(treepath,treelist[i]), collapse = ""))
    check<-name.check(tree, data)
    subtree <- drop.tip(tree, check$tree_not_data)
    ult_tree<-force.ultrametric(subtree,method="extend")
    ult_tree$edge.length <- ult_tree$edge.length * 100
    corr<-corPagel(1, phy=ult_tree, fixed = FALSE, form = ~Species)
    #Edit model each time loop is run
    pglsModel = gls(MFL ~  HL + Temp + PC1 , correlation = corr, data = data, method = "ML")
    cap <- capture.output(summary(pglsModel))
    lambda <- cap[11]
    PC <- cap[18]
    temp <- cap[17]
    HL <- cap[16]
    int <- cap[19]
    AIC <- cap[5]
    resid <- as.data.frame(pglsModel$residuals)
    temp.out <- rbind(temp.out, temp)
    HL.out <- rbind(HL.out, HL)
    PrecOrPC.out <- rbind(PrecOrPC.out, PC)
    lambda.out <- rbind(lambda.out, lambda)
    AIC.out <- rbind(AIC.out, AIC)
    tree.out <- rbind(tree.out, treelist[i])
    resid.out <- cbind(resid.out, resid)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
}

#Collect and collate model results using only one of the following 2-line snippets of code
#choose based on the variables in the model


out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Prec", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "PC1", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "AIC")

colnames(out.file) <- cols
out.file <- out.file[2:nrow(out.file), ]

#change outfile name as appropriate:
#write.table(out.file, file = "~/Desktop/NewAIC/HL_TempPC1.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

AICs <- cbind(tree.out, AIC.out)
cols <- c("Tree","AIC")
colnames(AICs) <- cols
AICs <- as.data.frame(AICs)
AICs <- AICs[2:length(AICs$AIC),]
AICs <- separate(AICs, col = "AIC", into = c("BL","BL2","AIC"), sep = "\\s")
AICs <- AICs[,c(1,4)]
#write.table(AICs, file = "~/Desktop/NewAIC/AIC_HL_TempPC1_Int.txt", row.names = FALSE,  col.names = TRUE, quote = FALSE, sep = "\t")
            
```

### Compare models of environmental variables vs. fang length
```{r}
#Reread and then merge all data AIC frames 
df1 <- read.table("~/Desktop/NewAIC/AIC_HLonly.txt",header = TRUE)
colnames(df1) <- c("Tree", "HL")
df2 <- read.table("~/Desktop/NewAIC/AIC_HL_Temp.txt",header = TRUE)
colnames(df2) <- c("Tree", "HL_Temp")
df3 <- read.table("~/Desktop/NewAIC/AIC_HL_Prec.txt",header = TRUE)
colnames(df3) <- c("Tree", "HL_Prec")
df4 <- read.table("~/Desktop/NewAIC/AIC_HL_PC1.txt",header = TRUE)
colnames(df4) <- c("Tree", "HL_PC1")


df_list <- list(df1, df2, df3, df4)

AIC_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  
AIC_merged <- AIC_merged[complete.cases(AIC_merged), ]

df <- ""
for(i in 1:length(AIC_merged$Tree)){
vals <- AIC_merged[i,]
vals <- as.numeric(vals[,2:5])
names(vals) <- colnames(AIC_merged[,2:5])
weights <- aicw(na.exclude(vals))
weights$model <- rownames(weights)
df <- rbind(df,weights)
}
df <- df[2:nrow(df),]
colnames(df) <- c("AIC","delta","model.weight","model")
df$model.weight <- as.numeric(df$model.weight)
df$model <- factor(df$model, levels = c("HL", "HL_Temp" ,  "HL_Prec" , "HL_PC1"))
ggboxplot(df, x = "model", y="model.weight", add = "jitter") + rotate_x_text(45)



#Make plot showing the results of the top model: MFL ~ HL + Temp
library(RColorBrewer)
cols = brewer.pal(3, "RdPu")
# Define colour pallete
# Use the following line with RColorBrewer
pal = colorRampPalette(cols)

# Rank variable for colour assignment
data$order = findInterval(data$Temp, sort(data$Temp))

# Make plot using pgls model of best tree from Alencar et al. 2016 for slope as defined above
corr<-corPagel(1, phy=ult_tree, fixed = FALSE, form = ~Species)
pglsModelpa <- gls(MFL ~ HL, correlation = corr, data=data, method = "ML")
plot(MFL ~ HL, data=data, pch=19, col=pal(nrow(data))[data$order],xlab="Species Mean Head Length (mm)", ylab = "Species Mean Fang Length (mm)")
#add pgls best fit line for HL vs MFL
abline(pglsModelpa, col = "black")

#Extract phylogenetic residuals from head length vs fang length to explore explanatory power of temperature
pglsModel_reduced <- gls(MFL ~ HL , correlation = corr,
                   data = data, method = "ML")

resid <- as.data.frame(pglsModel_reduced$residuals)
Temp <- data$Temp/10
resid_df <- cbind(resid,Temp)
colnames(resid_df)<-c("Residuals","Temperature")
resid_df$Species <- rownames(resid_df)

ols <- lm(Residuals ~ Temperature, data=resid_df)
plot(Residuals~I(Temperature/10), data=resid_df, xlab = "Temp (C)",ylab = "Residual Fang Length")
abline(ols, col = "black")
summary(ols)

```


### Analyzing models of diet-effects on fang length
#### Read in the data
```{r}
data1 <- read.table(file = "~/Desktop/FangData/FangDataFinal.txt", row.names=1,header = TRUE, sep = "\t")
data2 <- read.table(file = "~/Desktop/FangData/DietFang.txt", row.names=1,header = TRUE, sep = "\t")

data_merge <- merge(data2, data1,
                          by = 'row.names', all = TRUE)
#Remove rows with NAs
data <- data_merge[complete.cases(data_merge), ]
rownames(data) <- data$Row.names
colnames(data)[1] <- "Species"

data <- data[c(1,3:57,59:71),] #no fish eaters

diet_data <- data[,2:6] #no fish
diet_data_noZero <- cmultRepl(diet_data)
colnames(diet_data_noZero) <- c("mammalnz", "birdnz","reptilenz","amphibiannz","invertnz")
data <- cbind(data,diet_data_noZero)
clr_diet <- clr(diet_data_noZero)
colnames(clr_diet) <- c("clrmamm", "clrbird", "clrreptile","clramphibian","clrinvert")
data <- cbind(data,clr_diet)

#assess interclass diet item correlations after centered-log-ratio transformation
pairs(clr_diet, panel=function(x,y){
  points(x,y, pch=19)
  abline(lm(y~x), col='black', lty=3)
})


#create logit-transformed variables appropriate for linear modeling
data$logitMamm2 <- logit(data$mammalnz)
data$logitRep <- logit(data$reptilenz)
data$logitAmph <- logit(data$amphibiannz)
data$logitInv <- logit(data$invertnz)
data$logitMB <- logit(data$mammalnz + data$birdnz)
data$logitSoft <- logit(data$mammalnz +data$birdnz + data$amphibiannz)
data$logitHold <- logit(data$reptilenz  + data$amphibiannz + data$birdnz)
data$logitNotAmph <- logit(data$reptilenz + data$invertnz + data$mammalnz + data$birdnz)

data$Temp <- data$Temp/10 #get to normal units of celcius

#robust PCA for compositional data to extract principal components of dietary space
rpca <- pcaCoDa(diet_data_noZero)
summary(rpca)
#write.table(scores(rpca), "~/Desktop/NewSupps/DietPCA_scores.tsv", sep = "\t", quote = FALSE)
#write.table(loadings(rpca, cutoff=0), "~/Desktop/NewSupps/DietPCA_loadings.tsv", sep = "\t", quote = FALSE)
plot(rpca)
biplot(rpca)
dpc1 <- rpca$scores[,1]
dpc2 <- rpca$scores[,2]
data$dpc1 <- dpc1
data$dpc2 <- dpc2


#Get pruned tree
tree <- read.newick(file="~/Desktop/FangData/Vipertree_named.tre")
check<-name.check(tree, data)
check
subtree <- drop.tip(tree, check$tree_not_data)
ult_tree<-force.ultrametric(subtree,method="extend") #tree only not ultrametric due to rounding error
ult_tree$edge.length <- ult_tree$edge.length*100

```




### Run pgls models combining head length, temperature, and alternative "diet" constructs

```{r, eval=FALSE}

#substitute the following models 1 by 1 into the loop below
pgls_head_only <- gls(MFL ~ HL , correlation = corr,
                   data = data, method = "ML")
pgls_temp <- gls(MFL ~ HL + Temp, correlation = corr,
                      data = data, method = "ML")
pgls_mammal <- gls(MFL ~ HL + Temp + logitMamm , correlation = corr,
                      data = data, method = "ML")
pgls_reptile <- gls(MFL ~ HL + Temp + logitRep , correlation = corr,
                      data = data, method = "ML")
pgls_amphib <- gls(MFL ~ HL + Temp + logitAmph , correlation = corr,
                      data = data, method = "ML")
pgls_dpc1 <- gls(MFL ~ HL + Temp + dpc1  , correlation = corr,
                data = data, method = "ML")
pgls_dpc2 <- gls(MFL ~ HL + Temp + dpc2 , correlation = corr,
                 data = data, method = "ML")
pgls_soft <- gls(MFL ~ HL + Temp + logitSoft , correlation = corr,
                      data = data, method = "ML")

treepath <- "~/Desktop/FangPaper_Clean/Posterior_trees_vipers/Test/"
treelist <- dir(treepath, pattern ="tree")

out.file <- ""
HL.out <- ""
temp.out <- ""
diet.out <- ""
AIC.out <- ""
tree.out <- ""
lambda.out <- ""
resid.out <- ""
for(i in 1:length(treelist)) {
  tryCatch({
    tree <- read.newick(file = paste(c(treepath,treelist[i]), collapse = ""))
    check<-name.check(tree, data)
    subtree <- drop.tip(tree, check$tree_not_data)
    ult_tree<-force.ultrametric(subtree,method="extend")
    ult_tree$edge.length <- ult_tree$edge.length * 100
    corr<-corPagel(1, phy=ult_tree, fixed = FALSE, form = ~Species)
    ####Edit model each time loop is run#####
    pglsModel <- gls(MFL ~ HL + Temp + logitMamm , correlation = corr,
                     data = data, method = "ML")
    cap <- capture.output(summary(pglsModel))
    #lambda <- cap[11]
    diet <- cap[25]
    temp <- cap[24]
    HL <- cap[23]
    AIC <- MuMIn::AICc(pglsModel)
    resid <- as.data.frame(pglsModel$residuals)
    temp.out <- rbind(temp.out, temp)
    HL.out <- rbind(HL.out, HL)
    diet.out <- rbind(diet.out, diet)
    #lambda.out <- rbind(lambda.out, lambda)
    AIC.out <- rbind(AIC.out, AIC)
    tree.out <- rbind(tree.out, treelist[i])
    resid.out <- cbind(resid.out, resid)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
}

length(tree.out)
out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, diet.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "Diet", "AIC")
colnames(out.file) <- cols
out.file <- out.file[2:nrow(out.file), ]

#change outfile name as appropriate:
#write.table(out.file, file = "~/Desktop/NewAIC/NewDiet/HL_Temp_Mammals.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

AICs <- cbind(tree.out, AIC.out)
cols <- c("Tree","AIC")
colnames(AICs) <- cols
AICs <- as.data.frame(AICs)
AICs <- AICs[2:length(AICs$AIC),]

#write.table(AICs, file = "~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_Mammal.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

```

### Process AICc-based model comparisons to calculate Akaike weights (omega)

```{r}
#Reread and then merge all data AIC frames 
df1 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_only.txt",header = TRUE)
colnames(df1) <- c("Tree", "HL")
df2 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp.txt",header = TRUE)
colnames(df2) <- c("Tree", "HL_Temp")
df3 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_Mammal.txt",header = TRUE)
colnames(df3) <- c("Tree", "HL_Temp_Mammal")
df4 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_Amphbian.txt",header = TRUE)
colnames(df4) <- c("Tree", "HL_Temp_Amphib")
df5 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_Reptile.txt",header = TRUE)
colnames(df5) <- c("Tree", "HL_Temp_Reptile")
df6 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_softprey.txt",header = TRUE)
colnames(df6) <- c("Tree", "HL_Temp_Soft")
df7 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_dietPC1.txt",header = TRUE)
colnames(df7) <- c("Tree", "HL_Temp_dPC1")
df8 <- read.table("~/Desktop/NewAIC/NewDiet/AIC_HL_Temp_dietPC2.txt",header = TRUE)
colnames(df8) <- c("Tree", "HL_Temp_dPC2")

df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8)

AIC_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  

AIC_final <- subset(AIC_merged, !is.na(AIC_merged[,4]) == TRUE) #Get rows where model ran with the tree and mammal model
AIC_final <- subset(AIC_final, rowSums(!is.na(AIC_final[,2:9])) > 4) #Get trees where more than half the models converged

#write.table(AIC_final, file = "~/Desktop/NewAIC/NewDiet/AIC_merged.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

df <- ""
for(i in 1:length(AIC_final$Tree)){
vals <- AIC_final[i,]
vals <- as.numeric(vals[,2:9])
names(vals) <- colnames(AIC_final[,2:9])
weights <- aicw(na.exclude(vals))
weights$model <- rownames(weights)
df <- rbind(df,weights)
}
df <- df[2:nrow(df),]
df$w <- as.numeric(df$w)
df$model <- factor(df$model, levels = c("HL", "HL_Temp" ,  "HL_Temp_Mammal" , "HL_Temp_Amphib",  "HL_Temp_Reptile", "HL_Temp_soft","HL_Temp_dPC1","HL_Temp_dPC2"))
ggboxplot(df, x = "model", y="w", add = "jitter") + rotate_x_text(45)

```

### Run 'scaled variables' model using top model: MFL ~ HL + Temp + logitMammal

```{r}
#Get pruned tree
tree <- read.newick(file="~/Desktop/FangData/Vipertree_named.tre")
check<-name.check(tree, data)
check
subtree <- drop.tip(tree, check$tree_not_data)
ult_tree<-force.ultrametric(subtree,method="extend") #tree only not ultrametric due to rounding error
ult_tree$edge.length <- ult_tree$edge.length
#########################################

data$stdHead <- scale(data$HL)
data$stdTemp <- scale(data$Temp)
data$stdMamm <- scale(data$logitMamm)



#run pgls with head size and environment and extract residuals for visualization
corr <-corPagel(1, phy=ult_tree, fixed = FALSE, form = ~Species)

pgls_mammal <- gls(MFL ~ stdHead + stdTemp , correlation = corr,
                   data = data, method = "ML")
summary(pgls_mammal)

resid <- pgls_mammal$residuals
resid <- cbind(resid, data$DietCat, data$Specialist)
resid <- as.data.frame(resid)
colnames(resid) <- c("Residuals","Diet","DietSpec")
resid$Residuals <- as.numeric(resid$Residuals)
resid$DietSpec <- factor(resid$DietSpec, levels = c("Mspec","Rspec","Aspec","Bspec","Ispec","Gen"))
resid$Diet <- factor(resid$Diet, levels = c("M","R","A","B","I"))
ggviolin(resid, x = "DietSpec", y = "Residuals", fill = "Diet",
         palette = c("#638457", "#90E39A", "#DDF093", "#F6D0B1","#CE4760"),
          add = c("boxplot"), add.params = list(fill = "Diet"))

```


### Get Distribution of standardized effect sizes
Will use posterior distribution of trees from Alencar et al. 2016 for bootstrapped parameter estimates and standardized effect size comparison

```{r, eval=TRUE, message=FALSE, results='hide', warning=FALSE}
treepath <- "~/Desktop/FangPaper_Clean/Posterior_trees_vipers/Test/"
treelist <- dir(treepath, pattern ="tree")

out.file <- ""
HL.out <- ""
temp.out <- ""
diet.out <- ""
tree.out <- ""
lambda.out <- ""
resid.out <- ""
for(i in 1:length(treelist)) {
  tryCatch({
    tree <- read.newick(file = paste(c(treepath,treelist[i]), collapse = ""))
    check<-name.check(tree, data)
    subtree <- drop.tip(tree, check$tree_not_data)
    ult_tree<-force.ultrametric(subtree,method="extend")
    corr<-corPagel(1, phy=ult_tree, fixed = FALSE, form = ~Species)
    #Edit model each time loop is run
    pglsModel <- gls(MFL ~ stdHead + stdTemp + stdMamm , correlation = corr,
                     data = data, method = "ML")
    cap <- capture.output(summary(pglsModel))
    lambda <- cap[11]
    diet <- cap[18]
    temp <- cap[17]
    HL <- cap[16]
    AIC <- MuMIn::AICc(pglsModel)
    resid <- as.data.frame(pglsModel$residuals)
    temp.out <- rbind(temp.out, temp)
    HL.out <- rbind(HL.out, HL)
    diet.out <- rbind(diet.out, diet)
    lambda.out <- rbind(lambda.out, lambda)
    AIC.out <- rbind(AIC.out, AIC)
    tree.out <- rbind(tree.out, treelist[i])
    resid.out <- cbind(resid.out, resid)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
}
```

### Extract standardized parameter estimates from above loop results
```{r}
out.file <- cbind(tree.out, HL.out, temp.out, diet.out)
cols <- c("Tree", "HL", "Temp", "Diet")
colnames(out.file) <- cols
out.file <- out.file[2:nrow(out.file), ]

HL <- cbind(tree.out, HL.out)
cols <- c("Tree","HL")
colnames(HL) <- cols
HL <- as.data.frame(HL)
HL <- HL[2:length(HL$HL),]
HL <- separate(HL, col = "HL", into = c("A","B","C","D","E","Estimate","SE","F","Tval","G","Pval","H"), sep = "\\s")
HL$Estimate <- as.numeric(HL$Estimate)

Temp <- cbind(tree.out, temp.out)
cols <- c("Tree","Temp")
colnames(Temp) <- cols
Temp <- as.data.frame(Temp)
Temp <- Temp[2:length(Temp$Temp),]
Temp <- separate(Temp, col = "Temp", into = c("A","B","C","D","E","Estimate","SE","F","Tval","G","Pval","H"), sep = "\\s")
Temp$Estimate <- as.numeric(Temp$Estimate)

Percent_Mammals <- cbind(tree.out, diet.out)
cols <- c("Tree","Percent_Mammals")
colnames(Percent_Mammals) <- cols
Percent_Mammals <- as.data.frame(Percent_Mammals)
Percent_Mammals <- Percent_Mammals[2:length(Percent_Mammals$Percent_Mammals),]
Percent_Mammals <- separate(Percent_Mammals, col = "Percent_Mammals", into = c("A","B","C","D","E","Estimate","SE","F","Tval","G","Pval","H"), sep = "\\s")
Percent_Mammals$Estimate <- as.numeric(Percent_Mammals$Estimate)

effects <- cbind(HL$Estimate,Temp$Estimate, Percent_Mammals$Estimate)
colnames(effects) <- c("HL","Temp","Percent.Mammals")
melted_effects <- melt(effects)
colnames(melted_effects) <- c("Row","Variable","Coefficient")
melted_effects$Variable <- factor(melted_effects$Variable, levels = c("Percent.Mammals","Temp","HL"))
p1 <- ggboxplot(melted_effects, x = "Variable", y = "Coefficient", fill = "Variable",
         palette = c("#638457", "#90E39A", "#DDF093"),
         add = "jitter")
p1 + coord_flip()

mean(effects[,1]) #average standardized head length effect
mean(effects[,2]) #average standardized Temperature effect
mean(effects[,3]) #average standardized Percent.Mammals effect


mean <- mean(effects[,1])
sd <- sd(data$HL)
min <- min(effects[,1]) 
max <- max(effects[,1])
print(paste("Each",sd, "mm increase in head length is associated with between",min,"and",max,
            "mm increase in fang length"))


mean <- mean(effects[,2])
sd <- sd(data$Temp)
min <- min(effects[,2]) 
max <- max(effects[,2])
range <- max(data$Temp) - min(data$Temp)
sd_range <- range/sd
print(paste("Each",sd, "degrees C increase in mean annual temperature is associated with between",min,"and",max,
            "mm increase in fang length"))
print(paste("Given an average fang length of",round(mean(data$MFL), digits = 2), "mm, we expect species in the coolest areas to have fangs that are smaller by about", 
            round((min/mean(data$MFL))*sd_range*100, digits=1),"to",round((max/mean(data$MFL))*sd_range*100, digits = 1),
            "percent compared to the warmest areas"))


mean <- mean(effects[,3])
sd <- sd(data$logitMamm)
min <- min(effects[,3]) 
max <- max(effects[,3])
print(paste("A",sd, "change in the logit(percent mammals in diet) is associated with between",min,"and",max,
            "mm increase in fang length"))


#back-transform logits
middle_low <- (0.975*exp(-sd/2)-0.025) / (1+exp(-sd/2))
middle_high <- (0.975*exp(sd/2)-0.025) / (1+exp(sd/2))
percent_mamm_change <- middle_high - middle_low
print(paste("Shifting from a",middle_low,"to a",middle_high, "% mammalian diet is associated with between",min,"and",max,
            "mm increase in fang length"))

upper_low <- (0.975*exp(sd)-0.025) / (1+exp(sd))
upper_high <- (0.975*exp(sd*2)-0.025) / (1+exp(sd*2))
print(paste("Shifting from a",upper_low,"to a",upper_high, "% mammalian diet is associated with between",min,"and",max,
            "mm increase in fang length, thus changing an average fang by", 
            (min/mean(data$MFL))*100,"to",(max/mean(data$MFL))*100,"percent"))

max_logit_diet_diff <- abs(logit(0) - logit(1))
delta_logit_stdevs <- max_logit_diet_diff/sd
print(paste("Given an average fang length of",round(mean(data$MFL), digits = 2), "mm, we expect fangs in snakes that do not eat mammals to be smaller by between", 
            round((min/mean(data$MFL))*delta_logit_stdevs*100,digits = 1),"and",round((max/mean(data$MFL))*delta_logit_stdevs*100,digits=1),
            "percent compared to snakes that feed solely on mammals"))
```


###END OF DOCUMENT###




