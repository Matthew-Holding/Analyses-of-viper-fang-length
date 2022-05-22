# Analyses used to test biological correlates of fang length in Viperidae

### Load necessary R packages

```{r}
#Statistics
library(nlme)
library(ape)
library(phytools)
library(geiger)
library(rptR)
library(picante)
library(car)
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
```{r}
#Read in species mean fang and other morphological data for each species with species and rownames for tree pruning
data <- read.table(file = "FangDataFinal.txt", row.names=1,header = TRUE, sep = "\t")

#Reread data with species as column 1 for pgls analyses:
compdata<-read.table(file = "FangDataFinal.txt", sep = "\t", header = TRUE)
#log transform head and fang lengths
compdata$logHL <- log(compdata$HL)
compdata$logMFL <- log(compdata$MFL)

#Read name-edited version of the ultrametric viperid species tree from Alencar et al. 2016
tree <- read.newick(file="Vipertree_named.tre")
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
phyrma <- phyl.RMA(x = x, y=y, tree = ult_tree, method = "lambda", lambda = .84)
phyrma

#Run pgls() function implemented in nlme package to get a pgls regression slope
#First, fit Pagel's lambda correlation structure with pruned species tree
pagel <-corPagel(1, phy=ult_tree, fixed = FALSE, form =~Species)

#run pgls and summarize to extract model estimates
pglsModelpa <- gls(logMFL ~ logHL , correlation = pagel,
                   data = compdata, method = "ML")
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

```{r}
#Read in dataset of individual, adult snake head and fang lengths
data.ind <- read.table("Adults_allometry.tsv", sep = "\t", header = TRUE)
data.ind$logMFL <- log(data.ind$MFL)
data.ind$logHL <- log(data.ind$HLmm)
```

### Fit models of head length vs. fang length for each species with 10 or more individuals measured (Fig 1B)

```{r, eval=TRUE}
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
Bpgls <- 1.35
BRMA <- 1.47
abline(v=Bpgls,col="blue", lty=3, lwd=2)
abline(v=BRMA,col="red", lty=3, lwd=2)
```

### Species' repeatability of relative fang length metric

```{r}

data.ind <- read.table("~/Desktop/FangLengthEvolution/SizeControl/Adults_allometry.tsv", sep = "\t", header = TRUE)
rpt(FL_HL ~ (1|CORRECT_scientificName), grname = c("CORRECT_scientificName"), data = data.ind, datatype = c("Gaussian"),
              CI = 0.95, nboot = 10000, npermut = 1000, parallel = TRUE,
              ncores = 7, ratio = TRUE, adjusted = TRUE, expect = "meanobs",
              rptObj = NULL, update = FALSE)
```

### Calculate phylogenetic signal in relative fang length evolution across Viperidae and plot in phenogram

```{r}

#make a named vector of relative fang length data:
RelFangL <- compdata[,5] #ratio
names(RelFangL) <- compdata$Species
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


### Get localities from GBIF database for each species in the dataset and extract and process associated environmental data

```{r}
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
write.table(df, file = "~/Desktop/FangPaper_Clean/LatRecs/All_Rec2.txt", quote = FALSE, sep = "\t")


library(raster)
library(sp)
library(elevatr)
library(rgdal)

#Pull down 2.5 resolution WorldClim environmental data
r <- getData("worldclim",var="bio",res=2.5)
r <- r[[c(1:19)]]
names(r) <- c(1:19)

#Associate species' localities with environmental data rasters
data <- read.table(file = "All_Rec2.txt", header = TRUE, sep = "\t", row.names = NULL)
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
write.table(loadings(pca,cutoff=0),"env_pca_loadings.txt",
            quote=FALSE,sep = "\t")

#aggregate data by species for downstream analysis, extracting means and medians for each data type
dfp <- cbind.data.frame(pcadf,pc1, pc2, pc3)
all_med <- aggregate(.~ Species, data = dfp, median)
all_mean <- aggregate(.~ Species, data = dfp, mean)
write.table(all_med,"envir_values_median.txt",row.names = FALSE,quote=FALSE,sep = "\t")
write.table(all_mean,"envir_values_mean.txt",row.names = FALSE,quote=FALSE,sep = "\t")
```

### Fit models predicting mean fang length (MFL) of 199 species full dataset

```{r}
data <- read.table(file = "~/Desktop/FangPaper_Clean/Allometry/FangDataFinal.txt", row.names=1,header = TRUE, sep = "\t")
data$logMFL <- log(data$MFL)
data$logHL <- log(data$HL)

#Remove rows with NAs
data <- data[complete.cases(data), ]
data$Species <- rownames(data)

treepath <- "~/Desktop/FangPaper_Clean/Posterior_trees_vipers/Test/"
treelist <- dir(treepath, pattern ="tree")
```

#### Note that below this point code is written to be run on a single 'pglsModel' at a time while looping through a posterior set of 200 trees from Alencar et al. 2016. This portion of code *will not run* with a single click. Instead, specific pglsModel anew each time.

```
#model list:
#pglsModel = gls(MFL ~  HL, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Prec, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + PC1, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp + Prec, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp + PC1, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp + Prec + Temp*Prec, correlation = corr, data = data, method = "ML")
#pglsModel = gls(MFL ~  HL + Temp + PC1 + Temp*PC1, correlation = corr, data = data, method = "ML")

out.file <- ""
HL.out <- ""
temp.out <- ""
PrecOrPC.out <- ""
int.out <- ""
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
    pglsModel = gls(MFL ~  HL + Temp + PC1 + Temp*PC1, correlation = corr, data = data, method = "ML")
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
    int.out <- rbind(int.out, int)
    lambda.out <- rbind(lambda.out, lambda)
    AIC.out <- rbind(AIC.out, AIC)
    tree.out <- rbind(tree.out, treelist[i])
    resid.out <- cbind(resid.out, resid)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
}

#Collect and collate model results using only one of the following 2-line snippets of code
#choose based on the variables in the model
out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, PrecOrPC.out, int.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "Prec", "Interaction", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, PrecOrPC.out, int.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "PC1", "Interaction", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, PrecOrPC.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "Prec", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, PrecOrPC.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "PC1", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, PrecOrPC.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "Prec", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, PrecOrPC.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Prec", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Temp", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "Prec", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, temp.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "PC1", "AIC")

out.file <- cbind(tree.out, lambda.out, HL.out, AIC.out)
cols <- c("Tree", "lambda", "HL", "AIC")

colnames(out.file) <- cols
out.file <- out.file[2:nrow(out.file), ]

#change outfile name as appropriate:
write.table(out.file, file = "~/Desktop/NewAIC/HL_TempPC1_Int.txt", row.names = TRUE, 
            col.names = TRUE, quote = FALSE, sep = "\t")

AICs <- cbind(tree.out, AIC.out)
cols <- c("Tree","AIC")
colnames(AICs) <- cols
AICs <- as.data.frame(AICs)
AICs <- AICs[2:length(AICs$AIC),]
AICs <- separate(AICs, col = "AIC", into = c("BL","BL2","AIC"), sep = "\\s")
AICs <- AICs[,c(1,4)]
write.table(AICs, file = "~/Desktop/NewAIC/AIC_HL_TempPC1_Int.txt", row.names = FALSE, 
            col.names = TRUE, quote = FALSE, sep = "\t")

#write out residuals from only the top model when it is ran
write.table(resid.out, file = "~/Desktop/NewAIC/ResidualsOut.txt", row.names = TRUE, 
            col.names = TRUE, quote = FALSE, sep = "\t")

#Reread and then merge all data AIC frames 
df1 <- read.table("AIC_HLonly.txt",header = TRUE)
colnames(df1) <- c("Tree", "HL")
df2 <- read.table("AIC_HL_Temp.txt",header = TRUE)
colnames(df2) <- c("Tree", "HL_Temp")
df3 <- read.table("AIC_HL_Prec.txt",header = TRUE)
colnames(df3) <- c("Tree", "HL_Prec")
df4 <- read.table("AIC_HL_PC1.txt",header = TRUE)
colnames(df4) <- c("Tree", "HL_PC1")
df5 <- read.table("AIC_HL_TempPrecip.txt",header = TRUE)
colnames(df5) <- c("Tree", "HL_TempPrec")
df6 <- read.table("AIC_HL_TempPC1.txt",header = TRUE)
colnames(df6) <- c("Tree", "HL_TempPC1")
df7 <- read.table("AIC_HL_TempPrec_Int.txt",header = TRUE)
colnames(df7) <- c("Tree", "HL_TempPrec_Int")
df8 <- read.table("AIC_HL_TempPC1_Int.txt",header = TRUE)
colnames(df8) <- c("Tree", "HL_TempPC1_Int")

df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8)

AIC_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  
AIC_merged <- AIC_merged[complete.cases(AIC_merged), ]
write.table(AIC_merged, file = "~/Desktop/NewAIC/AIC_merged.txt", row.names = FALSE, 
            col.names = TRUE, quote = FALSE, sep = "\t")
#read in Akaike weights calculated in MS Excel:
weights <- read.table("weights.txt", header = TRUE)
df_melt <- melt(weights)
colnames(df_melt) <- c("Model", "Model Weight")
ggboxplot(df_melt, x = "Model", y="Model Weight", add = "jitter") + rotate_x_text(45)


#Make plot showing the results of the top model: MFL ~ HL + Temp*envPC1
library(RColorBrewer)
cols = brewer.pal(3, "RdPu")
# Define colour pallete
# Use the following line with RColorBrewer
pal = colorRampPalette(cols)

# Rank variable for colour assignment
data$order = findInterval(data$Temp, sort(data$Temp))
#Facilitate proper point sizing according to environmental PC1 scores
PC1_cex <- (data$PC1) + 10000
# Make plot using pgls model of best tree from Alencar et al. 2016 for slope as defined above
pglsModelpa <- gls(logMFL ~ logHL , correlation = pa,
                   data = compdata, method = "ML")
plot(MFL ~ HL, data=data, pch=19, col=pal(nrow(data))[data$order], cex=(PC1_cex)/10000)
#add pgls best fit line for HL vs MFL
abline(pglsModelpa, col = "black")


#Extract phylogenetic residuals from head length vs fang length to explore explanatory power of temperature
pglsModel_reduced <- gls(MFL ~ HL , correlation = pa,
                   data = data, method = "ML")

resid <- as.data.frame(pglsModel_red$residuals)
PC1 <- data$PC1
Temp <- data$Temp
resid_df <- cbind(resid,PC1,Temp)
colnames(resid_df)<-c("Residuals","envPC1","Temperature")
resid_df$Species <- rownames(resid_df)

ols <- lm(Residuals ~ I(Temperature/10), data=resid_df)
plot(Residuals~I(Temperature/10), data=resid_df, xlab = "Temp (C)",ylab = "Residual Fang Length")
abline(ols, col = "black")
summary(ols)

```


