# Analyses used to test biological correlates of fang length in Viperidae


library(nlme)

#Read in dataset of individual, adult snake head and fang lengths
```{r}
data <- read.table("Adults_allometry.tsv", sep = "\t", header = TRUE)
data$logMFL <- log(data$MFL)
data$logHL <- log(data$HLmm)
```

#fit relationships between log fang length and log head length for each species
fits <- nlme::lmList(logMFL ~ logHL | CORRECT_scientificName, data=data, )
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
Bpgls <- 1.36
BRMA <- 1.47
abline(v=Bpgls,col="blue", lty=3, lwd=2)
abline(v=BRMA,col="red", lty=3, lwd=2)
