#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(ggplot2)
library(gtable)
library(gridExtra)
library(grid)
library(dplyr)
library(maps)
library(rgeos)
library(sp)
library(maptools)
library(plyr)
library(lubridate)
library(shinyjs)
#######
#Set working directory
dir <- "C:/Users/Kyle/Desktop/Geog418FINAL/Kyle Harnadek/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR(dsn = "Pm25Sample.shp", layer = "Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))

#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 
#Read in dissemination tract shapefile:
census.tracts <- readOGR(dsn = "BC_DA.shp", layer = "BC_DA") 
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income \n In MVRD",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Income

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=10000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)

#ADD CODE SECTION:

####################
#Descriptive Statistics

meanIncome <- mean(income.tracts$Income, na.rm = TRUE)
meanPol <- mean(pm2.5$PM25, na.rm = TRUE)

sdIncome <- sd(income.tracts$Income, na.rm = TRUE)
sdPol <- sd(pm2.5$PM25, na.rm = TRUE)

medIncome <- median(income.tracts$Income, na.rm = TRUE)
medPol <- median(pm2.5$PM25, na.rm = TRUE)

modeIncome <- as.numeric(names(sort(table(income.tracts$Income), decreasing = TRUE))[1])
modePol <- as.numeric(names(sort(table(pm2.5$PM25), decreasing = TRUE))[1])

skewIncome <- skewness(income.tracts$Income, na.rm = TRUE)
skewPol <- skewness(pm2.5$PM25, na.rm = TRUE)

kurtIncome <- kurtosis(income.tracts$Income, na.rm = TRUE)
kurtPol <- kurtosis(pm2.5$PM25, na.rm = TRUE)

CoVIncome <- (sdIncome / meanIncome) *100
CoVPol <- (sdPol / meanPol) *100

normIncome <- shapiro.test(income.tracts$Income)$p.value
normPol <- shapiro.test(pm2.5$PM25)$p.value

#Create a table of descriptive stats

samples = c("Income", "Pollution") #Create an object for the labels
means = c(meanIncome, meanPol) #Create an object for the means
sd = c(sdIncome, sdPol) #Create an object for the standard deviations
median = c(medIncome, medPol) #Create an object for the medians
mode <- c(modeIncome, modePol) #Create an object for the modes
skewness <- c(skewIncome, skewPol) #Create an object for the skewness
kurtosis <- c(kurtIncome, kurtPol) #Create an object for the kurtosis
CoV <- c(CoVIncome, CoVPol) #Create an object for the CoV
normality <- c(normIncome, normPol) #Create an object for the normality PVALUE

##Check table values for sigfigs
normality <- signif(normality, digits = 2)
means <- round(means, 2)
sd <- round(sd, 2)
median <- round(median, 2)
mode <- round(mode, 2)
skewness <- round(skewness, 2)
kurtosis <- round(kurtosis, 2)

data.for.table = data.frame(samples, means, sd, median, mode, skewness, kurtosis, CoV, normality)

table1 <- tableGrob(data.for.table, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Income and Pollution Descriptive Stats", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)


grid.arrange(table1, newpage = TRUE)

#######################
#Global Moran's I

income.nb <- poly2nb(income.tracts) #Create neighbour object
income.net <- nb2lines(income.nb, coords=coordinates(income.tracts)) #create lines object with all previously defined Queen neighbour connnections
crs(income.net) <- crs(income.tracts) #Apply same coordinate system from income tract object


#Mapping Queens Case
tm_shape(income.tracts) + tm_borders(col='#DDEAF6') + 
  tm_shape(income.net) + tm_lines(col='red') +
  tm_layout(panel.labels = "Queen Weight")
#Create weight list
income.lw <- nb2listw(income.nb, zero.policy = TRUE, style = "W")
print.listw(income.lw, zero.policy = TRUE) #print weight list with previously defined parameters

#Morans test is performed with income tract object, selected variable, and weights matrix
mi <- moran.test(income.tracts$Income, income.lw, zero.policy = TRUE)
mi


moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(income.lw)


mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- (mI - eI) / sqrt(var)

if ( z > 1.96) {
  print("Positive correlation, reject null hypothesis")
} else if ( z < 1.96) {
  print("Negative correlation, reject null hypothesis")
} else {
  print("No significant spatial autocorrelation")
}

#####################


##################
#IDW

#IDW Interpolation
P.idw <- gstat::idw(PM25 ~ 1, pm2.5, newdata=grd, idp=4)
r       <- raster(P.idw)
r.m     <- mask(r, income.tracts)

idwmap <- tm_shape(r.m) + 
  tm_raster(n=10,palette = "seq", 
            title="IDW \nPredicted PM2.5 \n(in ppm)") + 
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE) +
  tm_layout(aes.palette = list(seq = "-RdBu"))

IDW.out <- vector(length = length(pm2.5))
for (i in 1:length(pm2.5)) {
  IDW.out[i] <- idw(PM25 ~ 1, pm2.5[-i,], pm2.5[i,], idp=4)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ pm2.5$PM25, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ pm2.5$PM25), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt( sum((IDW.out - pm2.5$PM25)^2) / length(pm2.5))




######################

pm2.5$x <- coordinates(pm2.5)[,1]
pm2.5$y <- coordinates(pm2.5)[,2]
#Point Pattern Analysis
zd <- zerodist(pm2.5)
zd

#if there are duplicates, remove them
pm2.5 <- remove.duplicates(pm2.5)

#create an "extent" object which can be used to create the observation window for spatstat
pm2.5.ext <- as.matrix(extent(pm2.5)) 

#observation window
window <- as.owin(list(xrange = pm2.5.ext[1,], yrange = pm2.5.ext[2,]))

#create ppp oject from spatstat
pm2.5.ppp <- ppp(x = pm2.5$x, y = pm2.5$y, window = window)
plot(pm2.5.ppp, par(oma=c(0,0,0,0)),
     main = "Station Locations")

#QUADRAT ANALYSIS
quads <- 20

qcount <- quadratcount(pm2.5.ppp, nx = quads, ny = quads)

plot(pm2.5.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red")

qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df <- plyr::count(qcount.df,'Freq')

##Change the column names so that x=number of points and f=frequency of quadrats with x point.
colnames(qcount.df) <- c("x","f")

sum.f.x2
sum.f.x2 <- sum(qcount.df$x * qcount.df$f^2)

M <- sum(qcount.df$f)

N <- sum(qcount.df$x * qcount.df$f)

sum.fx.2 <- sum((qcount.df$x * qcount.df$f)^2)


VAR <- (sum.f.x2 - (sum.fx.2 / M)) / (M - 1)

MEAN <- N/M

VMR <- VAR/MEAN



##Finally, perform the test statistic to test for the existence of a random spatial pattern.
chi.square = VMR*(M-1)
p = 1 - pchisq(chi.square, (M - 1))
p  

##Nearest Neighbour Distance
###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(pm2.5.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = mean(nndist(pm2.5.ppp, k=1))

#mean nearest neighbour for random spatial distribution

studyArea <- sum(area(income.tracts))
studyAreakm2 <- studyArea / 1000000
studyArea
pointDensity <- studyAreakm2 / sum(pm2.5.ppp$n)

r.nnd = 1/(2*sqrt(pointDensity))

d.nnd = 1.07453/(sqrt(pointDensity))

R = nnd/r.nnd

SE.NND <- 0.26136/(sqrt((sum(pm2.5.ppp$n))*pointDensity))

z = (nnd-r.nnd)/SE.NND
z

