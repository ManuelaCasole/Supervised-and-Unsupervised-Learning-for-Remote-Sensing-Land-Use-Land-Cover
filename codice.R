library(raster)
library(plotly)

# Blue
b2 <- raster('romacrop/B2.tif')
# Green
b3 <- raster('romacrop/B3.tif')
# Red
b4 <- raster('romacrop/B4.tif')
# Near Infrared (NIR)
b5 <- raster('romacrop/B5.tif')
b <- stack(b2,b3,b4,b5)

b6 <- raster('romacrop/B6.tif')


filenames <- paste0('romacrop/B', c(2,3,4,5,6,7), ".tif")
filenames
landsat <- stack(filenames)
landsat

#### image visualization 
par(mfrow = c(2,2))
plot(b2, main = "Blue", col = gray(0:100 / 100))
plot(b3, main = "Green", col = gray(0:100 / 100))
plot(b4, main = "Red", col = gray(0:100 / 100))
plot(b5, main = "NIR", col = gray(0:100 / 100))


####   RGB and   false color--- negative image
par(mfrow = c(1,2))
landsatRGB <- stack(b2, b3, b4)
plotRGB(landsatRGB, axes=TRUE, stretch="lin", main="Landsat True Color Composite")
landsatFCC <- stack(b5, b4, b3)
plotRGB(landsatFCC, axes=TRUE, stretch="lin", main="Landsat False Color Composite")

### renaming
names(landsat) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2')
names(landsat)



###normalized difference vegetation index
vi <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}
ndvi <- vi(landsat, 4, 3)
plot(ndvi, col = rev(terrain.colors(10)), main = "Landsat-ndvi")

######### histogram
# view histogram of data
{hist(ndvi,
      main = "Distribution of NDVI values",
      xlab = "NDVI",
      ylab= "Frequency",
      col = "wheat",
      xlim = c(-0.5, 1),
      breaks = 30,
      xaxt = 'n')
  axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))}


######## hide non green areas with v.i. <=0.4
veg <- reclassify(ndvi, cbind(-Inf, 0.4, NA))
plot(veg, main='Vegetation')

vegc <- reclassify(ndvi, c(-Inf, 0.25, 1, 0.25, 0.3, 2, 0.3, 0.4, 3, 0.4, 0.5, 4, 0.5, Inf, 5)) 
plot(vegc,col = rev(terrain.colors(4)), main = 'NDVI based thresholding')






######## PCA dimensionality reduction
set.seed(1)
sr <- sampleRandom(landsat, 10000)
pca <- prcomp(sr, scale = TRUE)
pca
summary(pca)
screeplot(pca)

pci <- predict(landsat, pca, index = 1:2)

plot(pci[[1]])

pc2 <- reclassify(pci[[2]], c(-Inf,0,1,0,Inf,NA))
par(mfrow = c(1,2))
plotRGB(landsatFCC, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", main =
          "Landsat False Color Composite")
plot(pc2, legend = FALSE)




########## UNSUPERVISED ###########

library(raster)
filenames <- paste0('romacrop/B', c(2,3,4,5,6,7), ".tif")
filenames
landsat8 <- stack(filenames)
landsat8
names(landsat8) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2')
ndvi <- (landsat8[['NIR']] - landsat8[['red']]) / (landsat8[['NIR']] + landsat8[['red']])
nr <- getValues(ndvi)

i <- !is.na(nr)
set.seed(99)
kmncluster <- kmeans(nr[i], centers = 4, iter.max = 200, nstart = 5, algorithm="Lloyd")
nr[i] <- kmncluster$cluster
knr <- setValues(ndvi, nr)


# Define a color vector for 4 clusters (learn more about setting the color later)
mycolor <- c("green","yellow","blue" ,"red")
par(mfrow = c(1,2))
plot(ndvi, col = rev(terrain.colors(10)), main = 'Landsat-NDVI')
plot(knr, main = 'Unsupervised k-means', col = mycolor )



############SUPERVISED##########




###################   data processing 
m <- 4
library(raster)
nlcd <- raster('romacrop/corine1.tif')
plot(nlcd)


km <- getValues(nlcd)
km[is.na(km)] <- 0
i <- !is.na(km)
set.seed(99)
kmncluster <- kmeans(km[i], centers = m, iter.max = 200, nstart = 5, algorithm="Lloyd")
km[i] <- kmncluster$cluster
nlcd <- setValues(nlcd, km)
plot(nlcd, col=c("blue","green","red" ,"yellow"))


nlcdclass <- c("Waterbodies", "Forest", "Artificial", "Agricultural")
classdf <- data.frame(classvalue1 = c(1:m), classnames1 = nlcdclass)

classcolor <- c("blue","green","red" ,"yellow")


# Now we ratify (RAT = "Raster Attribute Table") the ncld (define RasterLayer as a categorical variable). This is helpful for plotting.
nlcd2018 <- nlcd
nlcd2018 <- ratify(nlcd2018)
rat <- levels(nlcd2018)
levels(nlcd2018) <- rat

# Sampling
set.seed(99)
samp2018 <- sampleStratified(nlcd2018, size = 100, na.rm = TRUE, sp = TRUE)
table(samp2018$corine1)



##### recall spectral bands

filenames <- paste0('romacrop/B', c(2,3,4,5,6,7), ".tif")
filenames


landsat8 <- stack(filenames,ndvi)
landsat8

names(landsat8) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI')
names(landsat8)


# Extract the layer values for the locations
sampvals <- extract(landsat8, samp2018, df = TRUE)
# sampvals no longer has the spatial information. To keep the spatial information you use `sp=TRUE` argument in the `extract` function.
# drop the ID column
sampvals <- sampvals[, -1]
# combine the class information with extracted values
sampdata <- data.frame(classvalue = samp2018@data[,2], sampvals)






par(mfrow=c(1,1))
#################   classification tree
library(rpart)
# Train the model
cart <- rpart(as.factor(classvalue)~ blue + green + NIR + SWIR1 + SWIR2 + NDVI, data=sampdata, method = 'class', minsplit = 5)
# print(model.class)
# Plot the trained classification tree
plot(cart, uniform=TRUE, main="Classification Tree")
text(cart, cex = 0.8)


# Now predict the subset data based on the model; prediction for entire area takes longer time
pr2018 <- predict(landsat8, cart, type='class')
pr2018


plot(pr2018, col=c("blue","green","red" ,"yellow"),main = "Decision Tree classification of Landsat 8" )





library(dismo)
set.seed(99)
j <- kfold(sampdata, k = 5, by=sampdata$classvalue)
table(j)


x <- list()
for (k in 1:5) {
  train <- sampdata[j!= k, ]
  test <- sampdata[j == k, ]
  cart <- rpart(as.factor(classvalue)~., data=train, method = 'class', minsplit = 5)
  pclass <- predict(cart, test, type='class')
  # create a data.frame using the reference and prediction
  x[[k]] <- cbind(test$classvalue, as.integer(pclass))
}

y <- do.call(rbind, x)
y <- data.frame(y)
colnames(y) <- c('observed', 'predicted')


conmat <- table(y)

# change the name of the classes
colnames(conmat) <- classdf$classnames
rownames(conmat) <- classdf$classnames
conmat

# number of cases
n <- sum(conmat)
n

# number of correctly classified cases per class
diag <- diag(conmat)


# Overall Accuracy
OA <- sum(diag) / n
OA


# observed (true) cases per class
rowsums <- apply(conmat, 1, sum)
p <- rowsums / n

# predicted cases per class
colsums <- apply(conmat, 2, sum)
q <- colsums / n

expAccuracy <- sum(p*q)
kappa <- (OA - expAccuracy) / (1 - expAccuracy)
kappa

# Producer accuracy
PA <- diag / colsums


# User accuracy
UA <- diag / rowsums


outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
outAcc






# required packages.    random forest
library(party)
library(randomForest)
sampdata[,1] <- as.factor(sampdata[,1])
set.seed(12345)

# Create the forest.
output.forest <- randomForest(classvalue ~ blue + green + NIR + SWIR1 + SWIR2 + NDVI, 
                              data = na.omit(sampdata))


pr2018 <- predict(landsat8, output.forest, type='class')
pr2018
plot(pr2018, col=c("blue","green","red" ,"yellow"),main = "Random Forest classification of Landsat 8" )

# View the forest results.
print(output.forest) 

# Importance of each predictor.
print(importance(output.forest,type = 2))


conmat <- as.matrix(output.forest$confusion)

conmat<- conmat[,-5]
colnames(conmat) <- classdf$classnames
rownames(conmat) <- classdf$classnames
conmat <-as.table(conmat)



# number of cases
n <- sum(conmat)
n

# number of correctly classified cases per class
diag <- diag(conmat)


# Overall Accuracy
OA <- sum(diag) / n
OA


# observed (true) cases per class
rowsums <- apply(conmat, 1, sum)
p <- rowsums / n

# predicted cases per class
colsums <- apply(conmat, 2, sum)
q <- colsums / n

expAccuracy <- sum(p*q)
kappa <- (OA - expAccuracy) / (1 - expAccuracy)
kappa

# Producer accuracy
PA <- diag / colsums


# User accuracy
UA <- diag / rowsums


outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
outAcc