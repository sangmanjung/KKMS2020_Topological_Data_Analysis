library(TDA)
library(readxl)

## Data load
path_name <- 'insert your file path'
sheet_name <- excel_sheets(path = path_name)

tmp <- read_excel(path = path_name, sheet = sheet_name[1])
r <- nrow(tmp) # the number of observations
c <- ncol(tmp) # the number of variables
p <- length(sheet_name) # the number of subjects

X <- array(0,dim = c(r,c,p)) # data array
cx <- array(0,dim = c(c,c,p)) # the correlation-based distance

## Matrix representations
for (page in 1:p){
  tmp <- read_excel(path = path_name, sheet = sheet_name[page])
  X[,,page] <- data.matrix(tmp) # dataset
  cx[,,page] <- sqrt(1-abs(cor(X[,,page]))) # point clouds
}

## Single linkage dendrogram
MST <- array(0,dim = c(c-1,3,p)) # minimum spanning tree
par(mfrow = c(3,3))
for (i in 1:p){
  dd <- hclust(d = as.dist(cx[,,i]), method = 'single')
  MST[,,i] <- cbind(as.matrix(dd[['merge']]),as.matrix(dd[['height']]))
  plot(as.dendrogram(dd),main = sprintf('subject %s',i),
       xlab = 'the number of nodes',ylab = 'filtration')
}

## Persistent homology
maxdimension <- 2 # connected components, loops, and voids (b0, b1, and b2)
maxscale <- 2 # limit of the filtration
Rips_all <- array(0,dim = c(1,p)) # Vietoris-Rips filtration
for (i in 1:p){
  Rips_all[i] <- ripsDiag(X = cx[,,i],
                          maxdimension,maxscale,library = 'GUDHI',
                          dist = 'euclidean',printProgress = FALSE)
}

## Persistence diagram & Barcode
par(mfrow = c(3,3))
for (i in 1:p){
  plot(Rips_all[[i]], main = sprintf('subject %s',i)) # PDs
}
for (i in 1:p){
  plot(Rips_all[[i]], barcode = TRUE, main = sprintf('subject %s',i)) # barcodes
}

## Bottleneck distance
BT <- array(0,dim = c(p,p,maxdimension+1)) # pre-allocation
for (k in 1:(maxdimension+1)){
  for (i in 1:p){
    for (j in 1:p){
      BT[i,j,k] <- bottleneck(Rips_all[[i]],Rips_all[[j]],dimension = k-1)
    }
  } 
}

## Multidimensional scaling
mdsdim <- 2
MDS_all <- array(0,dim = c(p,mdsdim,maxdimension+1)) # pre-allocation
for (i in 1:(maxdimension+1)){
  MDS_all[,,i] <- cmdscale(d = BT[,,i],k = mdsdim) # MDS coordinates
}

par(mfrow = c(1,3)) # MDS plot
for (i in 1:(maxdimension+1)){
  plot(MDS_all[,,i],type = 'n',xlim = c(-max(MDS_all[,1,i]),max(MDS_all[,1,i])),
       ylim = c(-max(MDS_all[,2,i]),max(MDS_all[,2,i])),
       xlab = 'Dim1',ylab = 'Dim2',main = sprintf('Bottleneck for Betti-%d',i-1))
  segments(-1,0,1,0,lty = 'dotted')
  segments(0,-1,0,1,lty = 'dotted')
  # number 1~9 == subject 1~9
  text(MDS_all[,,i],rownames(MDS_all[,,i]),cex = 1.5,col = 'Blue')
}
