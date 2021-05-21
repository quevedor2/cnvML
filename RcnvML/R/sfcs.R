# Sweep space-filling curve
#  13  14  15  16
#   9  10  11  12
#   5   6   7   8
#   1   2   3   4
.sfc_sweep <- function(maxn){
  uids <- apply(expand.grid(c(0:maxn), 
                            c(0:maxn)), 
                1, paste, collapse="_")
  return(uids)
}

# Randomization of space-filling curve
#    7    8   13    5
#    6    2   12    1
#    4    3   15    9
#   16   11   14   10
.sfc_random <- function(maxn, seed=1234){
  set.seed(seed)
  grd <- expand.grid(c(0:maxn), 
                     c(0:maxn))
  grd <- grd[sample(c(1:nrow(grd)), size = nrow(grd), replace = FALSE),] # Randomization of UIDs
  
  uids <- apply(grd, 
                1, paste, collapse="_")
  return(uids)
}

# Scan space-filling curve
#   16   15   14   13
#    9   10   11   12
#    8    7    6    5
#    1    2    3    4
.sfc_scan <- function(maxn){
  grd <- expand.grid(c(0:maxn), 
                     c(0:maxn))
  
  # Reverse every other row of a sweep matrix
  for(i in c(0:maxn)[c(FALSE,TRUE)]){
    idx <- which(grd[,2] == i)
    grd[idx,1] <- rev(grd[idx,1])
  }
  
  uids <- apply(grd, 
                1, paste, collapse="_")
  return(uids)
}

# Diagonal space-filling curve
#    7   13   14   16
#    6    8   12   15
#    2    5    9   11
#    1    3    4   10
.sfc_diagonal <- function(maxn){
  minn <- 0
  
  compute_xy <- function(maxx, maxn, minx, dir='horizontal'){
    x <- minx
    while(minx < (maxn-1)){
      # left transition
      if(maxx < (maxn-1)) {
        maxx <- maxx+1
      } else {
        minx <- minx + 1
      }
      
      # right phase
      x <- switch(dir,
                  horizontal=c(x, c(minx:maxx)),
                  vertical=c(x, c(maxx:minx)))
      
      # right transition
      if(maxx < (maxn-1)){
        maxx <- maxx+1
      } else {
        minx <- minx + 1
      }
      
      # left phase
      x <- switch(dir,
                  horizontal=c(x, c(maxx:minx)),
                  vertical=c(x, c(minx:maxx)))
    }
    return(x)
  }
  x <- compute_xy(maxx=0, maxn=maxn, minx=minn, dir='horizontal')
  y <- compute_xy(maxx=0, maxn=maxn, minx=minn, dir='vertical')
  
  grd <- data.frame("x"=x, "y"=y)
  uids <- apply(grd, 
                1, paste, collapse="_")
  return(uids)
}

# Morton/Z-order space filling curve
#   11   12   15   16
#    9   10   13   14
#    3    4    7    8
#    1    2    5    6
#' @import morton
.sfc_morton <- function(maxn){
  grd <- fromMorton(c(0:(maxn^2-1)))
  grd <- do.call(cbind, grd)
  uids <- apply(grd, 
                1, paste, collapse="_")
  return(uids)
}


# Visualize the space-filling curves in text format
# Used to generate the comments for the sfc functions
.viz_mat <- function(uids){
  grd <- do.call(rbind, strsplit(uids, "_"))
  storage.mode(grd) <- 'integer'
  grd <- as.data.frame(grd)
  grd[,3] <- c(1:nrow(grd))
  n <- length(unique(sort(grd[,2])))
  
  demomat <- matrix(grd[order(grd[,'V2'], grd[,'V1']), 3], ncol=n)
  print(apply(t(demomat), 2, rev))
}


# n <- 4
# .viz_mat(.sfc_scan(n))
# .viz_mat(.sfc_random(n))
# .viz_mat(.sfc_sweep(n))
# .viz_mat(.sfc_morton(n))
# .viz_mat(.sfc_diagonal(n))
# 
# lv =4 
# pos <- t(sapply(0:(4^lv - 1), HilbertVis::hilbertCurvePoint, lv))
# colnames(pos) <- c('x', 'y')
# hc@POS = data.frame(x1 = pos$x[-n], y1 = pos$y[-n], x2 = pos$x[-1], y2 = pos$y[-1])
# 
