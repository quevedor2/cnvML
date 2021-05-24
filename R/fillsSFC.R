library(GenomicRanges)
library(RcnvML)

## Set up the initial reference ordered dataframe
outdir <- "~"
ord <-6
gbin_pos <- setupRefHcMatrix(order=ord)

## Details of the rectangles to overlap the input 2D matrix with
rcts <- list('row'=c('xleft'=1, 'xright'=64, 'ytop'=64, 'ybottom'=50),
             'square'=c('xleft'=30, 'xright'=50, 'ytop'=50, 'ybottom'=30),
             'column'=c('xleft'=20, 'xright'=30, 'ytop'=64, 'ybottom'=1))
rct_cols <- c('row'='maroon', 'square'='navyblue', 'column'='#8856a7')

## Cycle through each SFC and view the 1D -> 2D coverage
sfcs <- c('random', 'sweep', 'scan', 'diagonal', 'morton', 'hilbert')
lapply(names(rcts), function(id){
  rct <- rcts[[id]]
  rct_col <- rct_cols[id]
  
  pdf(file.path(outdir, paste0("sfc_", id, ".pdf")))
  # Plot the 2d space and where the box is located
  plot(0, type='n', xlim=c(0, max(gbin_pos$ord$x1)), ylim=c(0, max(gbin_pos$ord$y1)), 
       xaxt='n', yaxt='n', xlab='', ylab='')
  rect(xleft = rct['xleft']-1, ybottom = rct['ybottom']-1, 
       xright = rct['xright']-1, ytop = rct['ytop']-1, 
       col=rct_col, border = NA)
  
  par(mfrow=c(length(sfcs),1))
  # Cycle through each sfc and plot the 1D representation of that box
  na <- lapply(sfcs, function(sfc){
    if(sfc != 'hilbert'){
      bins <- mapSFC(sfc=sfc, hc_ord=gbin_pos$ord)
    } else {
      bins <- gbin_pos$ord
    }
    bins_spl <- splSquareFromBins(xleft = rct['xleft'], xright = rct['xright'], 
                                  ytop = rct['ytop'], ybottom = rct['ybottom'], bins)
    
    par(mar=c(2, 7, 1, 2.1))
    plot(0, type='n', xlim=c(0, max(bins$gord)+1), ylim=c(0,1), 
         yaxt='n', xaxt='n', xlab='', ylab='')
    axis(side = 2, at=0.5, labels=sfc, tick = FALSE, las=1)
    apply(bins_spl$`TRUE`, 1, function(i){
      rect(xleft = as.integer(i['gord'])-0.5, ybottom = 0, 
           xright = as.integer(i['gord'])+0.5, ytop = 1, 
           col=rct_col, border=NA)
    })
  })
  close.screen(all.screens=TRUE)
  dev.off()
})

