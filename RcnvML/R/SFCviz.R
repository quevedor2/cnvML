#' Visualizes the euclidean distance between 2D space-filling curves
#' 
#' @description Plots the euclidean distance matrix corresponding to 
#' the results from as.matrix(dist(data.frame(x1, y1))). 'frac_df' acts
#' as a supporting dataframe to label the chromosomes along the x and
#' y-axes.
#'
#' @param xdist distance matrix between a two-column dataframe
#' @param frac_df data frame to plot the chromosomes along the x and y-
#' axis with a corresponding colour
#' @importFrom assertthat assert_that
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap decorate_heatmap_body
#' @importFrom grid grid.segments
#' @importFrom grid grid.polygon
#' @export
#' 
plotEuclidDist <- function(xdist, frac_df){
  assert_that(is.matrix(xdist), ncol(xdist)==ncol(xdist),
              "'xdist' is not a properly formed square distance matrix")
  assert_that(is.data.frame(frac_df), ncol(frac_df)==3,
              all(c('start', 'end', 'col') %in% colnames(frac_df)),
              "'frac_df' must be a 3-column data frame of 'start', 'end', 'col'")
  
  ht = Heatmap(xdist, name = "dist", 
               cluster_rows = F, cluster_columns = F, 
               show_row_names = F, show_column_names = F, 
               heatmap_legend_param = list(title = "distance"))
  draw(ht, padding = unit(c(5, 5, 5, 2), "mm"))
  
  line_seg <- frac_df$start[-1]
  decorate_heatmap_body("dist", {
    grid.segments(x0=c((1-line_seg), rep(0, length(line_seg))), 
                  x1=c((1-line_seg), rep(1, length(line_seg))), 
                  y0=c(rep(0, length(line_seg)), line_seg),
                  y1=c(rep(1, length(line_seg)), line_seg), 
                  gp = gpar(lty = 2))
    apply(frac_df, 1, function(i, xpos, ypos){
      # side along the y-axis
      i1 <- as.numeric(i[1])
      i2 <- as.numeric(i[2])
      grid.polygon(x=c(xpos, xpos, xpos-0.03), y=c(i2, i1, i2+(i1-i2)/2),
                   gp=gpar(fill=i[3]))
      
      # top along the x-axis
      i1 <- 1-i1
      i2 <- 1-i2
      grid.polygon(x=c(i2, i1, i1+(i2-i1)/2), y=c(ypos, ypos, ypos+0.03), 
                   gp=gpar(fill=i[3]))
    }, xpos=0, ypos=1)
  })
}

