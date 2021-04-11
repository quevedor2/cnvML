#' A copy of HilbertCurve() that lets POS inputs
#'
#' @param s position that will be mapped as the start of the Hilbert curve. The value should a single numeric value.
#'   If it is a vector, the minimum is used.
#' @param e position that will be mapped as the end of the Hilbert curve. The value should a single numeric value.
#'   If it is a vector, the maximum is used.
#' @param level iteration level of the Hilbert curve. There will by ``4^level - 1`` segments in the curve.
#' @param mode "normal" mode is used for low ``level`` value and "pixel" mode is always used for high ``level`` value,
#'   so the "normal" mode is always for low-resolution visualization while "pixel" mode is used for high-resolution visualization.
#'   See 'details' for explanation.
#' @param reference whether add reference lines on the plot. Only works under 'normal' mode. The reference line
#'   is only used for illustrating how the curve folds.
#' @param reference_gp graphic settings for the reference lines. It should be specified by `grid::gpar`.
#' @param arrow whether add arrows on the reference line. Only works under 'normal' mode.
#' @param zoom Internally, position are stored as integer values. To better map the data to the Hilbert curve, 
#'   the original positions are zoomed
#'   according to the range and the level of Hilbert curve. E.g. if 
#'   the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4,
#'   the positions will be zoomed by ~x2000 so that values like 1.5, 1.555 can be mapped
#'   to the curve with more accuracy. You don't need to care the zooming thing, 
#'   proper zooming factor is calculated automatically.
#' @param newpage whether call `grid::grid.newpage` to draw on a new graphic device.
#' @param background_col background color.
#' @param background_border background border border.
#' @param title title of the plot.
#' @param title_gp graphic parameters for the title. It should be specified by `grid::gpar`.
#' @param start_from which corner on the plot should the curve starts?
#' @param first_seg the orientation of the first segment.
#' @param legend a `grid::grob` object, a `ComplexHeatmap::Legends-class` object, or a list of them.
#' @param padding padding around the Hilbert curve.
#' @param sfc_pos a dataframe object of x1, x2, y1, y2 positions according to POS 
#' 
#' @import grid
#' @import HilbertCurve
#' @import ComplexHeatmap
#' @import IRanges
#' @importFrom methods new
#' @importFrom HilbertVis hilbertCurve
#' @export
SfcCurve <- function(s, e, level = 4, mode = c("normal", "pixel"),
                       reference = FALSE, reference_gp = gpar(lty = 3, col = "#999999"), 
                       arrow = TRUE, zoom = NULL, newpage = TRUE, 
                       background_col = "transparent", background_border = NA, 
                       title = NULL, title_gp = gpar(fontsize = 16), 
                       start_from = c("bottomleft", "topleft", "bottomright", "topright"),
                       first_seg = c("horizontal", "vertical"), legend = list(), 
                       padding = unit(2, "mm"), sfc_pos=NULL){
  hc = new("HilbertCurve")
  level = as.integer(level)
  
  hc@data_range = c(min(s), max(e))
  
  if(s > e) {
    stop("`s` must be smaller than `e`.")
  }
  if(s < 0) {
    hc@OFFSET = s
  } else {
    hc@OFFSET = 0
  }
  s = HilbertCurve::hc_offset(hc, s)
  e = hc_offset(hc, e)
  
  if(is.null(zoom)) zoom = 10*(4^level-1)/(e - s)
  hc@ZOOM = zoom
  
  pos = HilbertCurve:::hilbertCurve(level)
  start_from = match.arg(start_from)[1]
  first_seg = match.arg(first_seg)[1]
  
  if(start_from == "bottomleft") {
    if(first_seg == "horizontal") {
      # default, no transformation
    } else if(first_seg == "vertical") {
      pos = HilbertCurve:::rotate(pos, degree = 90)
      pos = HilbertCurve:::flip(pos, direction = "horizontal")
    }
  } else if(start_from == "topleft") {
    if(first_seg == "horizontal") {
      pos = HilbertCurve:::rotate(pos, degree = 180)
      pos = HilbertCurve:::flip(pos, direction = "horizontal")
    } else {
      pos = HilbertCurve:::rotate(pos, degree = 270)
    }
  } else if(start_from == "bottomright") {
    if(first_seg == "horizontal") {
      pos = HilbertCurve:::flip(pos, direction = "horizontal")
    } else {
      pos = HilbertCurve:::rotate(pos, degree = 90)
    }
  } else if(start_from == "topright") {
    if(first_seg == "horizontal") {
      pos = HilbertCurve:::rotate(pos, degree = 180)
    } else {
      pos = HilbertCurve:::rotate(pos, degree = 270)
      pos = HilbertCurve:::flip(pos, direction = "horizontal")
    }
  }
  
  hc@start_from = start_from
  hc@first_seg = first_seg
  n = 4^level  # number of points
  
  breaks = round(seq(zoom(hc, s), zoom(hc, e), length = n)) 
  
  # n - 1 rows
  bins = IRanges(start = breaks[-length(breaks)], end = breaks[-1]) # number of segments
  
  # x and y correspond to the end point of each interval in ir
  hc@BINS = bins
  
  # position for two end points for every segment (e.g. interval in bins)
  if(!is.null(sfc_pos)){
    # If a SFC positions are given, swap them out with the default HilbertCurve 
    # positions to allow a different SFC:
    assert_that(is.data.frame(sfc_pos), 
                all(colnames(sfc_pos) %in% c('x1', 'x2', 'y1', 'y2')),
                msg="This is a copy of HilbertCurve() function and requires a mapping to
                visualize. As such, you need to feed in a dataframe of x1,x2,y1,y2 coords")
    hc@POS = sfc_pos[,c('x1', 'y1', 'x2', 'y2')]
  } else {
    # Use the default HilbertCurve code for position mapping
    hc@POS = data.frame(x1 = pos$x[-n], y1 = pos$y[-n], x2 = pos$x[-1], y2 = pos$y[-1])
  }
  
  hc@RANGE = c(breaks[1], breaks[length(breaks)])
  hc@LEVEL = level
  
  mode = match.arg(mode)[1]
  hc@MODE = mode
  
  
  if(newpage) grid.newpage()
  HilbertCurve:::increase_plot_index()
  
  # create a 2x2 layout
  if(length(title) == 0) {
    title_height = unit(0, "mm")
  } else {
    title_height = grobHeight(textGrob(title, gp = title_gp)) + unit(5, "mm")
  }
  
  .width = function(x) {
    if(inherits(x, "grob")) {
      grobWidth(x)
    } else if(inherits(x, "Legends")) {
      ComplexHeatmap:::width(x)
    } else {
      stop("Class not supported.")
    }
  }
  .height = function(x) {
    if(inherits(x, "grob")) {
      grobHeight(x)
    } else if(inherits(x, "Legends")) {
      ComplexHeatmap:::height(x)
    } else {
      stop("Class not supported.")
    }
  }
  .draw = function(x) {
    if(inherits(x, "grob")) {
      grid.draw(x)
    } else if(inherits(x, "Legends")) {
      ComplexHeatmap::draw(x)
    } 
  }
  
  if(length(legend) == 0) {
    legend_width = unit(0, "mm")
    legend_height = unit(0, "mm")
  } else {
    if(inherits(legend, "Legends")) {
      legend = list(legend)
    } else if(inherits(legend, "grob")) {
      legend = list(legend)
    } else if(inherits(legend, "list")) {
      
    } else {
      stop("`legend` should be a single legend or a list of legends.")
    }
    
    legend_width = max(do.call("unit.c", lapply(legend, .width))) + unit(4, "mm")
  }
  
  pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2, widths = unit.c(unit(1, "npc") - legend_width, legend_width), 
                                             heights = unit.c(title_height, unit(1, "npc") - title_height)),
                        name = paste0("hilbert_curve_", HilbertCurve:::get_plot_index(), "_global")))
  
  title_gp = HilbertCurve:::validate_gpar(title_gp, default = list(fontsize = 16))
  if(length(title) != 0) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.text(title, gp = title_gp)
    upViewport()
  }
  if(length(legend) != 0) {
    gap = unit(2, "mm")
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    # draw the list of legend
    legend_height = sum(do.call("unit.c", lapply(legend, .height))) + gap*(length(legend)-1)
    y = unit(0.5, "npc") + legend_height*0.5 
    for(i in seq_along(legend)) {
      pushViewport(viewport(x = unit(2, "mm"), y = y, height = .height(legend[[i]]), width = .width(legend[[i]]), just = c("left", "top")))
      .draw(legend[[i]])
      upViewport()
      y = y - gap - .height(legend[[i]])
    }
    
    upViewport()
  }
  
  size = unit(1, "snpc") - padding*2
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  
  pushViewport(viewport(name = paste0("hilbert_curve_", HilbertCurve:::get_plot_index()), xscale = c(-0.5, 2^level - 0.5), yscale = c(-0.5, sqrt(n)-0.5), width = size, height = size))
  grid.rect(gp = gpar(fill = background_col, col = background_border))
  
  reference_gp = HilbertCurve:::validate_gpar(reference_gp, default = list(lty = 3, col = "#999999"))
  if(mode == "normal") {
    if(reference) {
      grid.segments(hc@POS$x1, hc@POS$y1, hc@POS$x2, hc@POS$y2, default.units = "native", gp = reference_gp)
      
      # grid.points(hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))
      # grid.points(hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))
      
      # grid.text(round(unzoom(hc, start(bins)[1]), 2), hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#999999", fontsize = 10))
      # grid.text(round(unzoom(hc, end(bins)), 2), hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#999999", fontsize = 10))
      
      if(arrow) HilbertCurve:::grid_arrows(hc@POS$x1, hc@POS$y1, (hc@POS$x1+hc@POS$x2)/2, (hc@POS$y1+hc@POS$y2)/2, only.head = TRUE, arrow_gp = gpar(fill = reference_gp$col, col = NA))
    }
    upViewport(3)
  } else {
    background_col = background_col[1]
    background_col = col2rgb(background_col) / 255
    red = matrix(background_col[1], nrow = 2^level, ncol = 2^level)
    green = matrix(background_col[2], nrow = 2^level, ncol = 2^level)
    blue = matrix(background_col[3], nrow = 2^level, ncol = 2^level)
    hc@RGB$red = red
    hc@RGB$green = green
    hc@RGB$blue = blue
    
    HilbertCurve:::add_raster(hc@RGB)  # already jump to top vp
    
  }
  
  # message(strwrap("If your regions are mostly very small, consider to set 'mean_mode' to 'absolute' when using `hc_points()`, `hc_rect()` and `hc_layer()`."))
  
  return(invisible(hc))
}