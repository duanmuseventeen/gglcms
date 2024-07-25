#' Mass Spectrometry Dissociation Graphic
#'
#' @description This function draws plots of mass spectrometry dissociation based on ggplot2 system.
#'
#' @param dat a data frame, containing 4 columns named Compound,	m/z, Relative Intensity, Color.
#' @param img a character refers to the path of a picture of dissociation.
#' @param img.size size of dissociation picture.
#' @param heights a vector refers to heights of picture and dissociation graphic, respectively.
#' @param colors a vector used to customize color for refer and selected.
#'
#' @return a ggplot element
#' @export
#'
#' @examples
#' geom_MS2(MS2_sample)
geom_MS2 <- function(dat, img = "", img.size = 30, heights = c(1,1), colors = c("#222222", "#C4121A")){
  p <- ggplot(data = dat) +
    geom_segment(aes(x = `m/z`, y = 0, xend = `m/z`, yend = `Relative Intensity`, color = Color)) +
    ggrepel::geom_text_repel(aes(x = `m/z`, y = `Relative Intensity`, label = `m/z`, color = Color), max.overlaps = 15, vjust = 1) +
    labs(x = "m/z", y = "Relative Intensity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = colors) +
    guides(color = "none") +
    theme_classic() +
    theme(text = element_text(size = 12),
          plot.margin = ggplot2::margin(30,30,30,30),
          plot.title = element_text(hjust = .5)
    )
  if(img != ""){
    p2 <- ggplot() +
      geom_point_img(aes(x = 1, y = 1, img = img), size = img.size) +
      theme_void()
    p <- p2 + p + plot_layout(ncol = 1, heights = heights)
  }
  p
}
