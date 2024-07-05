#' Compare two Mass Spectrometry Dissociation Graphic
#'
#' @param data data frame, containing 5 columns named Compound,	m/z, Relative Intensity, Color and group. Color is grouped by "refer" and "selected", and group is grouped by "sample" and "refer".
#' @param img a character refers to the path of a picture of dissociation.
#' @param img.size size of dissociation picture.
#' @param heights a vector refers to heights of picture and dissociation graphic, respectively.
#' @param cutoff a threshold to filter fragment whose relative intensity under it
#' @param colors a vector used to customize color for sample-refer, sample-selected, refer-refer and refer-selected.
#'
#' @return a ggplot element
#' @export
#'
#' @examples
#' geom_disso2(ms_disso_compar, cutoff = 20)
geom_disso2 <- function(dat, img = "", img.size = 30, heights = c(1,1),
                        cutoff = 0,
                        colors = c("#222222", "#C4121A","#006eb3", "#C4121A")){
  dat <- dat %>%
    filter(`Relative Intensity` > cutoff) %>%
    mutate(`Relative Intensity` = ifelse(group == "reference",-`Relative Intensity`,`Relative Intensity`),
           Color = paste0(Color,"-",group)) %>%
    mutate(Color = factor(Color, levels = c("common-sample","feature-sample","common-reference","feature-reference")))
  p <- ggplot(data = dat) +
    geom_hline(yintercept = 0, color = "gray20") +
    geom_segment(aes(x = `m/z`, y = 0, xend = `m/z`, yend = `Relative Intensity`, color = Color)) +
    ggrepel::geom_text_repel(aes(x = `m/z`, y = `Relative Intensity`, label = `m/z`, color = Color), max.overlaps = 15, vjust = 1) +
    labs(x = "m/z", y = "Relative Intensity") +
    scale_y_continuous(expand = c(0,0), limits = c(-100,100)) +
    scale_color_manual(values = colors) +
    guides(color = "none") +
    theme_classic() +
    theme(text = element_text(size = 12),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
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
