#' Display Chromatographs from LC-MS
#'
#' @param data data frame. Contain X(min), Y(Counts), batch and Compound.
#' @param legend.position the default position of legends ("none", "left", "right", "bottom", "top", "inside"). A numeric vector of length two setting the placement of legends that have the "inside" position.
#' @param max.over for group label. 30 by default. Only available when number of group is smaller than 30.
#' @param guide.ncol column number of legend.
#' @param xlimit limits for x axis
#' @param x.nbreaks An integer guiding the number of major breaks for x axis. 5 by default.
#' @param label label for group.
#' @param batchlabel label for batch.
#' @param grid logical. Should grid be displayed?
#' @param shift logical. Shift chromatographs from different batches. Available when batch more than 1.
#' @param angle angle of shift.
#' @param percent minor degree for shift
#' @param integration logical. Should area of chromatographs be displayed?
#'
#' @return a ggplot element
#' @export
#'
#' @examples
#' tmp <- chromatograph %>% filter(batch %in% c("HB-1","HB-2"))
#' geom_chromat(tmp, shift = F, angle = 45)
#'
#' geom_chromat(tmp, shift = T, angle = 45, percent = 0.1)
geom_chromat <- function(data, legend.position = "none",
                         max.over = 30,
                         guide.ncol = 1,
                         xlimit = NA, x.nbreaks = 5,
                         label = F,
                         # xlab = "Retention Time (min)",
                         # ylab = "Intensity (cps)",
                         batchlabel = T, grid = F, shift = F, angle = 0, percent = 0.05,
                         integration = F){

  data <- data %>%
    mutate(group = paste0(batch,"-",Compound)) %>%
    mutate(group = factor(group))

  if(any(is.na(xlimit))){ xlimit = range(data$`X(min)`)
  }else{
    data <- data %>%
      filter(`X(min)` > xlimit[1] & `X(min)` < xlimit[2])
  }

  if(length(unique(data$group)) > 30) legend.position <- "none"

  if(shift){
    # 计算偏移
    ymax = max(data$`Y(Counts)`)

    xmax = max(data$`X(min)`)

    yshift = percent * ymax * cos(angle / 180 * pi)

    xshift = percent * xmax * sin(angle / 180 * pi)

    if(grid){
      grid.x <- seq(xlimit[1],xlimit[2],length.out = 5)
      grid.xend <- grid.x + (length(unique(data$batch)) - 1) * xshift
      grid.y <- rep(0, 5)
      grid.yend <- grid.y + (length(unique(data$batch)) - 1) * yshift
    }

    # 赋值
    data <- data %>%
      left_join(data.frame(
        batch = unique(data$batch),
        no = c(0:(length(unique(data$batch)) - 1)),
        stringsAsFactors = F),
        by = "batch") %>%
      mutate(x = no*xshift + `X(min)`,
             y = no*yshift + `Y(Counts)`) %>%
      group_by(group) %>%
      mutate(ymin = min(y))
    xlimit[2] <- xlimit[2] + (length(unique(data$batch)) - 1) * xshift

    data <- data

    # 作图
    p <- ggplot(data, aes(x = x, y = y, color = group, fill = group))
    if(grid){
      p <- p +
        annotate(geom = "segment", x = grid.x, xend = grid.xend, y = grid.y, yend = grid.yend, color = "gray90")
    }
    p <- p + geom_line()
    if(integration){
      p <- p +
        geom_ribbon(aes(ymin = ymin, ymax = y), alpha = 0.2)
    }
    if(label){
      data.label <- data %>%
        group_by(group) %>%
        filter(y == max(y))
      p <- p +
        ggrepel::geom_text_repel(data = data.label, max.overlaps = max.over,
                                 aes(x = x, y = y,label = group),
                                 show.legend=FALSE)
    }
    if(batchlabel){
      data.batch.label <- data %>%
        group_by(batch) %>%
        arrange(desc(x)) %>%
        filter(!duplicated(batch))
      p <- p +
        geom_text(data = data.batch.label, aes(x = x, y = y,label = batch),
                  show.legend=FALSE, col = "black", vjust = -0.1, hjust = -0.3)
    }
    p <- p +
      scale_x_continuous(limits = xlimit, n.breaks = x.nbreaks)
    if(legend.position == "none"){
      p <- p + guides(color = "none", fill = "none")
    }else{
      p <- p +
        guides(colour = guide_legend(ncol = guide.ncol), fill = "none"
               #shape = guide_legend(override.aes = list(lintype = "blank", size= 6))
        )
    }
    p <- p +
      # labs(x = xlab, y = ylab) +
      coord_cartesian(clip = "off") +
      theme_classic() +
      theme(text = element_text(size = 16),
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_blank(),
            plot.margin = ggplot2::margin(40,40,40,40),
            panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.position = legend.position)


  }else{
    data.label <- data %>%
      group_by(group) %>%
      filter(`Y(Counts)` == max(`Y(Counts)`))

    # 作图
    p <- ggplot(data, aes(x = `X(min)`, y = `Y(Counts)`, color = group)) +
      geom_line()
    if(integration){
      p <- p +
        geom_area(fill = group, alpha = 0.2)
    }
    if(label){
      p <- p +
        ggrepel::geom_text_repel(data = data.label, max.overlaps = max.over,
                                 aes(x = `X(min)`, y = `Y(Counts)`,label = group),
                                 show.legend=FALSE)
    }
    p <- p +
      scale_x_continuous(expand = c(0,0), n.breaks = x.nbreaks)
    if(legend.position == "none"){
      p <- p +
        guides(color = "none", fill = "none")
    }else{
      p <- p +
        guides(colour = guide_legend(ncol = guide.ncol), fill = "none"
               #shape = guide_legend(override.aes = list(lintype = "blank", size= 6))
        )
    }
    p <- p +
      # labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(text = element_text(size = 16),
            plot.margin = ggplot2::margin(40,40,40,40),
            panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.position = legend.position)
  }
  p
}
