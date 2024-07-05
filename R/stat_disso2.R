#' Compute Cosine Similarity of two Mass Spectrometry Dissociation Graphic
#'
#' @param ref data frame of reference MS dissocaiation.
#' @param sample data frame of sample MS dissocaiation.
#' @param ref.cutoff Remove fragments whose relative intensity under cut-off.
#' @param ppm To screen the same fragment of reference in sample, 10 by default.
#'
#' @return a number between -1 and 1.
#' @export
#'
#' @examples
#' stat_disso2(disso_reference, disso_sample, ref.cutoff = 30, ppm = 10000)
stat_disso2 <- function(ref, sample, ref.cutoff = 10, ppm = 10){
  ref <- ref %>% filter(`Relative Intensity` > ref.cutoff)

  vsample <- sample$`Relative Intensity`[
    lapply(ref$`m/z`, function(x){
      index <- which(abs((x - sample$`m/z`)) < (ppm * x / 1000000))
      if(length(index) == 0){index <- NA}
      return(index)
    }) %>% unlist]
  vsample[is.na(vsample)] <- 0
  return(lsa::cosine(ref$`Relative Intensity`,vsample) %>% as.numeric)
}
