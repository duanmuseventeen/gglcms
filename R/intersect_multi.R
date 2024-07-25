#' Intersect multiple elements
#'
#' @param list a list of elements to intersect
#'
#' @return in terms of input return vector or data.frame
#' @export
#'
#' @examples
#' intersect_multi(list = list(c(1:3),c(1:4),c(2:4)))
#'
#' intersect_multi(list = list(disso_sample, disso_sample, disso_sample))
#'
intersect_multi <- function(list){
  if(length(list) == 1){
    return(list)
  }else if(length(list) == 2){
    return(intersect(list[[1]], list[[2]]))
  }else  if(length(list) > 2){
    i <- length(list)
    list[[i-1]] <- intersect(list[[i]], list[[i-1]])
    list <- list[-i]
    intersect_multi(list)
  }
}


