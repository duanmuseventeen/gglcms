#' Compute similarity of two chromatographs
#'
#' @description To compute similarity of two chromatographs
#'
#' @param a Area or points of Chromatograph 1
#' @param b Area or points of Chromatograph 2
#' @param Na Number of peaks in Chromatograph 1
#' @param Nb Number of peaks in Chromatograph 2
#' @param Nab Number of peaks in both Chromatograph 1 and 2
#' @param method pearson, spearman, cosine, Overlap and Tanimoto methods are currently available
#'
#' @return a numeric value.
#' @export
#'
#' @examples
#' stat_simi(a,b, method = "pearson")
#' sata_simi(Na,Nb,Nab, method = "Tanimoto")
stat_simi <- function(a,b,Na,Nb,Nab,
                      method = "pearson"
){
  # c("pearson","spearman","cosine","Overlap","Tanimoto")
  if(method == "pearson"){
    cor.test(a,b,method = "pearson")
  }else if(method == "spearman"){
    cor.test(a,b,method = "spearman")
  }else if(method == "cosine"){
    lsa::cosine(a,b)
  }else if(method == "Overlap"){
    2*Nab/(Na+Nb)
  }else if(method == "Tanimoto"){
    Nab/(Na+Nb-Nab)
  }
}
