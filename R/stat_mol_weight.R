#' Compute molecular weight
#'
#' @param ... a vector or list of molecular formulas in character format
#'
#' @return list
#' @export
#'
#' @examples
#' stat_mol_weight(c("CH3COOH","NaCl","H2O","2H2O","90745","/"))
#'
#' stat_mol_weight(list(a = "CH3COOH", b = c("CH3CH2OH","KOH","Mg2+","[CH3COO]-")))

stat_mol_weight <- function(...){

  # library(rvest)
  #
  # urlb <- "https://iupac.qmul.ac.uk/AtWt/"
  #
  # annot <- read_html(urlb) %>%
  #   html_table(header = TRUE) %>% .[[2]] %>%
  #   mutate(`Atomic Wt` = `Atomic Wt` %>% str_remove_all(" ")) %>%
  #   mutate(`Atomic Wt` = `Atomic Wt` %>% str_remove_all("\\([0-9]+\\)")) %>%
  #   mutate(`Atomic Wt` = `Atomic Wt` %>% str_remove_all("\\[")) %>%
  #   mutate(`Atomic Wt` = `Atomic Wt` %>% str_remove_all("\\]")) %>%
  #   mutate(`Atomic Wt` = `Atomic Wt` %>% as.numeric) %>%
  #   mutate(nchar = nchar(Symbol)) %>%
  #   arrange(desc(nchar))

  annot <- data.frame(
    `At No` = c(0,0),
    Symbol = c("+","-"),
    Name = c("positive electron","negative electron"),
    `Atomic Wt` = c(0.0005489, -0.0005489),
    Notes = c("",""),
    nchar = c(1,1),
    stringsAsFactors = F,
    check.names = F
  ) %>%
    dplyr::bind_rows(
      annot
    )

  list <- list(...) %>% unlist

  lapply(list, function(x){
    raw <- x
    mw <- 0
    for (i in annot$Symbol) {
      if(i %in% c("+","-")){
        tmp <- str_extract_all(x, paste0('[0-9]*\\',i)) %>% unlist
        x <- str_replace_all(x, paste0('[0-9]*\\',i),"")
      }else{
        tmp <- str_extract_all(x, paste0(i,'[0-9]*')) %>% unlist
        x <- str_replace_all(x, paste0(i,'[0-9]*'),"")
      }

      if(length(tmp) != 0){
        num <- lapply(tmp, function(y){
          num <- str_extract_all(y, "[0-9]+") %>% unlist %>% as.numeric()
          num <- ifelse(length(num) == 0, 1, num)
          return(num)
        }) %>% do.call(sum, .)
        mw <- mw + num * annot$`Atomic Wt`[annot$Symbol == i]
      }
    }
    if(!(x %in% c("","[]","·"))) message(paste0("Unknown character, results for ",raw," may be not correct!"))
    return(mw)
  }) -> result
  result
}
