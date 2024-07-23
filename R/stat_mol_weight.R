#' Compute molecular weight
#'
#' @param ... a vector or list of molecular formulas in character format
#' @param atom_mass a data frame includes exact atomic mass
#' @param adduct If TRUE, compute molecular weight of adduct ions
#'
#' @return list
#' @export
#'
#' @examples
#' stat_mol_weight(c("CH3COOH","NaCl","H2O","2H2O","90745","/"))
#'
#' stat_mol_weight(list(a = "CH3COOH", b = c("CH3CH2OH","KOH","Mg2+","[CH3COO]-")))

stat_mol_weight <- function(..., atom_mass = annot, adduct = F){

  # 获取精确原子质量
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

  # 获取加和离子
  # https://hmdb.ca/spectra/ms/search
  # ADDUCT <- data.frame(
  #   ADDUCT = c("M+H",	"M-2H2O+H",	"M-H2O+H", "M-H2O+NH4", "M+Li", "M+NH4", "M+Na", "M+CH3OH+H",
  #              "M+K", "M+ACN+H", "M+2Na-H", "M+IsoProp+H", "M+ACN+Na", "M+2K-H", "M+DMSO+H",
  #              "M+2ACN+H", "M+IsoProp+Na+H", "2M+H", "2M+NH4", "2M+Na", "2M+3H2O+2H", "2M+K",
  #              "2M+ACN+H", "2M+ACN+Na", "M+2H", "M+H+NH4", "M+H+Na", "M+H+K", "M+ACN+2H", "M+2Na",
  #              "M+2ACN+2H", "M+3ACN+2H", "M+3H", "M+Na+2H", "M+2Na+H", "M+3Na", "M-H", "M-H2O-H",
  #              "M+F", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M+HAc-H", "M+Br", "M+TFA-H",
  #              "2M-H", "2M+FA-H", "2M+HAc-H", "3M-H", "M-2H", "M-3H"),
  #   CHARGE = c("1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+","1+",
  #             "1+","1+","1+","1+","1+","1+","1+","2+","2+","2+","2+","2+","2+","2+","2+","3+","3+",
  #             "3+","3+","1-","1-","1-","1-","1-","1-","1-","1-","1-","1-","1-","1-","1-","1-","2-","3-"),
  #   stringsAsFactors = FALSE
  # ) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT, "IsoProp", "C3H8O")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "2ACN", "ACNACN")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "3ACN", "ACNACNACN")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "ACN", "CH3CN")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "2H2O", "H4O2")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "3H2O", "H6O3")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "2H", "H2")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "3H", "H3")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "2Na", "Na2")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "3Na", "Na3")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "2K", "K2")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "ACN", "CH3CN")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "TFA", "CF3COOH")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "FA", "HCOOH")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "HAc", "CH3COOH")) %>%
  #   mutate(ADDUCT_FULL = str_replace_all(ADDUCT_FULL, "DMSO", "C2H6OS")) %>%
  #   mutate(ADDUCT_FULL = str_remove_all(ADDUCT_FULL, "[0-9]{0,1}M")) %>%
  #   mutate(ADDUCT_FULL = str_remove_all(ADDUCT_FULL, "[\\+\\-]")) %>%
  #   mutate(M = str_extract(ADDUCT, "[0-9]{0,1}M")) %>%
  #   mutate(M = ifelse(M == "M", 1, as.numeric(str_remove_all(M, "M")))) %>%
  #   mutate(`M/Z` = unlist(stat_mol_weight(ADDUCT_FULL)) + unlist(stat_mol_weight(CHARGE)))
  #
  # annot <- data.frame(
  #   `At No` = c(0,0),
  #   Symbol = c("+","-"),
  #   Name = c("positive electron","negative electron"),
  #   `Atomic Wt` = c(0.0005489, -0.0005489),
  #   Notes = c("",""),
  #   nchar = c(1,1),
  #   stringsAsFactors = F,
  #   check.names = F
  # ) %>%
  #   dplyr::bind_rows(
  #     annot
  #   )

  list <- list(...) %>% unlist

  lapply(list, function(x){
    raw <- x
    mw <- 0
    for (i in atom_mass$Symbol) {
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
        mw <- mw + num * atom_mass$`Atomic Wt`[atom_mass$Symbol == i]
      }
    }
    if(!(x %in% c("","[]","·"))) message(paste0("Unknown character, results for ",raw," may be not correct!"))
    return(mw)
  }) -> result

  names(result) <- list

  if(adduct){
    result <- lapply(result, function(x){
      ADDUCT <- ADDUCT %>%
        mutate(MW = x) %>%
        mutate(`ADDUCT M/Z` = MW*M + `M/Z`)
    })
  }

  result
}
