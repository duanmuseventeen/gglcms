#' Compare MS1 vs Reference Library
#'
#' @description This function compares MS1 derived from R package xcms vs reference library.
#'
#' @param REF directory of reference library file
#' @param Sample directory of sample file
#' @param Filter_Adduct the types of adduct ions should be considered. if "all", a total of 52 adduct ions will be computed
#' @param Filter_Charge the charges of adduct ions should be considered. if "all", all of adduct ions will be computed
#' @param error accepted error between reference and sample. Da by default.
#' @param use_ppm logical. if FALSE, the unit of error is Da, if TRUE, the unit of error is ppm.
#' @param atom_mass a data frame includes exact atomic mass. Average Molecular Weight derived from https://iupac.qmul.ac.uk/AtWt/ was used by default
#'
#' @return a list
#' @export
#'
#' @examples
#' no examples
stat_MS1 <-
  function(REF = "reference_path",
           Sample = "sample_path",
           Filter_Adduct = c("M+H","M-2H2O+H","M-H2O+H","M+Li","M+NH4","M+Na"),
           Filter_Charge = c("1+","2+","3+"),
           error = 0.1, use_ppm = F,
           atom_mass = annot){

    REF <- readxl::read_excel(REF) %>%
      filter(!duplicated(CAS))
    Sample <- data.table::fread(Sample)

    result <- data.frame(
      mzmed = NA,
      rtmed = NA,
      Compound = NA,
      Formula = NA,
      CAS = NA,
      ADDUCT = NA,
      `m/z` = NA,
      CHARGE = NA,
      ppm = NA,
      stringsAsFactors = F,
      check.names = F
    )
    n <- 1
    # require(progress)
    # pb <- progress_bar$new(total = nrow(Sample))
    for (i in 1:nrow(Sample)) {
      mzmed <- Sample$mzmed[i]
      for (j in 1:nrow(REF)) {
        if(any(Filter_Adduct == "all") & any(Filter_Charge == "all")){
          dict <- stat_mol_weight(REF$Formula[j], atom_mass = atom_mass, adduct = TRUE)[[1]]
        }else if(any(Filter_Adduct == "all") & any(Filter_Charge != "all")){
          dict <- stat_mol_weight(REF$Formula[j], atom_mass = atom_mass, adduct = TRUE)[[1]] %>%
            filter(CHARGE %in% Filter_Charge)
        }else if(any(Filter_Adduct != "all") & any(Filter_Charge == "all")){
          dict <- stat_mol_weight(REF$Formula[j], atom_mass = atom_mass, adduct = TRUE)[[1]] %>%
            filter(ADDUCT %in% Filter_Adduct)
        }else if(any(Filter_Adduct != "all") & any(Filter_Charge != "all")){
          dict <- stat_mol_weight(REF$Formula[j], atom_mass = atom_mass, adduct = TRUE)[[1]] %>%
            filter(ADDUCT %in% Filter_Adduct) %>%
            filter(CHARGE %in% Filter_Charge)
        }

        for (k in dict$ADDUCT) {
          delta <- abs(mzmed - dict$`ADDUCT M/Z`[dict$ADDUCT == k])
          ppm <- round((delta/ dict$`ADDUCT M/Z`[dict$ADDUCT == k]) * 10^6)
          flag <- ifelse(use_ppm, ppm < error, delta < error)
          if(flag){
            result[n,1:9] <- c(mzmed, Sample$rtmed[i],
                               REF$Compound[j],
                               REF$Formula[j],
                               REF$CAS[j],
                               k,
                               dict$`ADDUCT M/Z`[dict$ADDUCT == k],
                               dict$CHARGE[dict$ADDUCT == k],
                               ppm)
            n <- n + 1
          }
        }
      }
      # pb$tick()
    }

    return(result)
  }

