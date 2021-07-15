simPrescriptionData_atc <- function (n, max.prescriptions = 37, packages = list(list(c(200, 
                                                            30), c(400, 100), c(400, 300), c(500, 60)), list(c(750, 100), 
                                                                                                             c(750, 250), c(75, 500))), max.packages = 3, startDate = "1995-01-01") 
{
  pnr = eksd = NULL
  startDate <- as.Date(startDate)
  if (is.null(names(packages))) {
    atc <- c("R06AC01", "N01BB51", "D08AX01", "N02AA03", 
             "C09BA02", "B02BC", "L01AD02", "C08CA08", "L01AA02", 
             "J07BB", "R05CB15", "C09DA04", "L04AA32", "J05AE03", 
             "N04BD03", "M03AX01", "A12AX", "D06BB12", "R03DA03", 
             "A06AB06", "L03AA10", "N07BB05", "A09AA02", "G02AD04", 
             "H05BA01", "C03CB02", "C01EB17")
    glp1 <- c("A10BJ01", 
    "A10BJ02", 
    "A10BJ03", 
    "A10BJ04",  
    "A10BJ05",  
    "A10BJ06",  
    "A10BJ07")  
    dpp4 <- c("A10BH01",  
    "A10BH02",  
    "A10BH03",  
    "A10BH04",  
    "A10BH05",  
    "A10BH06",  
    "A10BH07",  
    "A10BH08",  
    "A10BH51",  
    "A10BH52")
    SGLT2i <- c("A10BK01", 
    "A10BK02",  
    "A10BK03",  
    "A10BK04",  
    "A10BK05",  
    "A10BK06",  
    "A10BK07")
    hypertension <- unlist(hypertensionATC)
    
    # dementia code, atc and ICD-10
    dementia <- c("N06DA01", "N06DA02", "N06DA03", "N06DA04", "N06DA05", "N06DA52", "N06DA53")
    #F01.50, F01.51, F02.80, F02.81, F03.90, F03.91, G30.0, G30.1, G30.8, G30.9, G31.01, G31.09, G31.1, G31.2, G31.81, G31.82, G31.83, G31.84, G31.85, G31.89, G31.9, G32.0
    
    atc <- c(atc, glp1, dpp4, SGLT2i, dementia, hypertension)
    
  }
  else {
    atc <- names(packages)
  }
  out <- data.table::rbindlist(lapply(1:n, function(i) {
    pat.i <- data.table::rbindlist(lapply(1:length(packages), 
                                          function(p) {
                                            pack <- unlist(packages[p], recursive = FALSE)
                                            M = sample(1:max.prescriptions, size = 1)
                                            sizes <- sapply(pack, "[", 2)
                                            if (length(sizes) == 1) 
                                              sizes <- rep(sizes, 2)
                                            else {
                                              if (M == 1) 
                                                sizes <- rep(sizes, 2)
                                              else sample(sizes, size = M, replace = TRUE)
                                            }
                                            strengths <- sapply(pack, "[", 1)
                                            if (length(strengths) == 1) 
                                              strengths <- rep(strengths, 2)
                                            else {
                                              if (M == 1) 
                                                strengths <- rep(strengths, 2)
                                              else strengths <- sample(strengths, size = M, 
                                                                       replace = TRUE)
                                            }
                                            out <- data.table::data.table(eksd = startDate + 
                                                                            rbinom(M, 1, 0.95) * floor(runif(M, 0, 5 * 
                                                                                                               365.25)), atc = sample(atc, size = M, replace = TRUE), 
                                                                          packsize = sample(sizes, size = M, replace = TRUE), 
                                                                          apk = sample(1:max.packages, size = M, replace = TRUE), 
                                                                          strnum = sample(strengths, size = M, replace = TRUE))
                                            out
                                          }))
    pat.i[, `:=`(pnr, i)]
  }))
  data.table::setkey(out, pnr, atc, eksd)
  data.table::setcolorder(out, c("pnr", "atc", "eksd", "strnum", 
                                 "packsize", "apk"))
  out
}