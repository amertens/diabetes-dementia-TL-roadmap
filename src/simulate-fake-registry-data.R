


#devtools::install_github("tagteam/heaven")

library(tidyverse)
library(heaven)
vignette('user-heaven')
browseVignettes('heaven')

# population
heaven::simPop(10)
# hospital admission
heaven::simAdmissionData(10)
# purchase of medicine
df<-heaven::simPrescriptionData(10)  

table(df$atc)

# http://medinfo.dk/sks/brows.php 
# http://www.medicinpriser.dk/ 
#   http://pro.medicin.dk/ 

  
  #Lookup codes:
  #https://www.whocc.no/atc_ddd_index/


#heaven::medicinMacro()