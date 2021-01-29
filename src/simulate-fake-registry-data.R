


#devtools::install_github("tagteam/heaven")

library(tidyverse)
library(heaven)
vignette('user-heaven')

# population
heaven::simPop(10)
# hospital admission
heaven::simAdmissionData(10)
# purchase of medicine
heaven::simPrescriptionData(10)  