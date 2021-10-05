setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

library(dplyr)
library(data.table)
library(ltmle)

dt_tmle <- readRDS("./simulation_ZW/dt_tmle_20210712.rds")

node_names <- c("age", "sex", 
                # "first_date_2nd_line", 
  expand.grid(c("C_", "A1_", "L_", "Y_"), 1:5) %>% apply(1, function(row) paste0(row, collapse = ""))
)
dt_use <- dt_tmle[, ..node_names]

for (x in grep("Y_", node_names)) dt_use[, (node_names[x]) := 1 - get(node_names[x])]
# for (x in grep("C_", node_names)) dt_use[, (node_names[x]) := 1 - get(node_names[x])]
ltmle(dt_use[sample(nrow(dt_tmle), 1000000, T), ], Anodes = grep("^A1_", node_names), Lnodes = grep("^L_", node_names), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
      abar = rep(1, 5))
