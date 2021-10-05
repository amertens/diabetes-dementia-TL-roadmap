# result_folder <- "/home/leo42k/Documents/projects/brc_output/ltmle/test/temp_output"
# result_folder <- paste0("/home/leo42k/Documents/projects/brc_output/202106/parallel_ltmle_202108_5000_glm", "/temp_output")
library(abind)
library(xtable)

# type <- "glm"
type <- "default"
# type <- "SL.hal9001"



identifier <- "ltmle_202108"
identifier_vec <- c()
size <- 5000; identifier_vec[1] <- paste0(identifier, "_", size, "_", type)
size <- 20000; identifier_vec[2] <- paste0(identifier, "_", size, "_", type)
size <- 200000; identifier_vec[3] <- paste0(identifier, "_", format(size, scientific = F), "_", type)

K <- 10

report_table <- list()

for (s in 1:3) {
  # result_folder <- paste0("/home/leo42k/Documents/projects/brc_output/202106/tmle3/parallel_", identifier_vec[s], "/temp_output")
  result_folder <- paste0("/home/leo42k/Documents/projects/brc_output/ltmle/parallel_", identifier_vec[s], "/temp_output")
  all_results <- lapply(1:K, function(i) {
    if(file.exists(paste0(result_folder, "/", i, ".RDS"))) {
      paste0(result_folder, "/", i, ".RDS") %>% readRDS %>% return
    } else {
      return(NULL)
    }
  })
  
  vec_time <- lapply(all_results, function(result) result$time_diff) %>% unlist
  vec_mem <- lapply(all_results, function(result) result$mem_diff) %>% unlist
  vec_mem <- vec_mem/1000000
  
  temp_table <- rbind(
    lapply(all_results, function(result) result$output$estimates[1]) %>% unlist %>% quantile(), 
    vec_time %>% quantile, 
    vec_mem %>% quantile  
  )
  
  report_table[[s]] <- temp_table
}

final_table <- report_table %>% lapply(function(u) u[2:3, ]) %>% abind(along = 1)
rownames(final_table) <- expand.grid(c("time", "mem"), c(5000, 20000, 200000)) %>% 
  apply(1, function(u) paste0(u, collapse = ", "))
final_table
round(final_table, 3) %>% xtable
