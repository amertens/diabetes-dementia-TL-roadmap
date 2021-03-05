

library(tidyverse)
library(nycflights13)

if(!require(DataExplorer)){install.packages("DataExplorer")}
library(DataExplorer)
data_list <- list(airlines, airports, flights, planes, weather)
plot_str(data_list)

merge_airlines <- merge(flights, airlines, by = "carrier", all.x = TRUE)
merge_planes <- merge(merge_airlines, planes, by = "tailnum", all.x = TRUE, suffixes = c("_flights", "_planes"))
merge_airports_origin <- merge(merge_planes, airports, by.x = "origin", by.y = "faa", all.x = TRUE, suffixes = c("_carrier", "_origin"))
final_data <- merge(merge_airports_origin, airports, by.x = "dest", by.y = "faa", all.x = TRUE, suffixes = c("_origin", "_dest"))

introduce(final_data)


plot_intro(final_data)


plot_missing(final_data)

final_data <- drop_columns(final_data, "speed")

plot_bar(final_data)

final_data[which(final_data$manufacturer == "AIRBUS INDUSTRIE"),]$manufacturer <- "AIRBUS"
final_data[which(final_data$manufacturer == "CANADAIR LTD"),]$manufacturer <- "CANADAIR"
final_data[which(final_data$manufacturer %in% c("MCDONNELL DOUGLAS AIRCRAFT CO", "MCDONNELL DOUGLAS CORPORATION")),]$manufacturer <- "MCDONNELL DOUGLAS"
final_data <- drop_columns(final_data, c("dst_origin", "tzone_origin", "year_flights", "tz_origin"))

plot_bar(final_data, with = "arr_delay")
plot_bar(final_data, by = "origin")


plot_histogram(final_data)

plot_correlation(na.omit(final_data), maxcat = 5L)

plot_correlation(na.omit(final_data), type = "c")
plot_correlation(na.omit(final_data), type = "d")

pca_df <- na.omit(final_data[, c("origin", "dep_delay", "arr_delay", "air_time", "year_planes", "seats")])
plot_prcomp(pca_df, variance_cap = 0.9, nrow = 2L, ncol = 2L)


create_report(final_data)

#To maximize the usage of this function, always supply a response variable 
#(if applicable) to automate various bivariate analyses. For example,
create_report(final_data, y = "arr_delay")
