#Prepare directories as needed and load in data from our GEO submission ----

dir.create(c("raw_data", "data_objects", "plots", "tabular_output", "chooseR",
             "cross_entropy_test", "cluster_annotation"))

for (i in seq_along(1:8)){
  mkdir(paste0("raw_data/sample", i))
  
  #add data download instructions here
}

