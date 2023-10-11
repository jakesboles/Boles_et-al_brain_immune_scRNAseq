dir.create("01_raw_data")

for (i in seq_along(1:8)){
  mkdir(paste0("01_raw_data/sample", i))
}