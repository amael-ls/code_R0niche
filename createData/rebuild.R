
#### Aim of prog: Reconstruct the sliced tree database from competition folder
# In order to fasten the computation of competition, I sliced the tree database.
# I now rebuild it from its slices

#### Load library
library(data.table)

#### Read data
ls_files = list.files("./competition")
(nbFiles = length(ls_files))

ls_slices = vector(mode = "list", length = nbFiles)

#### Rebuild
for (i in 1:nbFiles)
	ls_slices[[i]] = readRDS(paste0("./competition/", ls_files[i]))

database = rbindlist(ls_slices)

saveRDS(database, "./tree_sStar.rds")
