
#### Aim of prog: List the files to know which results was not saved
#

#### Load libraries
library(stringi)

#### Common variables
nb_arrays = 5000
ls_arrays = 1:nb_arrays

#### Check missing files
## List files
ls_files = list.files(path = "competition", pattern = ".rds")

## Extract number (array_id)
arrays_done = sort(as.integer(stri_sub(str = ls_files, from = 1,
	to = stri_locate(str = ls_files, regex = ".rds")[, 1] - 1)))

## Missing files
(missingFiles = ls_arrays[!(ls_arrays %in% arrays_done)])
(length(missingFiles))

if (length(missingFiles) == 0)
	print("There is no missing file")

#### Write in a file, to copy in the bash script
write(missingFiles, file = "missingFiles.txt",
	ncolumns = length(missingFiles),
    append = FALSE, sep = ",")
