
#### Aim of prog: Extract the regression slopes of the fixed effect (12 in total)
#

#### Load packages
library(rstanarm)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Common variables
ls_folders = list.files(path = "./", pattern = "array_[0-9]{1,}")
nbFolders = length(ls_folders)
chosenModel = "model_7.rds"

#### Extract and save parameters
for (folder in ls_folders)
{
	model = readRDS(paste0("./", folder, "/", chosenModel))
	fixef = model$coefficients[1:12]
	saveRDS(fixef, paste0("./", folder, "/fixef.rds"))
}
