
#### Aim: move the data because folders are really messy

library(stringi)
ls_folders = list.dirs(path = ".", recursive = FALSE)

for (folder in ls_folders)
{
	loadPath = paste0(folder, "/data/commondata/data0/")
	ls_files = list.files(path = loadPath)
	if (length(ls_files) != 0)
	{
		ls_rmFiles = list.files(path = folder)[!(list.files(path = folder) %in% c("data", "CONTENTS.txt"))]

		for (file in ls_files)
			file.rename(from = paste0(loadPath, "/", file), to = paste0(folder, "/", file))

		unlink(paste0(folder, "/", ls_rmFiles), recursive = TRUE)
		unlink(paste0(folder, "/data"), recursive = TRUE)
	}
}
