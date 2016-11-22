rm(list = ls())
# Prepare saved MCMCglmm models for GitHub - reduce the file size
## Take model files saved normally and compress while re-saving


# XXX EDIT THESE TWO ACCORDINGLY
setwd("< EDIT >")
outpath <- "< EDIT >"


# 6.4.2.2 Fixed explicit genetic group effects with Q (from `nadiv`)
load(file = "ggRegMC.RData")
origSize1 <- file.size("ggRegMC.RData")

save("ggRegMC", file = paste0(outpath, "ggRegMC.RData"), compress = "xz", compression_level = 9)
newSize1 <- file.size(paste0(outpath, "ggRegMC.RData"))

percReduction1 <- round(((newSize1 - origSize1) / origSize1)*100, 1)


# 6.4.2.3 Random implicit genetic group effects with A* (from `nadiv`)
load(file = "ggAstarRandMC.RData")
origSize2 <- file.size("ggAstarRandMC.RData")

save("ggAstarRandMC", file = paste0(outpath, "ggAstarRandMC.RData"), compress = "xz", compression_level = 9)
newSize2 <- file.size(paste0(outpath, "ggAstarRandMC.RData"))

percReduction2 <- round(((newSize2 - origSize2) / origSize2)*100, 1)

##############
percReduction1
percReduction2

