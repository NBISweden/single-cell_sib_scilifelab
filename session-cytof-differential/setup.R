library(HDCytoData)
fs <- Bodenmiller_BCR_XL_flowSet()

url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
for (fn in c(
    "PBMC8_metadata.xlsx", "PBMC8_panel_v3.xlsx", 
    "PBMC8_cluster_merging1.xlsx", "PBMC8_cluster_merging2.xlsx"))
    download.file(file.path(url, fn), destfile = fn, mode = "wb") 
