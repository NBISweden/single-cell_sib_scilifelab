## Dataset description:
2358 genes x 9288 cells.
Mammary epithelial cells from three independent studies:
spk (PumbedID_30089273): 2615 cells
vis (PubmedID_29158510): 3047 cells
wal (PubmedID_29225342): 3626 cells

## Cell Filtering:
Each study contains only cells that are assigned unambiguously to one of three major
mammary epithelial cell types; luminal_pogenitors luminal mature and basal:

                      spk  vis  wal
  luminal_progenitor 1097  694  729
  basal               475 1257  780
  luminal_mature     1043 1096 2117

## Gene Filtering:  
Datasets have been gene-filtered to only retain 2358 genes that are
A. Among the top 50% variable genes in **each of the 3 studies** AND
B. Present in >5% of the cells in **each of the 3 studies**.
