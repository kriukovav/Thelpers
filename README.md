# Repertoire-based mapping and time-tracking of helper T cell subsets in scRNA-Seq

## This is code to reproduce figures from the TCR-track manuscript

To map *your_query_seurat_object* to our reference scRNA-Seq dataset of helper T cells with Seurat, please follow these instructions:
1. download the reference (https://figshare.com/articles/dataset/full_reference_return_model_rds/23790393?file=45860319)
2. Run
   
```
intgr <- readRDS(path_to_reference dataset)

anchors <- FindTransferAnchors(reference = intgr, query = your_query_seurat_object, reference.reduction = "pca", k.anchor = 5)

query_scRNA2 <- MapQuery(anchorset = anchors, reference = intgr, query = your_query_seurat_object, refdata = "Subset", reference.reduction = "pca", reduction.model = "umap")
```
