# UIMP-phylo_pipeline

This is a simple 

## Download data

Go to https://datasciencehub.ifca.es/ and download the dataset

The dataset is available at: https://github.com/iirisarri/UIMP-phylo_pipeline/Conus_mito_nuclear.zip



## Concatenate alignment

```
uzip Conus_mito_nuclear.zip
cd Conus_mito_nuclear
FASconCAT-G_v1.02.pl -p -l -s -n -n
```

## Select best-fit evolutionary models for each gene
```
iqtree -s FcC_supermatrix.phy -m TESTONLY -mset mrbayes -nt 1
```

## Bayesian analysis

```
mb FcC_supermatrix.nex
```
