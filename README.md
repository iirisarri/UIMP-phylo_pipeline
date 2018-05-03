# Superbrief Phylogenomics 3h-CrashCourse

This is a very short hands-on course to introduce you to the awesomness of phylogenetic inference.

Nowadays we can infer the evolutionary history from large genomic datasets: phylogenomics.

We will learn a few basic concepts to get you started inferring trees. Note that this does not want to be a comprehensive account of available methods, models, pipelines, etc. The phylogenetic literature is extensive and one can find very good information elsewhere (see a few suggested sources at the very end)

# Data in phylogenomics

<<Data types.

Assembling (good-quality) phylogenomic datasets can be extremely painful and time-consuming. A large part of phylogenomics is spent in mining databases for new genomic data and an even larger part should be used in creating good phylogenomic datasets.

Luckily, we will skip this pain in the ass. 

# Homology, orthology and paralogy

<<explain

# A phylogenomic pipeline

As a toy example, we will use a dataset of several species of a cool marine snail (Conus), composed of nuclear and mitochondrial genes.
Conus are extremely [beautiful and cool animals](https://www.youtube.com/watch?v=zcBmMPJrrKk).

## Obtaining data


Go to https://datasciencehub.ifca.es/ and download the dataset

The dataset is available at: https://github.com/iirisarri/UIMP-phylo_pipeline/Conus_mito_nuclear.zip


## Multiple sequence alignment

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

## Bayesian tree inference

```
mb FcC_supermatrix.nex
```

## Maximum likelihood



# Some notes on systematic error

