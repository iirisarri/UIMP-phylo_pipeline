# Superbrief Phylogenomics 3h-CrashCourse

This is a very short hands-on course to introduce you to the awesomness of phylogenetic inference.

Nowadays we can infer the evolutionary history from large genomic datasets: phylogenomics.

We will learn a few basic concepts to get you started inferring trees. Note that this does not want to be a comprehensive account of available methods, models, pipelines, etc. The phylogenetic literature is extensive and one can find very good information elsewhere (see a few suggested sources at the very end)

A phylogenetic tree is a hypothesis of how our data (sequences) evolved. To select the most probable hypothesis we use statistics: probabilistic inference methods. These have explicit assumptions that can be tested and improved.

## Data in phylogenomics

<<Data types.

Assembling (good-quality) phylogenomic datasets can be extremely painful and time-consuming. A large part of phylogenomics is spent in mining databases for new genomic data and an even larger part should be used in creating good phylogenomic datasets.

Luckily, today we can skip this pain in the ass :-)

## Homology, orthology and paralogy

Homology: origin from common ancestry.
Orthology: homology derived from an speciation event.
Paralogy: homology derived from a duplication event.

## A phylogenomics pipeline

As a toy example, we will use a dataset of several species of a cool marine snail (*Conus*), composed of nuclear and mitochondrial genes.
*Conus* are extremely [beautiful and cool animals](https://www.youtube.com/watch?v=zcBmMPJrrKk).


## Obtaining data

Let's start by downloading the data.

Connect to [our server](https://datasciencehub.ifca.es/).

You can download the dataset from [this respository](https://github.com/iirisarri/UIMP-phylo_pipeline/Conus_mito_nuclear.zip)
Put the data in your preferred location, decompress and go to Alignments/Mito_and_nuclear.

```
wget https://github.com/iirisarri/UIMP-phylo_pipeline/Conus_mito_nuclear.zip
cd Conus_mito_nuclear/Alignments/Mito_and_nuclear
```

## Adding new sequences to our alignments

The easiest way of adding new sequences to a existing dataset is using BLAST. We will add the new transcriptome and mitogenome assembled for [*Conus ermineus*](https://en.wikipedia.org/wiki/Conus_ermineus).

Create the database with the new transriptome and mitogenome
```
makeblastdb -in your_assembly -dbtype nucl -parse_seqids
```
BLAST the nuclear genes against the transcriptome. Repeat it with mitochondrial genes and the newly assembled mitogenome
```
for f in *fas; do blastn -query $f -db my_database -num_descriptions 1 -num_alignments 1 > my_outfile.blastn
```
Extract hits
```
perl
```

## Multiple sequence alignment

We will align each gene file separately.

```
for f in *.fas; do mafft $f > out; mv out $f; done
```

## Concatenate alignment

Create a super-alignment by concatenating all gene files. We will use [FASconCAT](link), which will read in all \*.fas \*.phy or \*.nex files in the working directory and concatenate them (randomly) into a super alignment. Output format can be chosen and additional information is also printed out.

```
FASconCAT-G_v1.02.pl -p -l -s -n -n
```

## Select best-fit evolutionary models for each gene

Evolutionary models are mathematical descriptions of the process by which our data evolved. We need to assume a explicit model in probabilistic phylogenetic inference. 

Several evolutonary models are available. Which one to choose? The one that fits best our data is likely the one under which our data was generated. We will select the model that best fits our data.

In our case, we will select best-fit models for each gene separately. Among all possible models, we will first compare only those implemented in MrBayes, as it is our first analysis. This can be easily done with the iqtree software:

```
iqtree -s FcC_supermatrix.phy -m TESTONLY -mset mrbayes -nt 1
```

## Bayesian tree inference

Now we will start inferring actual trees, yay! We will start with Bayesian inference, because it takes longer to compute.

In MrBayes, all commands are read from the infile that also contains the alignments. So first of all we need to prepare the input file, which is in nexus format. To do this, we need to manually edit

```
mb FcC_supermatrix.nex
```

## Maximum likelihood



# Some notes on systematic error

