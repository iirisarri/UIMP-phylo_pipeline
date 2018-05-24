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

Key concepts to remember when dealing with sequence data.

**Homology**: origin from common ancestry
**Orthology**: homology derived from an speciation event
**Paralogy**: homology derived from a duplication event

## A phylogenomics pipeline

As a toy example, we will use a dataset of several species of a cool marine snail (Conidae), composed of nuclear and mitochondrial genes.
Cone snails are extremely [beautiful and cool animals](https://www.youtube.com/watch?v=zcBmMPJrrKk).

## Obtaining data

Let's start by downloading the data.

Connect to [our server](https://datasciencehub.ifca.es/).

You can download the dataset from [this respository](https://github.com/iirisarri/UIMP-phylo_pipeline/conidae_mito_nuclear.zip)
Put the data in your preferred location and decompress.
```
wget https://github.com/iirisarri/UIMP-phylo_pipeline/conidae_mito_nuclear.tar.gz
tar zxvf conidae_mito_nuclear.tar.gz
```
Create a new directory where we can work with our data without changing the original files. To start, copy the nuclear gene alignments and the new transcriptome you just assembled (or Trinity-assembled transcripts in `assemblies/Trinity_ermineus.fasta`).
```
cd conidae_mito_nuclear
mkdir phylosandbox
cp alignments/nuclear_genes/*fas phylosandbox
cp assemblies/Trinity_ermineus.fasta phylosandbox
cd phylosandbox
```

## Adding new sequences to our alignments

The easiest way of adding new sequences to a existing dataset is using BLAST. We will add the assembled transcripts of [*Conus ermineus*](https://en.wikipedia.org/wiki/Conus_ermineus) to the nuclear genes. To speed up this process, we will use a custom per script that takes a gene alignment \(query\) and the collection of sequences from the species of interest, in this case the transcriptome \(database\). It blasts each sequence in the query against each sequence in the database to find the best hit, then this hit is retrieved from the database and printed to a file (.hit.fa). The same could be done using the mitochondrial gene files and the collection of mitochondrial genes assembled for *C. ermineus*.

For the script to run, make sure to have all nuclear gene alignments and the Trinity assembly in the same folder. To avoid problems, simplify the headers of Trinity assemblies. Then, run the script in a for loop.
```
for f in Trinity*.fasta; do sed -E '/>/ s/ len=.+//g' $f > out; mv out $f; done
for f in *fas; do perl ../scripts/blastn_and_extract_hits.pl $f Trinity_ermineus.fasta; done 
```
The blastn_and_extract_hits.pl script should have produced one hits.fa file containing a single sequence, which should correspond to the best hit (most similar sequence found in the database). Check that a each file contains only a single best sequence.
Now, since we want beautiful names in the final alignments, let's change the default Trinity codes. Then we will be ready to append the new sequence to the existing alignment. Finally, we can clean the folder by removing intermediate files created by the perl script.
```
for f in *hits.fa; do sed -E '/>/ s/TRINITY.+/Chelyconus_ermineus/g' $f > out; mv out $f; done
for f in *fas; do  cat $f $f.hits.fa > $f.new.fas; done
rm *nogaps *blastn *hits.fa Trinity_ermineus.fasta*
```
NOTE: the best BLAST hit is not necessarily ortholog to the query sequences, but can be a paralog. Paralogs should always be removed prior to phylogenetic inference! One simple way could be to build single-gene trees and look for extremely long branches, which are probably paralgs. But this is not a trivial task and can become very time consuming for large datasets. Still, not a good enough reason to not do a proper job!

## Multiple sequence alignment

We will align each gene file separately. Ideally, nucleotides of protein-coding genes should be aligned at the codon level, using tools such as [TranslatorX](http://translatorx.co.uk/). For non-coding genes, nucleotides should be treated individually, as done by mafft:
```
for f in *.new.fas; do mafft $f > $f.mafft; done
```

## Alignment trimming

Some gene regions (e.g., fast-evolving regions) are difficult to align and thus positional homology is unceratin. It is unclear (probably problem-specific) whether trimming badly-aligned regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inferece. However, gently trimming very incomplete positions (e.g. with >80% gaps) reduces dataset time and thus speeds up computation time. 

Removing alignment positions with > 80% gaps.
```
for f in *mafft; do java -jar ../scripts/BMGE.jar -i $f -t DNA -op $f.trim.phy -h 1 -g 0.8; done
```
Alternatively, one can also remove high-entropy columns using de default settings (-m DNAPAM100:2 -g 0.2).
```
for f in *mafft; do java -jar ../scripts/BMGE.jar -i $f -t DNA -op $f.trim.phy; done
```
We will generate an output in *phylip* format that will allow us to quickly check that alignments look nice and clean using `less -S` :-)

## Concatenate alignment

Create a super-alignment by concatenating all gene files. We will use [FASconCAT](link), which will read in all \*.fas \*.phy or \*.nex files in the working directory and concatenate them (randomly) into a super alignment. Output format can be chosen and additional information is also printed out.

First of all, create a new directory and copy all relevant files (both the mitochondrial and nuclear gene alignments!)
```
mkdir concatenation
cd concatenation
mv ../*trim.phy .
cp ../../alignments/mitochondrial_genes/*fas .
perl ../../scripts/FASconCAT-G_v1.02.pl -p -l -s -n -n
```

## Select best-fit evolutionary models for each gene

Evolutionary models are mathematical descriptions of the process by which our data evolved. We need to assume a explicit model in probabilistic phylogenetic inference. 

Several evolutonary models are available. Which one to choose? The one that fits best our data is likely the one under which our data was generated. We will select the model that best fits our data.

In our case, we will select best-fit models for each gene separately. Among all possible models, we will first compare only those implemented in MrBayes, as it is our first analysis. This can be easily done with the iqtree software:
```
iqtree -s FcC_supermatrix.phy -spp FcC_supermatrix_partition.txt -m TESTONLY -mset mrbayes -nt 1
```

The best-fit models will be printed to screen (also available in the .log and .best_scheme.nex files).

## Bayesian tree inference
### Setting the analysis
Now we will start inferring actual trees, yay! We will start with Bayesian inference, because it takes longer to compute.

In MrBayes, the input file contains both the alignment and all necessary commands (at the end of the file). To prepare the input file, we will need to manually edit the *nexus* concatenated file (FcC_supermatrix_partition.nex, NOT FcC_supermatrix.nex). We will need to update the best-fit models based on our previous results (look for lset).

How to set up the models?
HKY and K2P correspond to nst=2
GTR and SYM correspond to nst=6
rate heterogeneity among sites is specified with rates=gamma (+G), rates=inv (+I) or rates=invgamma (+I+G)
Modify the lset line and add as many lset lines as required to specify all your models.

**Change number of generations?

Then, edit the line referring to MCMC by changing mcmc by mcmcp.

Once the input file is ready, let's simply run MrBayes. Good luck!!
```
mb FcC_supermatrix.nex
```
Setting mcmcp allow us to check whether the infile was correct. If so, we can start to run the MCMC chains by typing mcmc in the prompt.
### Checking convergence

After the runs are completed, we must check whether they converged. If not, then we should run the analysis longer or change the priors.

How do we know whether the runs have converged? A very simple way is to the the standard deviation of split frequencies, which should be below 0.01. If this is the case, then we are good to go.

Bear in mind that convergence is serious bussiness and one should normally use other tools, such as [AWTY](https://www.ncbi.nlm.nih.gov/pubmed/17766271) or [Tracer](http://tree.bio.ed.ac.uk/software/tracer/). But lest's leave this for next time.

### Summarizing the posterior

One the runs are converged, we can summarize the two MCMC chains to obtain our tree and posterior probabilities.

The begining of the MCMC chains is not very good (we say the sampled before reaching stationarity, i.e. "bad trees"). 

## Maximum likelihood



# Some notes on systematic error

