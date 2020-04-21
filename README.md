# Superbrief (3h) Phylogenomics CrashCourse

This is a very-short hands-on course to introduce you to the awesomeness of phylogenetic inference.

Nowadays we can infer the evolutionary history at unprecendented precision using large genomic datasets: phylogenomics.

We will learn the basic concepts to get you started inferring trees. Note that this does not intend to be a comprehensive account of methods, models, software, or pipelines. The phylogenetic literature is extensive and one can find good and more detailed information elsewhere (see a few suggested resources at the end).

A **phylogenetic tree** is a hypothesis of how our data (sequences) evolved. To select the most probable hypothesis we use statistics: probabilistic inference methods. These have explicit assumptions that can be tested and improved.

## Data in phylogenomics

The most commonly used data are DNA or amino acid sequences obtained through high-throughput sequencing. Various types of data types are commonly used in phylogenomics:

- Genomes
- Transcriptomes (RNA-seq)
- Exon capture/ anchored loci / UCEs
- RADseq
- Targeted amplicon sequencing
- Organellar genomes

Besides sequence data, other types of data can be derived from genomes, such as transposable element insertion events and other types of rare genomic changes.

Assembling (good-quality) phylogenomic datasets can be extremely time-consuming and painful. A large part of current phylogenomics practice is about mining databases for new genomic data and an assembling adn curating good phylogenomic datasets.

Luckily, today we can skip this pain in the ass :-)

## Homology, orthology and paralogy

Key concepts to remember when dealing with sequence data.

- **Homology**: origin from common ancestry
- **Orthology**: homology derived from an speciation event
- **Paralogy**: homology derived from a duplication event

## A phylogenomics pipeline

Due to time constraints, we cannot analyze a real genome-scale dataset consisting of hundreds or thousands of genes and species. Instead, we will use a multi-locus mitochondrial and nuclear dataset of several species of cool marine snails (family Conidae). This datset will allow us to practice the basic steps in the phylogenomics pipeline, given that the basic steps are essentially the same with larger datasets. Moreover, cone snails are extremely [beautiful and cool animals](https://www.youtube.com/watch?v=zcBmMPJrrKk).

## Obtaining the data

Let's start by downloading the data.

Connect to [our server](https://datasciencehub.ifca.es/) and clone this repository.
```
git clone https://github.com/iirisarri/UIMP-phylo_pipeline.git
```
Create a new directory `phylosandbox` where we can work with our data without changing the original files. To start, copy the nuclear gene alignments and the new transcriptome you just assembled for *Chelyconus ermineus* (or the Trinity-assembled transcripts in `assemblies/Trinity_ermineus.fasta`).
```
mkdir phylosandbox
cp alignments/nuclear_genes/*fas phylosandbox
cp assemblies/Trinity_ermineus.fasta phylosandbox
cd phylosandbox
```

## Adding new sequences to our alignments

We will add the new assembled transcripts from [*Chelyconus ermineus*](https://en.wikipedia.org/wiki/Conus_ermineus) to the set of nuclear genes present in our dataset.

The easiest way of adding new sequences to a existing dataset is using BLAST. To speed up this process, we will use a custom per script that takes a gene alignment (query) and the sequence collection from the species of interest (the transcriptome) and identifies *homologous* sequences from the database (hits). Then, these hits are retrieved from the database and printed to a new file (.hit.fa). The same could be done using the mitochondrial gene files and the collection of mitochondrial genes assembled for *C. ermineus*.

For the script to run, make sure to have all nuclear gene alignments and the Trinity assembly in the same folder. To avoid problems, simplify the headers of Trinity assemblies. Then, run the script in a for loop.
```
for f in Trinity*.fasta; do sed -E '/>/ s/ len=.+//g' $f > out; mv out $f; done
for f in *fas; do perl ../scripts/blastn_and_extract_hits.pl $f Trinity_ermineus.fasta; done 
```
The `blastn_and_extract_hits.pl` script should have produced one hits.fa file containing a single sequence per gene, which should correspond to the best hit (most similar sequence found in the database). Check that a each file contains only a single best sequence.

NOTE: `blastn_and_extract_hits.pl` requires BioPerl. If this is not installed in your system, you might have gotten some error messages. Do not worry! You can find the gene files with already added *C. ermineus* sequences in `alignments/nuclear_genes_complete/`.

Now, since we want beautiful names in the final alignments, let's change the default Trinity codes. Then, append the new sequences to each of the the existing alignment. Finally, we can clean the folder by removing intermediate files created by the perl script we no longer need.
```
for f in *hits.fa; do sed -E '/>/ s/TRINITY.+/Chelyconus_ermineus/g' $f > out; mv out $f; done
for f in *fas; do  cat $f $f.hits.fa > $f.new.fas; done
rm *nogaps *blastn *hits.fa Trinity_ermineus.fasta*
```
**NOTE**: the best BLAST hit is not necessarily ortholog to the query sequences, but can be a paralog. Paralogs should always be removed prior to phylogenetic inference! One simple way could be to build single-gene trees and look for extremely long branches, which are probably paralogs. But this is not a trivial task and can become very time consuming for large datasets. Still, not a good enough reason to not do a proper job!

## Multiple sequence alignment

At this point, we will have one file per gene, containing all previous species plus *C. ermineus* we just added. Now, we will align each gene file separately. Ideally, nucleotides of protein-coding genes should be aligned at the codon level, using tools such as [TranslatorX](http://translatorx.co.uk/). For non-coding genes, nucleotides should be treated individually. For simplicity, we will just do this using [MAFFT](https://mafft.cbrc.jp/alignment/server/):
```
for f in *.new.fas; do mafft $f > $f.mafft; done
```

## Alignment trimming

Some gene regions (e.g., fast-evolving regions) are difficult to align and thus positional homology is unceratin. It is unclear (probably problem-specific) whether trimming badly-aligned regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inferece. However, gently trimming very incomplete positions (e.g. with >80% gaps) reduces dataset length and thus speeds up computation time without significant information loss.

To trim alignment positions we can use [BMGE](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210) but several other software are also available.

Removing alignment positions with > 80% gaps.
```
for f in *mafft; do java -jar ../scripts/BMGE.jar -i $f -t DNA -op $f.trim.phy -h 1 -g 0.8; done
```
Alternatively, one can also remove high-entropy (likely fast-evolving) columns using de default settings (-m DNAPAM100:2 -g 0.2).
```
for f in *mafft; do java -jar ../scripts/BMGE.jar -i $f -t DNA -op $f.trim.phy; done
```
We will generate an output in *phylip* format that will allow us to quickly check that alignments look nice and clean using `less -S` :-)

## Concatenate alignment

Create a super-alignment by concatenating all gene files. We will use [FASconCAT](link), which will read in all \*.fas \*.phy or \*.nex files in the working directory and concatenate them (in a random order) into a super-alignment. Output format can be chosen and additional information is also printed out.

First of all, create a new directory and copy all relevant files (both the mitochondrial and nuclear gene alignments!)
```
mkdir concatenation
cd concatenation
mv ../*trim.phy .
cp ../../alignments/mitochondrial_genes/*fas .
perl ../../scripts/FASconCAT-G_v1.02.pl -p -l -s -n -n
```
Great!! Our dataset is ready to rock!!

## Maximum likelihood

Now we will start inferring actual trees, yay! Let's build a maximum likelihood (frequentist) estimate of the phylogeny. This is usually simpler and faster than Bayesian inference, among other things because we will not need to check convergence.

In this tutorial we will use [IQTREE](http://www.iqtree.org/), which we had used previously to infer best-fit evolutionary models. Running a maximum likelihood analysis with IQTREE is straightforward. We need to provide the input alignment (`-s`), gene partitions or coordinates (`-spp`), the model (`-m`), number of CPUs to use (`-nt`), and additional parameters. In this case, we will select again best-fit models according to the corrected Akaike Information Criterion (`-m TEST -merit AICc`). Compared to MrBayes, IQTREE allows us to use more models and this is preferable as it can improve model fit. In addition to the maximum likelihood tree, we will assess branch support using 1000 pseudoreplicates of ultrafast bootstrapping (`-bb 1000`). 
```
iqtree -s FcC_supermatrix.fas -spp FcC_supermatrix_partition.txt -m TEST -merit BIC -bb 1000 -nt 1 -pre partitioned
```

Congratulations! If everything went well, here you have the [maximum likelihood estimation of your phylogeny](https://www.youtube.com/watch?v=1FkhCQl2hRs&t=76s) (.treefile)! This can be visualized with FigTree. The numbers at branches (label) are ultrafast bootstrap proportions, which analogously to Bayesian posterior probabilites, inform us about the reliability of that branch. Values >70% can be trusted as robust.

## Select best-fit evolutionary models for each gene

Evolutionary models are mathematical descriptions of the process by which our data evolved. We need to assume a explicit model in probabilistic phylogenetic inference. 

Several evolutonary models are available. Which one to choose? The model that best fits our data is likely the one under which our data evolved. We will choose the best-fit model using a statistic criterion (Bayesian Information Content, BIC).

Let's select best-fit models for each gene separately. Among all possible models, we will first compare only those implemented in MrBayes, as it is our first analysis. This can be easily done with the [IQTREE](http://www.iqtree.org/) software:
```
iqtree -s FcC_supermatrix.fas -spp FcC_supermatrix_partition.txt -m TESTONLY -mset mrbayes -nt 1 -pre modelfinder
```

The best-fit models will be printed to screen (also available in the `.log` and `.best_scheme.nex` files).


## Bayesian tree inference

#### Setting the analysis

Now we will infer a Bayesian phylogenetic tree.

In [MrBayes](http://mrbayes.sourceforge.net/), the input file contains both the alignment and all necessary commands (at the end of the file). To prepare the input file, we will need to manually edit the *nexus* concatenated file (`FcC_supermatrix_partition.nex`, NOT `FcC_supermatrix.nex`). First, we need to specify the best-fit models inferred in the previous analysis (look for `lset` by the end of the file).

How to set up the models?
- HKY and K2P correspond to nst=2
- GTR and SYM correspond to nst=6
- Rate heterogeneity among sites is specified with rates=gamma (+G), rates=inv (+I) or rates=invgamma (+I+G)
Modify the lset line and add as many lset lines as required to specify all your models.

Then, edit the line referring to MCMC by changing `mcmc` to `mcmcp`. Now the input file is ready. If you had any problems formatting the nexus file, a input file can be found in `mrbayes/Conidae_supermatrix_partition.nex`.

Let's run MrBayes. Setting `mcmcp` and running MrBayes interactively (`-i`) allow us to check that the infile is correct before actually running the program.
```
mb -i FcC_supermatrix_partition.nex
```
If no errors are shown, we can start to run the MCMC chains by typing `mcmc` in the MrBayes prompt. Leave it running in the backgroud (`ctl+z` & `bg`) for a few hours.

#### Checking convergence

After the runs are completed, we must check whether the runs have converged. How do we know if they did? The simplest way is to check the maximum difference of standard deviation of split frequencies printed out by MrBayes, which should be below 0.01. If this is the case, then we are good to go. If not, then we should run the analysis longer, or change the priors.

Bear in mind that convergence is serious bussiness and in a serious world one should also use other tools, such as [AWTY](https://www.ncbi.nlm.nih.gov/pubmed/17766271) or [Tracer](http://tree.bio.ed.ac.uk/software/tracer/). But lest's leave this for next time :-)

#### Summarizing the posterior

Once the runs are converged, we can summarize the posterior estimates from the two MCMC chains to obtain our tree and posterior probabilities.

Before summarizing the MCMC chains, we must exclude the initial generations because these contain bad stimates of the paramters and trees (we say that MCMC chains did not reach stationarity). Usually discarding 10-25% of the initial cycles is enough (known as "burnin"), but the appropriate burnin proportion should be decided using the convergence tools mentioned above.

To summarize the MCMC chains in the stationary phase, we need to initialize again MrBayes. Make sure your nexus file has mcmcp and not mcmc, since running MrBayes with mcmc will overwrite the runs that just finished! Then, we summarize the parameters (`sump`) and trees (`sumt`) indicating the proportion of generations that we want to exclude as burnin (`relburnin`). If your analysis did not finish, a finished analysis can be found in the folder `mrbayes`.

```
mb -i FcC_supermatrix_partition.nex
 > sump relburnin 0.1
 > sumt relburnin 0.1
```
Congrats! You just generated your first [Bayesian estimate of phylogeny](https://www.youtube.com/watch?v=RMNwsdb5VU4)!

Download the consensus file (.con), which is a summary of the trees sampled by the two MCMC chains during the stationary phase. The tree can be visualized with [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or any other visualizing software. The numbers at branches (label) correspond to posterior probabilities, which inform us about the reliability of each branch in the phylogeny. Posterior probabilities >0.95 can be trusted as robust.

# Some notes on phylogenomics and systematic error

The use of genome-scale data drastically reduces stochastic error and thus produce very precise phylogeny inferences. However, phylogenomics might be more sensitive to violations of model assumptions. Such violations might be negligible in small-scale phylogenetics, but can be exacerbated when using hundreds or thousands of genes. This is known as **systematic error**, and in phylogenomic practice one *should* always ensure that the obtained phylogenies are robust to known sources of systematic biases. Systematic biases typically refer to evolutionary processes that are heterogeneous across the genome, between species, and in time. Often these reflect biological adaptations, such as the independent evolution of GC-rich genomes in heat-adapted bacteria.

Some of the above-mentioned heterogeneities can be accounted for with more sophisticated models than those used in this tutorial. Among the most commonly used ones are **mixture models** that account for the heterogeneity of evolution between genes and amino acids. These models typically fit the data better and thus produce more accurate estimates of phylogeny. Since these models are more complex, they require many more parameters to be estimated from the data and they are computationally more expensive.

With increasing sizes of genomic datasets and the the use of more and more complex models (to account for systematic error) the best is to find a trade-off between alignment size and the sophistication of the model.

Because large datasets increase the precision of phylogeny, commonly used **branch support** values (non-parametric bootsrapping and even more Bayesian posterior probabilities) are usually **inflated**. Therefore, high support from such measures in phylogenomics is only partly informative. More stringent measurements are being developed, such as [gene jackknifing](https://www.nature.com/articles/s41559-017-0240-5) or [quartet-based measures](https://onlinelibrary.wiley.com/doi/abs/10.1002/ajb2.1016). The classic on-parametric bootstrapping also requires a lot of computational effort for large datasets, and thus other faster approximate measures are being developed, including the ultra-fast bootstrapping used above.

Finally, a major realization of the genomic era has been that phylogenies reconstructed from different genes or genomic regions differ from each other and/or from the global **species tree**. The reasons can be biological or methodological. Methods that reconcile these incongruences are also commonly used, e.g. [ASTRAL](https://github.com/smirarab/ASTRAL).

# Useful resources / Further reading

- **Tree thinking – essential in bioinformatics**. [The tree-thinking challenge](The tree-thinking challenge) by D. Baum et al. Science (2005) and [Tree thinking for all biology: the problem with reading phylogenies as ladders of progress](https://onlinelibrary.wiley.com/doi/full/10.1002/bies.20794) by KE Omland et al. BioEssays (2008).
- Inferencia filogenética, Chapter 9 in [Bioinformática con Ñ](https://www.scribd.com/doc/231270078/Bioinformatica-con-N) by F. Abascal, I. Irisarri and R. Zardoya.
- [High-throughput genomic data in systematics and phylogenetics](https://www.annualreviews.org/doi/10.1146/annurev-ecolsys-110512-135822) by EM Lemmon and AR Lemmon AREES (2013)
- [Phylogenomics – An Introduction](https://www.springer.com/gb/book/9783319540627) by C. Bleidorn (2017).
- [The Phylogenetic Handbook 2nd Ed.](http://www.cambridge.org/us/academic/subjects/life-sciences/genomics-bioinformatics-and-systems-biology/phylogenetic-handbook-practical-approach-phylogenetic-analysis-and-hypothesis-testing-2nd-edition?format=PB&isbn=9780521730716#mzxhCgA9pXliYy7j.97) edited by P. Lemey et al. (2009)
