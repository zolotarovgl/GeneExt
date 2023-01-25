# GeneExt - Gene extension for improved scRNA-seq data counting   

# TODOs:
* mapping module  
* orphan peaks 
* clean temporary directory  
* make it possible to output gtf files for bed inputs  

* __cellranger mock gtf__
* __orphan peaks__  
* metaplot - coverage relative to TES
* add a proper tag to the orphan peaks source or ditch tagging alltogether 

Try filtering orphan peaks by coverage or something else.



# Table of Contents
1. [Installation](#installation)
2. [Manual](#Manual)


# Installation  

Dependencies:  
* `macs2`  - not needed if you already have the peaks  
* `bedtools`
* `samtools`
* python `gffutils`  
* python `numpy`  

These dependencies can be installed with `conda`: 

```
# create environment
conda create -n geneext
conda activate genext
# install dependencies
conda install -c bioconda -c conda-forge gffutils bedtools numpy macs2 samtools 
```

# Manual  

## Motivation   

Non-model species often have incomplete annotations of their 3-prime untranslated regions (3'-UTRs). At the same time, some of the most popular single-cell RNA sequencing methods are biased towards 3' ends of mRNA molecules. In result, this creates a bias in gene counting for genes with missing 3'-UTRs:   
!['Gene_counting'](./img/gene_counting.png)

To mitigate this effect, `GeneExt` will try to extend the genes in your reference genome by __using the scRNA-seq dataset itself__ (or any 3'-biased type of transcriptomics data). This approach should increase the number of UMI counts registered per gene.  

## How it works  
1. `GeneExt` is supposed to accept alignment files of reads from any 3'-end biased single-cell or bulk RNA-seq protocol. It will then call peaks from this data using `macs2` software and will try to extend genes to the peaks downstream. Alternatively, if you already have peaks you want to extend the genes to (e.g. somehow determined mRNA cleavage site coordinates), having `macs2` installed is not necessary. 
2. For every gene, the most downstream peak will be chosen as a new mRNA cleavage site. The maximal distance from a gene to a peak is controlled by an '-m' parameter (see below).  
3. After genes are extended, `GeneExt` will write an a file which can be used to build a genome reference (e.g. with `cellranger mkref`).  

## Important parameters   

### Maximum extension length   

`-m` parameter specifies the maximum distance the gene is allowed to be extended for. Setting `-m` to larger values will almost always result in longer extensions of genes and thus more reads counted per gene.  
However, the genome annotation will surely __contain unannotated genes__. In such cases, you may actually __misassign the reads to the gene they don't belong to__:  
!['Missing_gene'](./img/missing_gene.png)

Thus, instead of setting `-m` to unrealistically big values, we advice setting it to something biologically meaningful (e.g 1x-2x of median length of gene) and to use it along with calling "orphan peaks" ( `--orphan` option, vis "What are the orphan peaks?").

## Input & Output   

`GeneExt` accepts following input formats and converts them to the output:   
* `bed` &rarr; `crgtf`   
* `gff` &rarr; `crgtf`,`gtf`  
* `gtf` &rarr; `crgtf`,`gtf`  

Output:  

In general, `GeneExt` will try to output a properly formatted `gtf` file that can be used as an input to `cellranger mkref`. However, since `gtf` files vary in their attributes, this may not always be possible ( see "Input debugging").  
For such cases, it is also possible to output a "mock" cellranger gtf file (`crgtf`) with only gene ranges labeled as exons. This file can also be accepted by `cellranger`.  

Note:  
1. In "mock" `gtf` file, every gene will be present as a single feature of a type "exon". This format disregards exon/intron structure of the genes which makes it unusable for analyses which depend on this structure (e.g. RNA-velocity). 
2. If genes are provided in a `bed` file, then the output will always be the mock file.  


## Extension modes  

Depending on the application (or rather your taste), you may want to add extensions differently.  
Currently, the following modes are available:   
* `new_transcript`  
* `new_exon`  

![Extension modes](./img/extension_modes.png)

Note: in general, it doesn't matter which type of extension you choose - it will not affect read counting.  

## What are the "orphan peaks"?  

[t.b.a.]

# Tutorial  

Try executing the tool with the sample data:  

```python geneext.py -g sample/genome_sample.gtf -b sample/sample_10x.bam -m 10000 -o genome_extended.gtf```

Note, here `-m` option specifies the maximal distance to which a gene is allowed to be extended (10kb in this case).   
If the program finishes without an aardvark showing up, you will find the resulting annotation `genome_extended.gtf` along with a temporary directory (`tmp/` by default) which stores intermediate files useful for debugging.  

Now, you can also try running the tool with an `--orphan` setting:

```python geneext.py -g sample/genome_sample.gtf -b sample/sample_10x.bam -m 10000 --orphan -o genome_extended_orphan.gtf```





[t.b.a.]





## Available options   


```man
usage: geneext.py [-h] [-b B] [-g G] [-p P] [-o O] [-m M] [-inf INF]
                  [-ouf OUF] [-t T] [-v V] [-j J] [-e E]

Extend genes in 3' direction using single-cell RNA-seq data

optional arguments:
  -h, --help  show this help message and exit
  -b B        Input .bam file.
  -g G        Genome .gtf/.gff/.bed file.
  -p P        Peaks .bed file. Incompatible with -b.
              If provided, extension is performed using specified peaks coordinates.
              Can be seful in cases of FLAM-seq / Nano-3P-seq data or when manual filtering of the peaks is needed.
  -o O        Output annotation.
  -m M        Maximal distance for gene extension. [10000]
  -inf INF    Input genes file format, if None, will be guessed from a file extension.
  -ouf OUF    Output file format, if not given, will be guessed from a file extension.
  -t T        Temporary directory. [tmp]
  -v V        Verbosity level. 0,[1],2
  -j J        Number of cores for samtools. [1]
  -e E        How to extend the gene (only for .gff/.gtf files) [new_mrna]
                * new_transcript - creates a new transcript feature with the last exon extended
                * new exon - creates an extended last exon
```