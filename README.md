!['header'](./img/logo.png)

GeneExt takes as input scRNA-seq mapped reads and a gene annotation file (GTF or GFF, any version) and outputs an extended gene annotation file for improved scRNA-seq transcript quantification.

# Installation  

Tool dependencies can be installed with `conda` or `mamba`: 

```
# create environment
conda create -n geneext
conda activate geneext
# install dependencies
conda install -c bioconda -c conda-forge gffutils bedtools numpy macs2 samtools pysam rich pandas 
```

# Test run
Once dependencies are installed, try running `GeneExt` with sample data:  
```
python geneext.py -g test_data/annotation.gtf -b test_data/alignments.bam -o result.gtf
```
Please, first test whether the tool works properly on the test data! 

For details on how to obtain alignment file, please refer to the manual [manual](Manual.md).   
If problems persist, don't hesitate to contact the authors.