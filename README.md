!['header'](./img/logo.png)

GeneExt is a tool to improve genome annotations to increase the number of UMIs counted in single-cell RNA-sequencing data (e.g. 10x Chromium).

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
python geneext.py -g [genome .bed/.gtf/.gff] -b [10x.bam] -o [output name]
```
__CAVE:__ Please, first test whether the tool works properly on the test data! 

For details and troubleshooting, please refer to the manual [manual](Manual.md).   
If problems persist, don't hesitate to contact the authors.