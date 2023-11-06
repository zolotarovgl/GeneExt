!['header'](./img/logo.png)

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
__CAVE:__ Please, first test whether the tool works properly on test data! 

For details and troubleshooting, please refer to the manual. If problems persist, don't hesitate to contact the authors.