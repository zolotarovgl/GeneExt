# TODOs:
- [x] coverage filtering for the peaks  
  - how to filter the peaks? - so far, by coverage  
- [x] solve the problem with non-unique orphan peak naming  
- [x] proper input/output parsing   
- [x] peak coverage distributions check  
   - results are not clear
- [x] add description of the output files    
- [x] default max size - gene median length  
- [x] speedup coverage computation? 
- [x] .bam file subsampling to a manageable amount of reads 
- [x] verbosity 1: simple output; verbosity 2 : a lot of output  
- [x] filering based on the GENE-OVERLAPPING peaks!
- [x] filtering by the mean coverage     
- [x] Add reporting:   
  - [x] distance from the closest peaks  
- [x] update the function guessing the extension     
- [x] remove big temporary files nog `--keep`  
- [x] input gff -> output gtf  
- [x] solve `bedtools coverage` RAM problem for large datasets: replace with `pysam`   
- [x] check whether read fraction for subsampling works properly - it doesn't: mapped reads vs all of the reads  
- [x] mapping statistics estimation    
- [x] make mapping stats consistent  
- [x] add the function that adds the genes if not found in the file 
  - [x] test with the genomes available in the lab.   
- [x] orphan peaks should be outersected with genic regions   
- [x] __Fix extension to the overlapping peaks__  
  - [x] Distance to the END of the peak  
- [x] add the number of peaks reporting   
- [x] speedup gtf fixing function  
- [x] add gene number reporting in gene adding function  
- [x] add mRNA-> transcript conversion as a part of file fixing 
- [x] __add orphan peak mapping rate estimation separately (otherwise, doesn't make much sense)__
- [x] fork to sebelab
- [x] disable checks if `--estimate` is set
- [x] __gene overextension on the same strand__
- [x] non-extended genes are not written anymore
- [x] update mRNA ranges as well.   
- [x] Solve the same for `.gff` inputcat  
- [x] solve missing CDS features 
- [x] clipped extensions should be written down
- [x] exons from an artificial transcript should have their transcriptIDs updated as well 
- [x] remove peaks overlapping more than one gene  
- [x] __coverage parallelization__  
- [x] __clip antisense extensions if extended into another gene__  - maybe by introducing 2 exon clipping modes?
- [x] fix single-gene contigs removal bug  
- [x] if `--onlyfix` is set, only try fixing the annotation 
- [x] fix empty gtf/gff attributes during parsing 
- [x] gene overextension bug fix - minus strand    
- [x] retain only particular gene features - avoid double quotes     
- [x] add 5' clipping logs - it should be evident which genes have been omitted from the analysis - added
- [x] fix and error with read subsampling  
- [ ] solve the cases of gene overlaps post-extension (when extension of both genes results in the overlap)
- [ ] remove format guessing from the name - it's very confusing for the file parsing  
- [ ] fix missing mitochondrial genes   
- [ ] fix case where there are multiple transcripts per gene without gene id => select the longest one  
- [ ] clear logging of the 5'-clipping events - affected 5' and 3' genes  
- [ ] fix coverage estimation - should be more precise  


   


Input annotation fixes:  
- [ ] add adding missing transcript features if missing   

Minimal cellranger output:  
- [ ] to output `crgtf` files for bed inputs  - gene, transcript, exon
- [ ] parallelize the coverage computation step    

Reporting:
- [ ] fix report path error when called outside of the directory   
- [ ] no report if there is no alignment     
- [ ] add log file  
- [ ] clean all unnecessary files 
- [ ] more explicit intermediate file naming   

Extension modes  
- [ ] replace transcript   
- [ ] replace exon  

Performance:
- [x] 5' gene clipping     
- [ ] try out `gffread` standardized output files, make sure it's compatible (can be accepted by genext)     
- [ ] add extension by continuous coverage   

Orphan peaks:  
- [x] `helper.add_orphan` should be split into getting the orphan peaks and adding them to allow for __peak merging__ if requested   
- [x] make sure cellranger accepts the file with orphan peaks  
- [ ] orphan peak linkage:
  - [x] simple distance merging     
  - [ ] splice junctions if a splice junctions file is provided  
- [x] rename orphan peak clusters  
- [ ] add proper orphan peak linkage - e.g. into the transcripts not full ranges:
  - [ ] linked orphan peaks should represent exons belonging to the transcript cluster.
- [x] __Bug: prevent peaks from clustering across the genes - only allow orphan peak clustering in the intergenic regions__
- [ ] Rename the peaks: UTOs -> Unknown Transcribed Objects, UTS Unknown/Unidentified Transcribed Segments

Manual: 
- [x] add peak filtering manual        
- [x] add orphan peaks manual   
- [ ] add the section covering gene overlaps  
- [ ] __change read direction in the figures__   
- [ ] __add 5'-clipping figure__ 