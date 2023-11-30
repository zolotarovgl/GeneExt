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
- [x] add proper tempdir naming - don't re-run bam subsampling if a temporary file is found!
- [ ] solve the cases of gene overlaps post-extension (when extension of both genes results in the overlap)
- [ ] fix missing mitochondrial genes   
- [ ] fix case where there are multiple transcripts per gene without gene id => select the longest one  
- [ ] clear logging of the 5'-clipping events - affected 5' and 3' genes  
- [ ] fix coverage estimation - should be more precise  
- [ ] Tadh subsampling bug
- [ ] add .bed to gtf / gff conversion 



   


Input annotation fixes:  
- [x] add adding missing transcript features if missing   

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
- [x] replace transcript   
- [x] replace exon  

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

# Most pressing 
- [x] introduce intermediate files pickup - the bam shouldn't be split every time 
- [x] the format should be guessed from the output1 name and, if none, should be inferred from the input file name 
- [x] the output should be compatible with STARsolo, alevin-fry, cellranger
- [x] replace --keep with --purge for intermediate file keeping 
- [x] explain better what a --mean-coverage means and make it default 
- [x] peplace peakp with --peak_perc
- [x] make --estimate default and remove the option  
- [x] add force option - rerunning everything each time
- [x] introduce picking the longest isoform per gene  
- [x] add reporting missing transcripts per gene  
- [x] add message report for vanilla execution: number of genes extended etc. 
- [x] fix when no extensions are reported
- [x] update the manual  
- [x] add intergenic reads estimation for vanilla execution  
- [x] __MAJOR BUG: missing gene entries after extension__
- [x] __MAJOR BUG: -strand lastexon written as last for gtf file__     
- [x] fix the temporary directory name 
- [x] __weird bug when ran from different location?__ - should use proper conda env
- [x] __MAJOR BUG: kept peaks used for extension__ - can't reproduce so far  
- [ ] problem with gene-to-transcript assignment in gtf file  
- [ ] change orphan peak ids / names to somethehing meaningful  
- [ ] replace test data with sometihng not from octopus
- [ ] __intergenic mapping estimate reporting bug__  
- [ ] find a way to log everything into a log file  
- [ ] make sure re-running is working properly - fixed genome should be picked up  
- [ ] genome correction module should be a separate part of the helper and the main script   
- [ ] selecting the longest transcript doesn't work 
- [ ] gff -> gtf 
- [ ] gff3 extension should be considered as gff otherwise it's a problem
- [ ] check for dependencies in the preflight checks  
- [ ] re-organize manual to go through the steps described in the paper 


