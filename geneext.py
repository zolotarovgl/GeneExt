######################################################################################################## 
import argparse
from argparse import RawTextHelpFormatter

import os
import subprocess
import re


parser = argparse.ArgumentParser(description="Extend genes in 3' direction using single-cell RNA-seq data",formatter_class=RawTextHelpFormatter)
parser.add_argument('-b', default= None,help = 'Input .bam file.')
parser.add_argument('-g', default= None,help = 'Genome .gtf/.gff/.bed file.' ,required = True)
parser.add_argument('-p', default= None,help = 'Peaks .bed file. Incompatible with -b.\nIf provided, extension is performed using specified peaks coordinates.\nCan be seful in cases of FLAM-seq / Nano-3P-seq data or when manual filtering of the peaks is needed.')
parser.add_argument('-o', default = None, help = 'Output annotation.',required = True)
parser.add_argument('-m', default = None, help = 'Maximal distance for gene extension.\nIf not set, median length of gene (genomic span!) is used.')
parser.add_argument('-inf', default = None, help = 'Input genes file format, if None, will be guessed from a file extension.')
parser.add_argument('-ouf', default = None, help = 'Output file format, if not given, will be guessed from a file extension.')
parser.add_argument('-t', default = str('tmp'), help = 'Temporary directory. [tmp]')
parser.add_argument('-tag', default = str('GE'), help = 'Tag to be added to the fake gene source and IDs so these can be easily identified downstream. [GE]')
parser.add_argument('-v', default = int(1), help = 'Verbosity level. 0,[1],2')
parser.add_argument('-j', default = '1', help = 'Number of cores for samtools. [1]')
parser.add_argument('-e', default = 'new_transcript', help = 'How to extend the gene (only for .gff/.gtf files) [new_mrna]\n\t* new_transcript - creates a new transcript feature with the last exon extended\n\t* new exon - creates an extended last exon')
parser.add_argument('--orphan',action='store_true', help = 'Whether to add orphan peaks')
parser.add_argument('--peakp',default = 25, help = 'Coverage threshold (percentile of macs2 genic peaks coverage). [1-99, 25 by default].\nThis parameter allows to filter out the peaks based on the coverage BEFORE gene extension. All peaks called with macs2 are required to have a coverage AT LEAST as Nth percentile of the peaks falling within genic regions.')
#parser.add_argument('--orphanp',default = 25, help = 'NOT IMPLEMENTED!\nCoverage threshold (percentile of orphan peaks coverage). [0-100, 75 by default].\nThis parameter allows to filter out the orphan peaks based on the coverage (AFTER gene extension).')
parser.add_argument('--subsamplebam',default = None, help = 'If set, will subsample bam to N reads before the peak calling. Useful for large datasets. Bam file should be indexed.\nDefault: None')
parser.add_argument('--report', action='store_true', help = 'Use this option to generate a report.')
#parser.add_argument('--debug', action='store_true', help = 'Maximum verbosity for debugging')
args = parser.parse_args()


########### Arguments ################
######################################
bamfile = args.b
tempdir = args.t
verbose = args.v
peaksfile = args.p
genefile = args.g 
outputfile = args.o
maxdist = args.m
extension_mode = args.e
threads = args.j
tag = args.tag

# peak coverage percentile:
coverage_percentile = args.peakp

scriptloc = os.path.dirname(os.path.realpath(__file__))

def error_print():
    os.system('cat %s/geneext/err.txt' % scriptloc)
    quit()

# pipeline settings:
do_mapping = False
do_macs2 = False 
do_report = args.report # add the option!
do_orphan = args.orphan
do_subsample = args.subsamplebam is not None


#######################################################################
if peaksfile is None and bamfile is None:
    print("Please, specify either alignment [-b] or peaks file [-p]!")
    error_print()

if outputfile is None:
    print('Please, specify the output file [-o]!')
    error_print()

#if maxdist is None:
#    print('Please, specify the maximum length of gene extension [-m]!')
#    error_print()

if bamfile is not None:
    if os.path.isfile(bamfile):
        do_macs2 = True
        print('Alighment file ... OK')
    else:
        print('Specified alignment file does not exist!')
        error_print()

elif peaksfile is not None:
    if os.path.isfile(peaksfile):
        do_macs2 = False
        print('Found a peaks file, skipping macs2 ...')
    else:
        print('Specified peaks file does not exist!\nPlease, specify either a valid .bam file for macs2 or a peaks file.')
        error_print()

if peaksfile is not None and bamfile is not None:
    print('Please, specify either a .bam file with reads [-b] or a peaks file [-p] but not both at the same time!')
    error_print()

if genefile is None:
    print('Missing genome annotation file [-g]!')
    error_print()

elif not os.path.isfile(genefile):
    print('Genome annotation file .... DOES NOT EXIST!')
    error_print()
else:
    print('Genome annotation file .... OK')

# set temporary directory:
if verbose:
    print("Temporary directory: %s" % tempdir)
if not os.path.exists(tempdir):
   os.makedirs(tempdir)
   if verbose:
        print('Directory created: %s' % tempdir)

########### Functions ################
######################################

from geneext import helper


########### Main #####################
######################################
def parse_input_output_formats():
    # guess the format of input annotation 
    infmt = helper.guess_format(genefile)
    outfmt = 'gtf'
    if verbose:
        print('Input: %s, guessed format: %s\nOutput: %s, guessed format: %s' % (genefile,infmt,outputfile,outfmt))
    return(infmt,outfmt)


def run_peakcalling():
    helper.split_strands(bamfile,tempdir,verbose = verbose,threads = threads)
    helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
    helper.run_macs2(tempdir+'/' + 'minus.bam','minus',tempdir,verbose = verbose)
    helper.collect_macs_beds(outdir = tempdir,outfile = rawpeaks,verbose = verbose)

#def outersect_peaks(genefile = None,peaksfile = None,outputbed = None,verbose = None):
#"""Take macs peaks and outersect them with the genome annotation file"""
#    # convert to bed:
# CAVE! Do you need to convert the input file to bed? 
#    #genefile_bed = tempdir + '/' + 'genes.bed'
#    #helper.gxf2bed(infile = genefile,outfile = genefile_bed)
#    #helper.outersect(inputbed_a = genefile_bed,inputbed_b = peaksfile,outputbed=peaksfilt,by_strand = True, verbose = verbose)
#    helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)


def run_orphan(infmt,outfmt):
    # Get the orpan peaks not assigned to any of the gene and add them to the genome
    # aha, you have to do a second round of outersection but this time also regardless of the strand 
    # to be conservativ
    if infmt != 'bed':
        # remove orphan peaks overlapping extended regions:  
        genefile_ext_bed = tempdir + '/' + 'genes_ext.bed'
        orphan_bed = tempdir + '/' + 'orphan.bed'
        helper.gxf2bed(infile = outputfile,outfile = genefile_ext_bed,featuretype = 'gene')
        helper.outersect(inputbed_a = peaksfilt,inputbed_b=genefile_ext_bed,outputbed = orphan_bed,by_strand = False,verbose = verbose)
        
        # Now, add the orphan peaks:
        if infmt == 'gtf' and outfmt == 'gtf':
            helper.add_orphan_peaks(infile = outputfile,peaksbed = orphan_bed,fmt = 'gtf',tmp_outfile = tempdir + '/' + 'orphan_toadd.' + outfmt,tag = tag)
        elif infmt == 'gff' and outfmt == 'gff':
            helper.add_orphan_peaks(infile = outputfile,peaksbed= orphan_bed,fmt = 'gff')
        else:
            print("Don't know how to add orphan peaks!")
    else:
        print("Don't know how to add orphan peaks!")


#     genefile_bed = tempdir + '/' + 'genes.bed'
#     helper.gxf2bed(infile = genefile,outfile = genefile_bed,featuretype = 'gene')
#     helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile_bed,outputbed=peaksfilt,verbose = verbose)
#else:
#    print("Filtering peaks.")
#    helper.outersect(inputbed_a = genefile,inputbed_b = peaksfile,outputbed=peaksfilt,verbose = verbose)


#####################################

if __name__ == "__main__":
    # parse input and output formats - run the pipeline accordingly
    infmt,outfmt = parse_input_output_formats()
    if not infmt in ['bed','gff','gtf']:
        print('Unknown input format!')
        quit()

    # if -m is not set, get a median gene size:
    if not maxdist:
        maxdist = helper.get_median_gene_length(inputfile = genefile,fmt = infmt)
        if verbose:
            print('Maximum allowed extension length is not set, getting median size of the gene - %s bp.' % str(maxdist))

# 0. MAPPING - not implemented     
    if do_mapping:
        raise(NotImplementedError())
# 0.1 BAM SUBSAMPLING  
    if do_subsample:
        print('======== Running subsampling ========')
        subsampled_bam = tempdir + '/subsampled.bam' 
        nsubs = int(args.subsamplebam)
        # check here if it's an integer
        helper.subsample_bam(inputbam = bamfile,outputbam = subsampled_bam,nreads = nsubs)
        # now, replace for downstream:
        bamfile = subsampled_bam
# 1. MACS2
    if do_macs2:
        print('======== Runnig macs2 ========')
        peaksfile = tempdir + '/' + 'allpeaks.bed'
        helper.split_strands(bamfile,tempdir,verbose = verbose,threads = threads)
        helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
        helper.run_macs2(tempdir+'/' + 'minus.bam','minus',tempdir,verbose = verbose)
        helper.collect_macs_beds(outdir = tempdir,outfile = peaksfile,verbose = verbose)

    else:
        print('Skipping macs2. Running gene extension with %s and %s.' % (peaksfile,genefile))

# 2. Remove peaks overlapping genes:
    #peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
    #print('Removing peaks overlapping genes ... => %s' % peaksfilt)
    ##outersect_peaks(genefile = genefile, peaksfile = peaksfile, outputbed = peaksfilt, verbose = verbose)
    #helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)

# 3. If used macs to call the peaks, filter the peaks by the coverage:
# REPLACE WITH A PIPELINE FUNCTION 
    if do_macs2:    
        print('Calculating peak coverage...')
        covfile = tempdir + '/' + 'allpeaks_coverage.bed'
        #tmp_peakfile = tempdir + '/' + 'allpeaks_noov.bed_tmp'
        # copy filtered file with peaks
        #os.system('cp %s %s' % (peaksfilt,peaksfilt + '_tmp'))

        # compute coverage for all the peaks:
        helper.get_coverage(inputbed_a=peaksfile,input_bam = bamfile,outputfile = covfile,verbose = True)
    
        
        # get the peaks overlapping genes:
        print('Getting genic peaks ...')
        genicpeaksfile = tempdir + '/genic_peaks.bed'
        helper.outersect(inputbed_a = covfile, inputbed_b = genefile,outputbed=genicpeaksfile,by_strand = True,verbose = verbose)
        if not coverage_percentile:
            print('Coverage percentile is not set - retaining all the peaks...')
            count_threshold = 0
        else:
            # get coverage percentile from genic peaks:
            count_threshold = helper.get_coverage_percentile(inputfile = genicpeaksfile,percentile = coverage_percentile,verbose = verbose)
            # hell, you also have to replace the column 
            if verbose:
                print('%s-th coverage percentile for %s is %s reads. Filtering out the peaks below this value...' % (coverage_percentile,genicpeaksfile,str(count_threshold)))
        print('Filtering by coverage ...')


        # get peaks not overlapping the genes and filter them by coverage
        peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
        peaksfiltcov = peaksfilt.replace('.bed','_fcov.bed')
        print('Removing peaks overlapping genes ... => %s' % peaksfilt)
        #outersect_peaks(genefile = genefile, peaksfile = peaksfile, outputbed = peaksfilt, verbose = verbose)
        helper.outersect(inputbed_a = covfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)
        helper.filter_by_coverage(inputfile = peaksfilt,outputfile = peaksfiltcov,threshold = count_threshold,verbose = True)
        quit()

    else:
        # Simply remove peaks overlapping genes:
        peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
        print('Removing peaks overlapping genes ... => %s' % peaksfilt)
        #outersect_peaks(genefile = genefile, peaksfile = peaksfile, outputbed = peaksfilt, verbose = verbose)
        helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)

# 3. Extend genes 
# # REPLACE WITH A PIPELINE FUNCTION 
    helper.extend_genes(genefile = genefile,peaksfile = peaksfilt,outfile = outputfile,maxdist = int(maxdist),temp_dir = tempdir,verbose = verbose,extension_type = extension_mode,infmt = infmt,outfmt = outfmt,tag = tag)
    #helper.extend_genes(genefile,peaksfile,outputfile,int(maxdist),tempdir,verbose,extension_mode,tag = tag)
    if do_orphan:
        run_orphan(infmt = infmt,outfmt = outfmt)
    if do_report:
        helper.do_report()


# but coverage calculation is only possible if one has alignment data 