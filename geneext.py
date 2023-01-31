######################################################################################################## 
import argparse
from argparse import RawTextHelpFormatter

import os
import subprocess
import re


parser = argparse.ArgumentParser(description="Program: GeneExt (Extend genes in 3' direction using single-cell RNA-seq data)\nVersion: 0.7",formatter_class=RawTextHelpFormatter)
parser.add_argument('-g', default= None,help = 'Genome .gtf/.gff/.bed file.' ,required = True) 
parser.add_argument('-b', default= None,help = 'Input .bam file.')
parser.add_argument('-p', default= None,help = 'Peaks .bed file. Incompatible with -b.\nIf provided, extension is performed using specified peaks coordinates.\nCan be seful in cases of FLAM-seq / Nano-3P-seq data or when manual filtering of the peaks is needed.') 
parser.add_argument('-o', default = None, help = 'Output annotation.\n\n\n================ Additional arguments ================\n',required = True)
parser.add_argument('-m', default = None, help = 'Maximal distance for gene extension.\nIf not set, a median length of gene (genomic span!) is used.')
parser.add_argument('-inf', default = None, help = 'Input genes file format, if None, will be guessed from a file extension.')
parser.add_argument('-ouf', default = None, help = 'Output file format, if not given, will be guessed from a file extension.')
parser.add_argument('-t', default = str('tmp'), help = 'Temporary directory. [tmp]')
parser.add_argument('-tag', default = str('GE'), help = 'Tag to be added to the fake gene source and IDs so these can be easily identified downstream. [GE]')
parser.add_argument('-v', default = int(1), help = 'Verbosity level. 0,[1],2')
parser.add_argument('-j', default = '1', help = 'Number of cores for samtools. [1]')
parser.add_argument('-e', default = 'new_transcript', help = 'How to extend the gene (only for .gff/.gtf files) [new_mrna]\n\t* new_transcript - creates a new transcript feature with the last exon extended\n\t* new exon - creates an extended last exon')
parser.add_argument('--orphan',action='store_true', help = 'Whether to add orphan peaks')
parser.add_argument('--mean_coverage', action='store_true', help = 'Whether to use mean coverage for peak filtering.\nMean coverage = [ # mapping reads]/[peak width].')
parser.add_argument('--peakp',default = 25, help = 'Coverage threshold (percentile of macs2 genic peaks coverage). [1-99, 25 by default].\nAll peaks called with macs2 are required to have a coverage AT LEAST as N-th percentile of the peaks falling within genic regions.\nThis parameter allows to filter out the peaks based on the coverage BEFORE gene extension.')
#parser.add_argument('--orphanp',default = 25, help = 'NOT IMPLEMENTED!\nCoverage threshold (percentile of orphan peaks coverage). [0-100, 75 by default].\nThis parameter allows to filter out the orphan peaks based on the coverage (AFTER gene extension).')
parser.add_argument('--subsamplebam',default = None, help = 'If set, will subsample bam to N reads before the peak calling. Useful for large datasets. Bam file should be indexed.\nDefault: None')
parser.add_argument('--report', action='store_true', help = 'Use this option to generate a PDF report.')
parser.add_argument('--keep', action='store_true', help = 'Use this to keep .bam and other temporary files in the a temporary directory. Useful for debugging.')
parser.add_argument('--estimate', action='store_true', help = 'NOT IMPLEMENTED\nWhether to estimate intergenic read proportion.\nUseful for quick checking intergenic mapping rate.')
args = parser.parse_args()

callcmd = 'python ' + os.path.basename(__file__) + ' '+ " ".join(["-"+str(k)+' '+str(v) for k,v in zip([arg for arg in vars(args)],[getattr(args,arg) for arg in vars(args)]) if v ])
print(callcmd)
#print([k+" "+v for k,v in zip([arg for arg in vars(args)],[getattr(args,arg) for arg in vars(args)])  )

########### Arguments ################
######################################
bamfile = args.b
tempdir = args.t
verbose = int(args.v)
peaksfile = args.p
genefile = args.g 
outputfile = args.o
maxdist = args.m
extension_mode = args.e
threads = args.j
tag = args.tag

# peak coverage percentile:
mean_coverage = args.mean_coverage
coverage_percentile = args.peakp

scriptloc = os.path.dirname(os.path.realpath(__file__))

def error_print():
    os.system('cat %s/geneext/err.txt' % scriptloc)
    quit()

# pipeline settings:
do_mapping = False
do_macs2 = False  
do_orphan = args.orphan
do_subsample = args.subsamplebam is not None
do_estimate = args.estimate and bamfile
do_clean = not args.keep
do_orphan_merge = False  



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
if verbose > 0:
    print("Temporary directory: %s" % tempdir)
if not os.path.exists(tempdir):
   os.makedirs(tempdir)
   if verbose > 0:
        print('Directory created: %s' % tempdir)

do_report = args.report and do_macs2

########### Functions ################
######################################

from geneext import helper


########### Main #####################
######################################
def parse_input_output_formats():
    # guess the format of input annotation 
    infmt = helper.guess_format(genefile)
    outfmt = 'gtf'
    if verbose > 0:
        print('Input: %s, guessed format: %s\nOutput: %s, guessed format: %s' % (genefile,infmt,outputfile,outfmt))
    return(infmt,outfmt)


def run_peakcalling():
    helper.split_strands(bamfile,tempdir,verbose = verbose,threads = threads)
    helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
    helper.run_macs2(tempdir+'/' + 'minus.bam','minus',tempdir,verbose = verbose)
    helper.collect_macs_beds(outdir = tempdir,outfile = rawpeaks,verbose = verbose)

def run_orphan(infmt,outfmt,verbose,merge = False):
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
            helper.add_orphan_peaks(infile = outputfile,peaksbed = orphan_bed,fmt = 'gtf',tmp_outfile = tempdir + '/' + 'orphan_toadd.' + outfmt,tag = tag,verbose = verbose)
        elif infmt == 'gff' and outfmt == 'gff':
            helper.add_orphan_peaks(infile = outputfile,peaksbed= orphan_bed,fmt = 'gff',verbose=verbose)
        else:
            print("Don't know how to add orphan peaks!")
    else:
        print("Don't know how to add orphan peaks!")

# Reporting functions 
def generate_report():
    """This function will take as an input gene extensions and plot the distributions"""
    # script.R maxdist quant closest_gene_file allpeaks_cov_file allpeaks_noov_file extension_file
    #args = c('10000','.25','tmp/genes_peaks_closest','allpeaks_coverage.bed','tmp/allpeaks_noov.bed','tmp/extensions.tsv')
    os.system('Rscript geneext/report.r %s %s %s %s %s %s %s' % (maxdist,coverage_percentile/100,tempdir + '/_genes_peaks_closest',covfile,peaksfilt,tempdir+'/extensions.tsv',str(verbose)))

def clean_tmp(tempdir = None):
    # clean temporary directory of big files 
    toremove = [tempdir+'/'+x for x in os.listdir(tempdir) if '.bam' in x or x[0] == '_']
    for file in toremove:
        if verbose > 0:
            print("Removing %s" % file)
        os.remove(file)

#####################################

if __name__ == "__main__":
    print('======== Preflight checks ======================')
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
    if verbose > 0:
        print('Checks done.')

# 0. MAPPING - not implemented     
    if do_mapping:
        print('======== Running mapping =======================')
        raise(NotImplementedError())
# 0.1 BAM SUBSAMPLING  
    if do_subsample:
        print('======== Running subsampling ===================')
        subsampled_bam = tempdir + '/subsampled.bam' 
        nsubs = int(args.subsamplebam)
        if not os.path.isfile(bamfile + '.bai'):
            if verbose > 0:
                print('Indexing %s' % bamfile)
            helper.index_bam(bamfile,verbose = verbose,threads=threads)
        # check here if it's an integer
        helper.subsample_bam(inputbam = bamfile,outputbam = subsampled_bam,nreads = nsubs,verbose = verbose,threads=threads)
        # now, replace for downstream:
        if verbose > 0:
            print('Subsampling done.')
        bamfile = subsampled_bam

# 1. MACS2
    if do_macs2:
        print('======== Running macs2 =========================')
        peaksfile = tempdir + '/' + 'allpeaks.bed'
        helper.split_strands(bamfile,tempdir,verbose = verbose,threads = threads)
        helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
        helper.run_macs2(tempdir+'/' + 'minus.bam','minus',tempdir,verbose = verbose)
        helper.collect_macs_beds(outdir = tempdir,outfile = peaksfile,verbose = verbose)
        if verbose > 0:
            print('Macs2 done.')
    else:
        if verbose > 0:
            print('Skipping macs2. Running gene extension with %s and %s.' % (peaksfile,genefile))

# 3. If used macs to call the peaks, filter the peaks by the coverage:
# REPLACE WITH A PIPELINE FUNCTION 
    if do_macs2:    
        print('======== Filtering macs2 peaks =================')
        covfile = tempdir + '/' + 'allpeaks_coverage.bed'
        # compute coverage for all the peaks:
        # check if bam file is indexed:
        if not os.path.isfile(bamfile + '.bai'):
            if verbose > 0:
                print('Indexing %s' % bamfile)
            helper.index_bam(bamfile,verbose = verbose,threads=threads)
        if verbose > 0:
            print('Computing coverage ...')
        helper.get_coverage(inputbed_a=peaksfile,input_bam = bamfile,outputfile = covfile,verbose = verbose,mean = mean_coverage)
        # get the peaks overlapping genes:
        if verbose > 0:
            print('Getting genic peaks ...')
        genicpeaksfile = tempdir + '/genic_peaks.bed'
        helper.outersect(inputbed_a = covfile, inputbed_b = genefile,outputbed=genicpeaksfile,by_strand = True,verbose = verbose)
        if not coverage_percentile:
            if verbose > 0:
                print('Coverage percentile is not set - retaining all the peaks...')
            count_threshold = 0
        else:
            # get coverage percentile from genic peaks:
            count_threshold = helper.get_coverage_percentile(inputfile = genicpeaksfile,percentile = coverage_percentile,verbose = verbose)
            # hell, you also have to replace the column 
            if verbose > 0:
                print('%s-th %s coverage percentile for %s is %s %s. Filtering out the peaks below this value...' % (coverage_percentile,'mean' if mean_coverage else '' ,genicpeaksfile,str(count_threshold),'reads/base' if mean_coverage else 'reads'))
        # get peaks not overlapping the genes and filter them by coverage
        peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
        peaksfiltcov = peaksfilt.replace('.bed','_fcov.bed')
        if verbose > 0:
            print('Removing peaks overlapping genes ... => %s' % peaksfilt)
        helper.outersect(inputbed_a = covfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)
        helper.filter_by_coverage(inputfile = peaksfilt,outputfile = peaksfiltcov,threshold = count_threshold,verbose = True)
    else:
        # Simply remove peaks overlapping genes:
        peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
        if verbose > 0:
            print('Removing peaks overlapping genes ... => %s' % peaksfilt)
        #outersect_peaks(genefile = genefile, peaksfile = peaksfile, outputbed = peaksfilt, verbose = verbose)
        helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)

# 3. Extend genes 
    print('======== Extending genes =======================')
    helper.extend_genes(genefile = genefile,peaksfile = peaksfilt,outfile = outputfile,maxdist = int(maxdist),temp_dir = tempdir,verbose = verbose,extension_type = extension_mode,infmt = infmt,outfmt = outfmt,tag = tag)
    quit()
    if do_orphan:
        print('======== Adding orphan peaks ===================')
        run_orphan(infmt = infmt,outfmt = outfmt,verbose = verbose,merge = do_orphan_merge)
    if do_report:
        print('======== Creating report =======================')
        generate_report()
        if verbose > 0:
            print('Report: report.pdf')
    if do_estimate:
        print('======== Estimating intergenic mapping =========')
        genicbed = tempdir + '/genic.bed'
        intergenicbed = tempdir + '/intergenic.bed'
        chrsizesfile = tempdir + '/chr_sizes.tab'

        print('%s:' % genefile)
        helper.get_genic_beds(genomeanno=genefile,genomechr=chrsizesfile,verbose = verbose,infmt = infmt,genicbed=genicbed,intergenicbed=intergenicbed)
        Ntot, Ngen, Nigen = helper.estimate_mapping(bamfile = bamfile,genicbed= genicbed,intergenicbed=intergenicbed,threads=threads,verbose = verbose)
        print('Total mapped reads: %s\nGenic reads: %s (%s %%)\nIntergenic reads: %s (%s %%)' % (str(Ntot),str(Ngen),str(round(Ngen/Ntot*100,2)),str(Nigen),str(round(Nigen/Ntot*100,2))))        
        
        print('%s:' % outputfile)
        helper.get_chrsizes(tempdir = tempdir, bamfile = bamfile, outfile = chrsizesfile, verbose = verbose)
        helper.get_genic_beds(genomeanno=outputfile,genomechr=chrsizesfile,verbose = verbose,infmt = infmt,genicbed=genicbed,intergenicbed=intergenicbed)
        Ntot, Ngen, Nigen = helper.estimate_mapping(bamfile = bamfile,genicbed= genicbed,intergenicbed=intergenicbed,threads=threads,verbose = verbose)
        print('Total mapped reads: %s\nGenic reads: %s (%s %%)\nIntergenic reads: %s (%s %%)' % (str(Ntot),str(Ngen),str(round(Ngen/Ntot*100,2)),str(Nigen),str(round(Nigen/Ntot*100,2))))        
        

    if do_clean:
        print('======== Cleaning temporary directory ==========')
        clean_tmp(tempdir = tempdir)
    print('======== Done ==================================')


# but coverage calculation is only possible if one has alignment data 
