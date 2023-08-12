######################################################################################################## 
import argparse
from argparse import RawTextHelpFormatter
from geneext import helper
import logging
import os
import sys


parser = argparse.ArgumentParser(description=
"""
Program: GeneExt (Extend genes in 3' direction using single-cell RNA-seq data)
Version: 1.0
""",
                           formatter_class=RawTextHelpFormatter)
parser.add_argument('-g', default= None,help = 'Genome .gtf/.gff/.bed file.' ,required = True) 
parser.add_argument('-b', default= None,help = 'Input .bam file.')
parser.add_argument('-p', default= None,help = 'Peaks .bed file. Incompatible with -b.\nIf provided, extension is performed using specified coordinates in .bed format.\n(Can be seful in cases of FLAM-seq / Nano-3P-seq data or when manual filtering of the peaks is needed.)') 
parser.add_argument('-o', default = None, help = 'Output annotation file.\n\n\n================ Additional arguments ================\n',required = False)
parser.add_argument('-m', default = None, help = 'Maximal distance for gene extension.\nIf not set, a median length of gene (genomic span!) is used.')
parser.add_argument('-inf', default = None, help = 'Input genes file format, if None, will be guessed from a file extension.')
parser.add_argument('-ouf', default = None, help = 'Output file format, if not given, will be guessed from a file extension.')
parser.add_argument('-t', default = str('tmp'), help = 'Temporary directory. [tmp]')
parser.add_argument('-tag', default = str('GeneExt'), help = 'Tag to be added to the fake gene source and IDs so these can be easily identified downstream. [GE]')
parser.add_argument('-v', default = int(1), help = 'Verbosity level. 0,[1],2,3')
parser.add_argument('-j', default = '1', help = 'Number of cores for parallelization. [1]')
parser.add_argument('-e', default = 'new_transcript', help = 'How to extend the gene (only for .gff/.gtf files) [new_mrna]\n\t* new_transcript - creates a new transcript feature with the last exon extended\n\t* new exon - creates an extended last exon')
parser.add_argument('--clip_mode',default = 'sense',help = 'How to treat gene extension overlaps.\nsense - default,restrict overlaps into downstream genes on the same strand\nboth - restrict overlaps regardless of the strand.')
parser.add_argument('--clip5', action='store_true', help = "Use this to clip 5' overlaps between genes. The downstream gene will be clipped.\nCAVE: use carefully if the genome contains many overlapping genes (e.g. mitochondrial genomes).\n\n\n================ Orphan peaks ================\n")
parser.add_argument('--orphan',action='store_true', help = 'Whether to add orphan peaks')
parser.add_argument('--orphan_maxdist', default = int(10000), help = 'Orphan peak merging: Maximum distance between orphan peaks to merge. [100000]')
parser.add_argument('--orphan_maxsize', default = None, help = 'Orphan peak merging: Maximum size of an orphan peak cluster. Defalt: 2 x [median gene length, bp]')
parser.add_argument('--mean_coverage', action='store_true', help = 'Whether to use mean coverage for peak filtering.\nMean coverage = [ # mapping reads]/[peak width].')
parser.add_argument('--peakp',default = 25, help = 'Coverage threshold (percentile of macs2 genic peaks coverage). [1-99, 25 by default].\nAll peaks called with macs2 are required to have a coverage AT LEAST as N-th percentile of the peaks falling within genic regions.\nThis parameter allows to filter out the peaks based on the coverage BEFORE gene extension.\n\n\n================ Miscellaneous ================\n')
parser.add_argument('--subsamplebam',default = None, help = 'If set, will subsample bam to N reads before the peak calling. Useful for large datasets. Bam file should be indexed.\nDefault: None')
parser.add_argument('--report', action='store_true', help = 'Use this option to generate a PDF report.')
parser.add_argument('--keep', action='store_true', help = 'Use this to keep .bam and other temporary files in the a temporary directory. Useful for troubleshooting.')
parser.add_argument('--estimate', action='store_true', help = 'Use this to just estimate intergenic read proportion.\nUseful for quick checking intergenic mapping rate.')
parser.add_argument('--nomerge', action='store_true', help = 'Do not merge orphan peaks based on distance.')
parser.add_argument('--onlyfix', action='store_true', help = 'If set, GeneExt will only try to fix the annotation, no extension is performed')



########### Parse Arguments ################
############################################

########### Pipeline Functions ################
###############################################

# set a logger for logging:
class Logger(object):
    def __init__(self,logfilename):
        self.terminal = sys.stdout
        self.log = open(logfilename, "a")
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass    

def pipeline_error_print(x=None):
    print("\n"+x + '\n')
    os.system('cat %s/geneext/err.txt' % scriptloc)
    quit()

#######################################################################################

def parse_input_format():
    if args.inf is not None:
        infmt = args.inf
        if verbose > 0:
            print('Input: %s; Specified format: %s'  % (genefile,infmt))
    else:
        infmt = helper.guess_format(genefile)
        if verbose > 0:
            print('Input: %s, guessed format: %s' % (genefile,infmt))
    if not infmt in ['bed','gff','gtf']:
        pipeline_error_print('Unknown input format "%s"!\nPlease either make sure the format is written properly or specify it yourself using --ouf option.' % infmt)
    return(infmt)
def parse_output_format():
    if args.ouf is not None:
        outfmt = args.ouf
        if verbose > 0:
            print('Output: %s; Specified format: %s' % (outputfile,outfmt))
    else:
        # guess an output format from the name extension
        outfmt = helper.get_extension(outputfile)
        if verbose > 0:
            print('Output: %s, guessed format: %s' % (outputfile,outfmt))
    if not outfmt in ['bed','gff','gtf']:
        pipeline_error_print('Unknown output format "%s"!\nPlease either make sure the format is written properly or specify it yourself using --ouf option.' % outfmt)
    return(outfmt)      

def compare_infmt_outfmt(infmt,outfmt):
    if infmt != outfmt:
        pipeline_error_print('Please, make sure the input and the output have the same format!')

##########################################################################

def run_peakcalling():
    helper.split_strands(bamfile,tempdir,verbose = verbose,threads = threads)
    helper.run_macs2(tempdir+'/' + 'plus.bam','plus',tempdir,verbose = verbose)
    helper.run_macs2(tempdir+'/' + 'minus.bam','minus',tempdir,verbose = verbose)
    helper.collect_macs_beds(outdir = tempdir,outfile = rawpeaks,verbose = verbose)


############# Orphan peaks functions ##############
def get_orphan(genefile = None,genefile_ext_bed = None,peaks_bed = None,orphan_bed = None,infmt = None, outfmt = None,verbose = False,merge = False):
        """Store peaks not falling within new genic regions as orphan peaks """
        if infmt != 'bed':
            # remove orphan peaks overlapping extended regions:  
            helper.gxf2bed(infile = genefile,outfile = genefile_ext_bed,featuretype = 'gene')
            helper.outersect(inputbed_a = peaks_bed,inputbed_b=genefile_ext_bed,outputbed = orphan_bed,by_strand = False,verbose = verbose,f = 0.0001)
        else:
            print("Don't know how to add orphan peaks!")


# Pipeline: orphan peaks 
def run_orphan():
    genefile_ext_bed = tempdir + '/' + 'genes_ext.bed'
    orphan_bed = tempdir + '/' + 'orphan.bed'
    # if no bam file provided - no coverage filtering - all peaks not overlapping the genes. Else - coverage-filtered peaks. 
    get_orphan(genefile = outputfile, genefile_ext_bed= genefile_ext_bed,peaks_bed = peaksfilt,orphan_bed = orphan_bed,infmt = infmt, outfmt = outfmt, verbose = verbose)
    print('Orphan peaks: orphan peaks generated.')
    if do_orphan_merge:
        print('Orphan peaks: merging by distance.')
        orphan_merged_bed  = tempdir + '/' + 'orphan_merged.bed'
        helper.merge_orphan_distance(orphan_bed = orphan_bed,orphan_merged_bed = orphan_merged_bed,genic_bed = genefile_ext_bed,tempdir = tempdir,maxdist = orphan_maximum_distance,maxsize = orphan_maximum_size, verbose = verbose)
        print('Orphan peaks: merged peaks - %s' % orphan_merged_bed)
        helper.add_orphan_peaks(infile = outputfile,peaksbed=orphan_merged_bed,fmt = outfmt,verbose=verbose,tag = tag) 
    else:
        print('Orphan peaks: no merging -> directly adding orphan peaks.')
        helper.add_orphan_peaks(infile = outputfile,peaksbed=orphan_bed,fmt = outfmt,verbose=verbose,tag = tag)  


# Reporting functions 
def generate_report():
    """This function will take as an input gene extensions and plot the distributions"""
    with open(tempdir + '/_callcmd','w') as outfile:
        outfile.write(callcmd + '\n')
    cmd = 'Rscript %s/geneext/report.r %s %s %s %s %s %s %s %s' % (scriptloc,maxdist,str(int(coverage_percentile)/100),tempdir + '/_genes_peaks_closest',covfile,peaksfilt,tempdir+'/extensions.tsv',str(verbose),outputfile)
    if verbose > 1:
        print('Running:\n%s' % cmd)
    os.system(cmd)


def clean_tmp(tempdir = None):
    # clean temporary directory of big files 
    toremove = [tempdir+'/'+x for x in os.listdir(tempdir) if '.bam' in x or  x[0] == '_']
    for file in toremove:
        if verbose > 0:
            print("Removing %s" % file)
        os.remove(file)

def report_stats(filename=None,Ntot=None,Nmap=None,Ngen=None,Nigen=None,Norph=None):
    """Report mapping statistics in a way similar to cellranger"""
    o = '%s:\nTotal reads: %s\nMapped reads: %s (total: %s %%)\nGenic reads: %s (total: %s %%; mapped: %s %%)\nOrphan peaks: %s (total: %s %%; mapped: %s %%)\nIntergenic reads: %s (total: %s %%; mapped: %s %%)' % (filename,str(Ntot),str(Nmap),str(round(Nmap/Ntot*100,2)),str(Ngen),str(round(Ngen/Ntot*100,2)),str(round(Ngen/Nmap*100,2)),str(Norph),str(round(Norph/Ntot*100,2)),str(round(Norph/Nmap*100,2)),str(Nigen),str(round(Nigen/Ntot*100,2)),str(round(Nigen/Nmap*100,2)))
    return(o)


def estimate_mapping(tempdir = None,bamfile = None,infmt = None,threads = 1, verbose = False):
        # TODO: add orphan peak estimate in case there are orphan peaks? 
        # if orphan peaks present, estimate maping in the orphan peaks as well. 
        genicbed = tempdir + '/genic.bed'
        intergenicbed = tempdir + '/intergenic.bed'
        chrsizesfile = tempdir + '/chr_sizes.tab'

        helper.get_chrsizes(tempdir = tempdir, bamfile = bamfile, outfile = chrsizesfile, verbose = verbose)
        helper.get_genic_beds(genomeanno=genefile,genomechr=chrsizesfile,verbose = verbose,infmt = infmt,genicbed=genicbed,intergenicbed=intergenicbed)
        Ntot, Nmap, Ngen, Nigen = helper.estimate_mapping(bamfile = bamfile,genicbed= genicbed,intergenicbed=intergenicbed,threads=threads,verbose = verbose)
        old_rep = report_stats(genefile,Ntot,Nmap,Ngen,Nigen)
        helper.get_genic_beds(genomeanno=outputfile,genomechr=chrsizesfile,verbose = verbose,infmt = infmt,genicbed=genicbed,intergenicbed=intergenicbed)
        Ntot, Nmap, Ngen, Nigen = helper.estimate_mapping(bamfile = bamfile,genicbed= genicbed,intergenicbed=intergenicbed,threads=threads,verbose = verbose)
        new_rep = report_stats(outputfile,Ntot,Nmap,Ngen,Nigen)
        print(old_rep)
        print(new_rep)

        with open(tempdir + '/mapping_stats.txt','w') as outfile:
            outfile.write(old_rep + '\n')
            outfile.write(new_rep + '\n')
        outfile.close()

def estimate_mapping(tempdir = None,bamfile = None,genefile = None,infmt = None,threads = 1, verbose = False,orphanbed = None):
        # TODO: add orphan peak estimate in case there are orphan peaks? 
        # if orphan peaks present, estimate maping in the orphan peaks as well. 
        genicbed = tempdir + '/genic.bed'
        intergenicbed = tempdir + '/intergenic.bed'
        chrsizesfile = tempdir + '/chr_sizes.tab'

        # prepare the files
        helper.get_chrsizes(tempdir = tempdir, bamfile = bamfile, outfile = chrsizesfile, verbose = verbose)
        helper.get_genic_beds(genomeanno=genefile,genomechr=chrsizesfile,verbose = verbose,infmt = infmt,genicbed=genicbed,intergenicbed=intergenicbed)

        Ntot = helper.count_reads(bamfile=bamfile,bed = None,flags = '',threads=threads,verbose = verbose)
        Nmap = helper.count_reads(bamfile=bamfile,bed = None,flags = '-F 4',threads=threads,verbose = verbose)
        Ngen = helper.count_reads(bamfile=bamfile,bed = genicbed,flags = '',threads=threads,verbose = verbose)
        Nigen = Nmap - Ngen

        if orphanbed:
            # compute for genic regions
            genicnoorphanbed = tempdir + '/genic_no_orphan.bed'
            helper.outersect(inputbed_a=genicbed,inputbed_b=orphanbed,outputbed = genicnoorphanbed,verbose = verbose)
            Ngen = helper.count_reads(bamfile=bamfile,bed = genicnoorphanbed,flags = '',threads=1,verbose = verbose)
            Norph = helper.count_reads(bamfile=bamfile,bed = orphanbed,flags = '',threads=1,verbose = verbose)
        else:
            Norph = 0
        return(Ntot,Nmap,Ngen,Nigen,Norph)

def run_estimate(tempdir = None,bamfile = None,genefile = None,outputfile = None,infmt = None,threads = 1, verbose=False,orphanbed = None,onlyestimate = True):
    statsfile = tempdir + "/mapping_stats.txt"
    if onlyestimate:
        Ntot,Nmap,Ngen,Nigen,Norph = estimate_mapping(tempdir = tempdir,bamfile = bamfile,genefile = genefile,infmt = infmt,threads = threads, verbose =verbose,orphanbed = None)
        old_rep = report_stats(genefile,Ntot,Nmap,Ngen,Nigen,Norph)
        print(old_rep)
    else:
        Ntot,Nmap,Ngen,Nigen,Norph = estimate_mapping(tempdir = tempdir,bamfile = bamfile,genefile = genefile,infmt = infmt,threads = threads, verbose = verbose,orphanbed = None)
        old_rep = report_stats(genefile,Ntot,Nmap,Ngen,Nigen,Norph)
        Ntot,Nmap,Ngen,Nigen,Norph = estimate_mapping(tempdir = tempdir,bamfile = bamfile,genefile = outputfile,infmt = infmt,threads = threads, verbose = verbose,orphanbed = orphanbed)
        new_rep = report_stats(outputfile,Ntot,Nmap,Ngen,Nigen,Norph)
        print(old_rep)
        print(new_rep)
    with open(statsfile,'w') as ofile:
        ofile.write(old_rep)
        if not onlyestimate:
            ofile.write(new_rep)
    ofile.close()

    # depends on whether it was called only multple files or ont 




########### Main #####################
######################################

if __name__ == "__main__":
    
    args = parser.parse_args()
    bamfile = args.b
    tempdir = args.t
    verbose = int(args.v)
    peaksfile = args.p
    genefile = args.g 

    outputfile = args.o
    
    maxdist = args.m
    extension_mode = args.e
    threads = int(args.j)
    tag = args.tag
    clip_mode = args.clip_mode

    # peak coverage percentile:
    mean_coverage = args.mean_coverage
    coverage_percentile = args.peakp

    # pipeline execution:
    do_mapping = False
    do_macs2 = args.b is not None  
    do_subsample = args.subsamplebam is not None
    do_estimate = args.estimate
    do_clean = not args.keep
    do_report = args.report and do_macs2
    do_fix_only = args.onlyfix

    do_orphan = args.orphan
    do_orphan_merge =  do_orphan and not args.nomerge

    # Orphan merging defaults:
    orphan_maximum_distance =int(args.orphan_maxdist)
    orphan_maximum_size = int(args.orphan_maxsize) if args.orphan_maxsize else None

    # Clipping 5'-overlaps 
    do_5clip = args.clip5
    #do_5clip = True

    # Logging 
    if outputfile:
        logfilename = outputfile + ".GeneExt.log"
    else:
        logfilename = "GeneExt.log"
    sys.stdout = Logger(logfilename)

    scriptloc = os.path.dirname(os.path.realpath(__file__))
    callcmd = 'python ' + os.path.basename(__file__) + ' '+ " ".join(["-"+str(k)+' '+str(v) for k,v in zip([arg for arg in vars(args)],[getattr(args,arg) for arg in vars(args)]) if v ])
    print(callcmd)
    # print the values of the pipeine settings


    # If fix_only set - report only doing genome annotation fixes:
    if do_fix_only:
        print('================================================')
        print("CAVE: --onlyfix is set. Only fixing the genome annotation file %s, no extension will be performed!" % (args.g))
    print('======== Preflight checks ======================')

    if not do_fix_only:
        # parse input and output formats - run the pipeline accordingly
        if peaksfile is None and bamfile is None:
            pipeline_error_print("Please, specify either alignment [-b] or peaks file [-p]!")

        if bamfile is not None:
            if os.path.isfile(bamfile):
                print('Alighment file ... OK')
            else:
                pipeline_error_print('Specified alignment file does not exist!')

    if not do_estimate:
        # check output file:
        if outputfile is None and not do_estimate:
            pipeline_error_print('Please, specify the output file [-o]!')

        elif peaksfile is not None:
            if os.path.isfile(peaksfile):
                print('Found a peaks file, skipping peak calling ...')
            else:
                pipeline_error_print('Specified peaks file does not exist!\nPlease, specify either a valid .bam file for macs2 or a peaks file.')

        if peaksfile is not None and bamfile is not None:
            pipeline_error_print('Please, specify either a .bam file with reads [-b] or a peaks file [-p] but not both at the same time!')

        if genefile is None:
            pipeline_error_print('Missing genome annotation file [-g]!')

        elif not os.path.isfile(genefile):
            pipeline_error_print('Genome annotation file .... DOES NOT EXIST!')
        else:
            print('Genome annotation file .... OK')


    else:
        # just check the bam file:
        if bamfile is None:
            pipeline_error_print('Can not estimate mapping without an alignment file!')
    # set temporary directory:
    if verbose > 0:
        print("Temporary directory: %s" % tempdir)
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)
        if verbose > 0:
                print('Directory created: %s' % tempdir)
    else:
        print('Temporary directory exists. Overwriting!')

    ####################################################
    # Parse input / output formats
    # The output format should also be parsed separately 
    infmt = parse_input_format()
    outfmt = parse_output_format()
    if(outfmt == 'bed'):
        pipeline_error_print("I can't output .bed yet!")
    compare_infmt_outfmt(infmt,outfmt)
    ####################################################
    # Check and fix the input file
    if infmt in ['gff','gtf']:
        features = helper.get_featuretypes(genefile)
        if not 'transcript' in features:
            print('Genome annotation warning: Could not find "gene" features in %s! Trying to fix ...' % genefile)
            if 'mRNA' in features:
                fpref = '.'.join(genefile.split('/')[-1].split('.')[:-1])
                fext = genefile.split('/')[-1].split('.')[-1]
                genefilewmrna = tempdir + '/' + fpref + '_mRNA2transcript.' + fext
                #genefilewmrna = tempdir + '/' + genefile.split('/')[-1].replace('.' + infmt,'_mRNA2transcript.' + infmt)
                print('Found "mRNA" features - renaming as transcripts ...')
                helper.mRNA2transcript(infile = genefile,outfile = genefilewmrna, verbose = verbose)
                genefile = genefilewmrna
            else:
                fpref = '.'.join(genefile.split('/')[-1].split('.')[:-1])
                fext = genefile.split('/')[-1].split('.')[-1]
                genefilewmrna = tempdir + '/' + fpref + '_addtranscripts.' + fext
                #genefilewmrna = tempdir + '/' + genefile.split('/')[-1].replace('.' + infmt,'_addtranscripts.' + infmt)
                helper.add_transcript_features(infile = genefile, outfile = genefilewmrna,verbose = verbose)
        if not 'gene' in features:
            print('Could not find "gene" features in %s! Trying to fix ...' % genefile)
            fpref = '.'.join(genefile.split('/')[-1].split('.')[:-1])
            fext = genefile.split('/')[-1].split('.')[-1]
            genefilewgenes = tempdir + '/' + fpref + '_addgenes.' + fext
            #genefilewgenes = tempdir + '/' + genefile.split('/')[-1].replace('.' + infmt,'_addgenes.' + infmt)
            helper.add_gene_features(infile = genefile,outfile = genefilewgenes,infmt = infmt,verbose = verbose)
            print('Fix done, annotation with gene features: %s' % genefilewgenes )
            genefile = genefilewgenes
        # Fix 5'overlaps 
        if do_5clip:
            print("Clipping 5' overlaps in genes using %s cores..." % str(threads))
            # rename the file properl y
            fpref = '.'.join(genefile.split('/')[-1].split('.')[:-1])
            fext = genefile.split('/')[-1].split('.')[-1]
            genefile5clip = tempdir + '/' + fpref + '_5clip.' + fext
            helper.clip_5_overlaps(infile = genefile,outfile = genefile5clip,threads = threads,verbose = verbose)
            print("Fixed 5' overlaps in genes: %s -> %s" % (genefile,genefile5clip))
            genefile = genefile5clip
    print('Checks done.')

    ##################################################
    if not do_fix_only:
        # parse input file format: 
        if not do_estimate:    
            # if -m is not set, get a median gene size:
            if not maxdist:
                maxdist = helper.get_median_gene_length(inputfile = genefile,fmt = infmt)
                if verbose:
                    print('Maximum allowed extension length is not set, getting median size of the gene - %s bp.' % str(maxdist))
            # if maximum size for orphan peak is not set, set it to the median gene size:
            if do_orphan_merge:
                if not orphan_maximum_size:
                    orphan_maximum_size = 2 * helper.get_median_gene_length(inputfile=genefile,fmt = 'gff')

        # 0. MAPPING - not implemented     
            if do_mapping:
                print('======== Running mapping =======================')
                raise(NotImplementedError())
        # 0.1 BAM SUBSAMPLING  
            if do_subsample:
                print('======== Running subsampling ===================')
                subsampled_bam = tempdir + '/subsampled.bam' 
                nsubs = int(args.subsamplebam)
                if bamfile:
                    if not os.path.isfile(bamfile + '.bai'):
                        if verbose > 0:
                            print('Indexing %s' % bamfile)
                        helper.index_bam(bamfile,verbose = verbose,threads=threads)
                # check here if it's an integer
                helper.subsample_bam(inputbam = bamfile,outputbam = subsampled_bam,nreads = nsubs,verbose = verbose,threads=threads)
                # now, replace for downstream:
                if verbose > 0:
                    print('Indexing %s' % bamfile)
                helper.index_bam(subsampled_bam,verbose = verbose,threads=threads)
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
                if not os.path.isfile(bamfile + '.bai') and 'subsample' in bamfile:
                    if verbose > 0:
                        print('Indexing %s' % bamfile)
                    helper.index_bam(bamfile,verbose = verbose,threads=threads)
                if verbose > 0:
                    print('Computing coverage ...')
                helper.get_coverage(inputbed_a=peaksfile,input_bam = bamfile,outputfile = covfile,verbose = verbose,mean = mean_coverage,threads = threads)
                # get the peaks overlapping genes:
                if verbose > 0:
                    print('Getting genic peaks ...')
                genicpeaksfile = tempdir + '/genic_peaks.bed'
                helper.intersect(inputbed_a = covfile, inputbed_b = genefile,outputbed=genicpeaksfile,by_strand = True,verbose = verbose)
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
                # f=1 will filter out only the peaks fully contained within genes!!!
                helper.outersect(inputbed_a = covfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose,f = 1)
                helper.filter_by_coverage(inputfile = peaksfilt,outputfile = peaksfiltcov,threshold = count_threshold,verbose = True)
                peaksfilt = peaksfiltcov
            else:
                # Simply remove peaks overlapping genes:
                peaksfilt = tempdir + '/' + 'allpeaks_noov.bed'
                if verbose > 0:
                    print('Removing peaks overlapping genes ... => %s' % peaksfilt)
                #outersect_peaks(genefile = genefile, peaksfile = peaksfile, outputbed = peaksfilt, verbose = verbose)
                helper.outersect(inputbed_a = peaksfile,inputbed_b = genefile,outputbed=peaksfilt,by_strand = True, verbose = verbose)

        # 3. Extend genes 
            print('======== Extending genes =======================')
            helper.extend_genes(genefile = genefile,peaksfile = peaksfilt,outfile = outputfile,maxdist = int(maxdist),temp_dir = tempdir,verbose = verbose,extension_type = extension_mode,infmt = infmt,outfmt = outfmt,tag = tag,clip_mode = clip_mode)
            
        # 4. Add orphan peaks
            if do_orphan:
                print('======== Adding orphan peaks ===================')
                #run_orphan(infmt = infmt,outfmt = outfmt,verbose = verbose,merge = do_orphan_merge)
                run_orphan()
            if do_report:
                print('======== Creating PDF report =======================')
                generate_report()

        # 5. Estimate intergenic mapping
            print('======== Estimating intergenic mapping =========')
            if bamfile:
                if do_orphan:
                    run_estimate(tempdir = tempdir,bamfile = bamfile,genefile = genefile,outputfile = outputfile,infmt = infmt,threads = threads, verbose=verbose,orphanbed = tempdir + '/orphan_merged.bed',onlyestimate = False)
                else:
                    run_estimate(tempdir = tempdir,bamfile = bamfile,genefile = genefile,outputfile = outputfile,infmt = infmt,threads = threads, verbose=verbose,orphanbed = None,onlyestimate = False)
            else:
                print("No bamfile specified - omitting mapping estimation")
            if do_clean:
                print('======== Cleaning temporary directory ==========')
                clean_tmp(tempdir = tempdir)
        elif bamfile:
            print('--estimate is set. Skipping extension, estimating mapping rates for %s with %s' % (genefile,bamfile))
            run_estimate(tempdir = tempdir,bamfile = bamfile,genefile = genefile,outputfile = outputfile,infmt = infmt,threads = threads, verbose=verbose,orphanbed = None,onlyestimate = True)
        print('======== Done ==================================')

