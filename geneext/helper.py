import subprocess
import os
import re
import numpy as np
import gffutils
import pysam
import multiprocessing
from functools import partial
import time

# macs2_helper

def split_strands(bamfile,outdir,verbose = False,threads = 1):
    """Using samtools, separate input bam file by strands"""
    cmd='samtools view  -@ %s -F 16 %s -b > %s' % (threads,bamfile,outdir + '/plus.bam')
    if verbose > 1:
        print('Running:\n\t%s' % cmd) 
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if verbose > 1:
        print(ps)
    cmd='samtools view  -@ %s -f 16 %s -b > %s' % (threads,bamfile,outdir + '/minus.bam')
    if verbose > 1:
        print('Running:\n\t%s' % cmd) 
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if verbose > 1:
        print(ps)

def run_macs2(bamfile,peaks_prefix,outdir,verbose = False):
    """This function launches MACS2 to call peaks from a .bam file"""
    cmd = ("macs2","callpeak","-t", bamfile ,"-f", "BAM", "--keep-dup", "20","-q", "0.01" , "--shift", "1" ,"--extsize", "20", "--broad", "--nomodel", "--min-length", "30", "-n",peaks_prefix,"--outdir", outdir)
    try:
        if verbose > 1:
            print('Running:\n\t%s' % ' '.join(cmd))
        ps = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        ps.check_returncode()
    except subprocess.CalledProcessError as e:
        print ( "\n########## macs2 FAILED ##############\nreturn code: ", e.returncode, "\nOutput: ", e.stderr.decode("utf-8") )
        raise

def collect_macs_beds(outdir,outfile,verbose = False):
    """This script collects macs2 output into a single bed file of peaks"""
    cmd = "cat %s/plus_peaks.broadPeak %s/minus_peaks.broadPeak | cut -f 1-6 | awk 'BEGIN{OFS=__\\t__}{if($4~/plus/){$6=__+__}else{$6=__-__};print $0}' | bedtools sort -i - > %s" % (outdir,outdir,outfile)
    cmd = cmd.replace('__','"')
    if verbose > 1 :
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if verbose > 1:
        print('Done collecting beds: %s' % (outfile))


def index_bam(infile,verbose = False,threads = 1):

    cmd = 'samtools index -@ %s %s' % (threads,infile)
    if verbose>1:
        print('Running:\n\t%s' % cmd)
    os.system(cmd)

# Coverage functions   
def count_reads_in_region(region, aln):
    read_count = aln.count(contig=region.chrom, start=region.start, stop=region.end, region=None, until_eof=False, read_callback='nofilter', reference=None, end=None)
    return read_count

def compute_mean_coverage(region, aln):
    read_count = aln.count(contig=region.chrom, start=region.start, stop=region.end, region=None, until_eof=False, read_callback='nofilter', reference=None, end=None)
    mean_coverage = read_count / (region.end - region.start)
    return mean_coverage


################################################################
from typing import List
from multiprocessing import Pool
import time
from functools import partial

def process_region(region, aln, mean):
    if mean:
        read_count = compute_mean_coverage(region, aln)
    else:
        read_count = count_reads_in_region(region, aln)
    return (region, read_count)

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def func(regions, input_bam, mean):
    aln = pysam.AlignmentFile(input_bam, "rb")
    return [(region, process_region(region=region, aln=aln, mean=mean)) for region in regions]


def get_coverage(inputbed_a, input_bam, outputfile, verbose=0, mean=False, threads=1):
    start = time.time()
    bed = parse_bed(inputbed_a)
    # split into chunks for threads
    chunks = [x for x in split(bed, threads)]
    print("mean chunk size: %s" % np.mean(np.array([len(x) for x in chunks])))
    with open(outputfile, "w") as fout:
        if verbose > 1:
            print(f"Computing coverage for {str(len(bed))} peaks with {threads} threads...")
        with multiprocessing.Pool(threads) as pool:
            for result in pool.imap_unordered(partial(func, input_bam=input_bam, mean=mean), chunks):
                for region, (region_obj, read_count) in result:
                    fout.write(
                        "\t".join(
                            [
                                region_obj.chrom,
                                str(region_obj.start),
                                str(region_obj.end),
                                region_obj.id,
                                "0",
                                region_obj.strand,
                                str(read_count),
                            ]
                        )
                        + "\n"
                    )
    end = time.time()
    e = end - start
    print('Done computing coverage: %s' % (outputfile))
    print(f"{round(e/len(bed)*1000,2)} ms per region")
###############################################################################

def get_coverage_percentile(inputfile = None,percentile = None, verbose = False):
    """Given an input bed file with a coverage, get a coverage percentile"""
    percentile = int(percentile)
    if percentile > 0:
        if verbose > 0:
                print('Getting a %s-th percentile ...' % percentile)
        cmd = "cut -f 7 %s  | sort -n | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p=(c/100*%s); p=p%%1?int(p)+1:p; print a[c-p-1]}'" % (inputfile,str(100-percentile))
        if verbose > 1:
            print('Running:\n\t%s' % cmd)
        ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        out = ps.communicate()[0].decode("utf-8").rstrip()
        return(out)
    else:
        return('0')


def filter_by_coverage(inputfile = None,outputfile = None,threshold = None,verbose = False):
    """
    Filter bed file by the last column
    """
    cmd = "awk '$NF>=%s' %s | cut -f 1-7 > %s" % (str(threshold),inputfile,outputfile)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    # Report how many peaks have been retained:
    if verbose > 1:
        ps = subprocess.Popen('wc -l %s' % outputfile,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        n = ps.communicate()[0].decode("utf-8").rstrip().split(' ')[0]
        ps = subprocess.Popen('wc -l %s' % inputfile,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        N = ps.communicate()[0].decode("utf-8").rstrip().split(' ')[0]
        print('Retained %s/%s (%s %%)' % (str(n),str(N),str(round(int(n)/int(N)*100,2))))

def outersect(inputbed_a,inputbed_b,outputbed,by_strand = True,verbose = False,f = 0.0001):
    """This function returns non-overlapping peaks"""
    if by_strand:
        strand = '-s'
    else:
        strand = ''
    cmd = "bedtools intersect -f %s -a %s -b %s %s -v > %s" % (str(f),inputbed_a,inputbed_b,strand,outputbed)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        

def intersect(inputbed_a,inputbed_b,outputbed,by_strand = True,verbose = False):
    """This function returns non-overlapping peaks"""
    if by_strand:
        strand = '-s'
    else:
        strand = ''
    cmd = "bedtools intersect -wa -a %s -b %s %s > %s" % (inputbed_a,inputbed_b,strand,outputbed)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

# Parse helper
def get_extension(filepath):
    fmt = os.path.splitext(filepath)[-1][1:]
    return(fmt)

def guess_format(filepath,verbose = False):
    """Guess file exension from a header"""
    #with open(filepath) as file:
    #    head = [next(file).rstrip() for x in range(5)]
    with open(filepath) as file:
        head = next(file).rstrip()
        firstchar = head[0]
        while firstchar == '#':
            head = next(file).rstrip()
            firstchar = head[0]
    file.close()
    
    
    nfields = len(head.split('\t'))
    if verbose > 1:
        print('%s fields in file' % nfields)

    if nfields == 9:
        if verbose > 1:
            print('File has 9 fields - .gff/.gtf?')
        if 'ID=' in head.split('\t')[8]:
            return('gff')
        elif '_id "' in head.split('\t')[8]:
            return('gtf')
        else:
            raise(ValueError('Can not determine the format of the file (counted %s fields) - is it gtf/gff?\nUse -inf/-ouf to specify the format!\n\n%s\n\n' % (nfields,head)))
    elif nfields in [6,7]:
        return('bed')
    else:
        raise(ValueError('Can not determine the format of the file (counted %s fields) - is it gtf/gff or bed?\nUse -inf/-ouf to specify the format!\n\n%s\n\n' % (nfields,head)))

def parse_bed(infile):
    with open(infile) as file:
        lines = [line.rstrip().split('\t') for line in file if not '#' in line]
        regs = [Region(chrom = x[0],start = int(x[1]),end = int(x[2]),id = str(x[3]),strand = str(x[5]),score = 0) for x in lines]
    file.close()
    return(regs)
    

def parse_gff(infile,featuretype = None):
    def gff_get_ID(x):
    # get ID from .gff 9th field
        m = re.search('ID=(.+?);', x)
        # Aha! Newly
        if m:
            return(m.group(1))
        else:
            return(None)

    with open(infile) as file:
            lines = [line.rstrip().split('\t') for line in file if not '#' in line]
            if not featuretype:
                regs = [Region(chrom = x[0],start = int(x[3]),end = int(x[4]),id = gff_get_ID(str(x[8])),strand = str(x[6])) for x in lines]
            else:
                regs = [Region(chrom = x[0],start = int(x[3]),end = int(x[4]),id = gff_get_ID(str(x[8])),strand = str(x[6])) for x in lines if x[2] == featuretype]
    if len(regs) == 0:
        raise(ValueError('No %s found in file %s!' % (featuretype if featuretype else 'lines',infile)))
    else:
        return(regs)


def parse_gtf(infile,featuretype = None):
    def gtf_get_ID(x):
        """x should be a string from the 9th field of a .gtf file. """
        o = [y.strip() for y in x.split(';') if len(y.strip())]
        d = {k.split(' ')[0]:k.split(' ')[1].replace('"','') for k in o}
        if 'ID' in d.keys():
            return(d['ID'])
        elif 'gene_id' in d.keys():
            return(d['gene_id'])
        else:
            raise(ValueError('No gene ID found in file %s:\n%s' % (infile,''.join(x))))

    with open(infile) as file:
            lines = [line.rstrip().split('\t') for line in file if not '#' in line]
            if not featuretype:
                regs = [Region(chrom = x[0],start = int(x[3]),end = int(x[4]),id = gtf_get_ID(str(x[8])),strand = str(x[6])) for x in lines]
            else:
                regs = [Region(chrom = x[0],start = int(x[3]),end = int(x[4]),id = gtf_get_ID(str(x[8])),strand = str(x[6])) for x in lines if x[2]==featuretype]
    if len(regs) == 0:
        raise(ValueError('No %s found in file %s!' % (featuretype if featuretype else 'lines',infile)))
    else:
        return(regs)


def _guess_format(filepath,fmt = None,featuretype = None):
    """Guess file extension"""
    if not fmt:
        # get extension from file name:
        extension = os.path.splitext(filepath)[-1][1:]
    else:
        extension = fmt
    # check if bed:
    with open(filepath) as myfile:
        # skip lines with comments;
        head = next(myfile)
        while '#' in head:
            head = next(myfile)
    n_fields = len(head.split('\t'))
    if extension == 'bed':
        if n_fields in [6,7]:
            return('bed')
        else:
            raise ValueError('Guessed %s extension from %s. Number of fields is wrong (%s).' % (extension,filepath,n_fields))
    if extension == 'gtf' or extension == 'gff':
        if n_fields == 9:
            if extension == 'gff':
                return('gff')
            else:
                return('gtf')
        else:
            raise ValueError('Guessed %s extension from %s. Number of fields is wrong (%s).' % (extension,filepath,n_fields))
    else:
        raise ValueError('Unknown file extension in file %s. Please, rename your file: <filename>[.bed/.gff/.gtf]' % (filepath))

def check_ext_read_file(filepath,featuretype = None):
    """Get file extension and load a file as a list of Region objects."""
    extension = guess_format(filepath)
    if not extension in ['gff','bed','gtf']:
        raise(ValueError("Please, use bed / gff / gtf file!"))
    else:
        if extension == 'bed':
            return(parse_bed(filepath))
        elif extension == 'gff':
            return(parse_gff(filepath,featuretype))
        else:
            return(parse_gtf(filepath,featuretype))

def write_bed(outfile,regs):
    """Write a set of Region objects into the bed file"""
    with open(outfile,'w') as file:
        for reg in regs:
            file.write('\t'.join([reg.chrom,str(reg.start),str(reg.end),reg.id,str(reg.score),reg.strand]) + '\n')  


def gxf2bed(infile,outfile,featuretype = None):
    """This function loads the gff/gtf file and returns a bed file"""
    print(infile)
    regs = check_ext_read_file(infile,featuretype=featuretype)
    write_bed(outfile,regs)


# Region helper

class Region:

    def __init__(self, chrom,start, end, strand, id = None,score = None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.id = id
        if score:
            self.score = score
        else:
            self.score = 0 
 
    def show(self):
        return('%s\t%s\t%s\t%s\t%s' % (self.chrom,self.start,self.end,self.id,self.strand))


    def length(self):
        return(self.end-self.start)

    @staticmethod
    # this method checks if 2 regions are overlapping 
    def is_overlapping(region_a,region_b):
        return(not((region_a.chrom == region_b.chrom and region_a.end < region_b.start) or (region_a.chrom == region_b.chrom and  region_b.end < region_a.start)))

    @staticmethod
    # CAVE: there may be a problem here if 2 regions are overlapping - distance would be defined for non-overlapping regions 
    # only!
    def get_distance(region_a,region_b):
        if not Region.is_overlapping(region_a,region_b):
            if region_a.chrom == region_b.chrom:
                if region_a.strand and region_b.strand == '+':
                # compare the ends 
                    return(abs(region_a.end - region_b.end))
                elif region_a.strand and region_b.strand == '-':
                    return(abs(region_a.start - region_b.start))
                else:
                    # Not on the same strand 
                    return(None)
            else:
                return(None)
        else:
            return(None)
    @staticmethod
    def get_closest(query,target):
        """Return the region in target list closes to the query"""
        distances = {y.id:Region.get_distance(query,y) for i in range(len(target)) for y in target if Region.get_distance(query,y)}
        min_dist = min(distances,key = distances.get)
        if min_dist:
            #return({x.id:x for x in target}[min_dist])
            return(min_dist)
    @staticmethod
    # check if the regions are downstream:
    def a_is_downstream_b(region_a,region_b):
        if not Region.is_overlapping(region_a,region_b):
            if region_a.strand == region_b.strand == '+':
                return(region_a.start > region_b.end)
            elif region_a.strand == region_b.strand == '-':
                return(region_a.end < region_b.start)
            else:
                return(None)
        else:
            return(None)
    @staticmethod
    # check if the regions are downstream:
    def a_is_upstream_b(region_a,region_b):
        if not Region.is_overlapping(region_a,region_b):
            if region_a.strand == region_b.strand == '+':
                return(region_a.end < region_b.start)
            elif region_a.strand == region_b.strand == '-':
                return(region_a.start > region_b.end)
            else:
                return(None)
        else:
            return(None)

    @staticmethod
    # get closest downstream region
    def get_closest_downstream(query,target):
        downstream = [x for x in target if Region.a_is_downstream_b(x,query)]
        if downstream:
            return(Region.get_closest(query,downstream))
        else:
            return(None)
    @staticmethod
    # get closest upstream region
    def get_closest_upstream(query,target):
        upstream = [x for x in target if Region.a_is_upstream_b(x,query)]
        if upstream:
            return(Region.get_closest(query,upstream))
        else:
            return(None)


# GXF helper 
def gffutils_import_gxf(filepath,verbose = False):
    if verbose > 0:
        print('\tgffutils: creating a database in memory (may take a while for a big .gff/.gtf)...')
    db = gffutils.create_db(filepath, ':memory:',disable_infer_genes=True,disable_infer_transcripts=True, merge_strategy = 'create_unique',transform = gffutils_transform_func,keep_order = True)
    if verbose > 1:
        print('======== Features loaded: ======================')
        for a,b in zip(db.featuretypes(),[db.count_features_of_type(x) for x in db.featuretypes()]):
            print(a + ':' + str(b))
    return(db)

def replace_gff_gtf(x):
    return(x.replace(';;',';').replace('=',' "').replace(';','"; ')+'"')

def str_gtf(f,attributes = ['gene_id','transcript_id','three_prime_ext']):
    # creates a feature line from gffutils feature
    # clean attributes - retain only gene_id, transcript_id
    ats = [x for x in f.attributes if x in attributes]
    at = '; '.join(['%s "%s"' % (k,v)  for k,v in zip([x for x in ats],[f[x][0] for x in ats])])
    o = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.end),'.',f.strand,'.',at])

    return(o)

def str_gff(f,attributes= ['ID','Parent']):
    # creates a feature line from gffutils feature
    # clean attributes - retain only gene_id, transcript_id
    ats = [x for x in f.attributes if x in attributes]
    at = '; '.join(['%s=%s' % (k,v)  for k,v in zip([x for x in ats],[f[x][0] for x in ats])])
    o = '\t'.join([f.chrom,f.source,f.featuretype,str(f.start),str(f.end),'.',f.strand,'.',at])

    return(o)

def extend_gff(db,extend_dictionary,output_file,extension_mode,tag,verbose = False,infmt = None,outfmt = None):
    """ This function will extend .gff file by either adding new mRNA or by extending the last exon\n
    extension_type - the way in which to handle .gff/.gtf file
        - new_transcript - creates new fake mRNA with the last exon extended to match the peak 
        - new_exon - just adds an artificial last exon with an extension 
        - replace_transcript - not implemented, replaces an original transcript with the one with extended last exon 
    
    
    """


    with open(output_file, 'w') as fout:
        cnt=0
        for d in db.directives:
            # write down output 
            fout.write('##{0}\n'.format(d))
        for i,feature in enumerate(db.features_of_type("gene")):
            if infmt == 'gff':
                if not 'ID' in feature.attributes:
                    print(feature)
                    exit("no gene_id attribute found for gene %s! ( see the feature line above)" % feature.id)
            elif infmt == 'gtf':
                if not 'gene_id' in feature.attributes:
                    print(feature)
                    exit("no gene_id attribute found for gene %s! ( see the feature line above)" % feature.id)    
            n_exons = len([x for x in db.children(db[feature.id],featuretype='exon')])
            n_transcripts = len([x for x in db.children(db[feature.id],featuretype = 'transcript')])
            if n_exons and feature.id in extend_dictionary.keys(): 
                # dictoinary with written exons:
                written_exons = []
                if verbose > 2:
                    print(feature.id)
                    print("%s strand." % db[feature.id].strand)
                    print("%s transcripts found" % n_transcripts )
                    print("%s exons found" % n_exons)
                    print('Find mRNA with the most downstream exon ...')

                    if n_exons == 0 or n_transcripts == 0:
                        print('No exons or transcripts found in the file!\nFeature types found in the annotation file:\n')
                        print('\n'.join([x for x in db.featuretypes()]))
                        continue
                # Identify the most downstream exon per gene 
                if db[feature.id].strand == '+':
                    # update gene bound:
                    # gene end is probably not updated properly
                    feature.end = db[feature.id].end + extend_dictionary[feature.id]
                    max_end = max([f.end for f in db.children(db[feature.id],featuretype='exon')])
                    last_exon = [x for x in db.children(db[feature.id],featuretype='exon') if x.end == max_end ][0]
                    last_exon.end = last_exon.end + extend_dictionary[feature.id]
                    if verbose > 2:
                        print("Gene end chage: %s --> %s; %s" % (str(last_exon.end),str(last_exon.end+extend_dictionary[feature.id]),str(extend_dictionary[feature.id])))
                elif db[feature.id].strand == '-':
                    feature.start = db[feature.id].start - extend_dictionary[feature.id]
                    max_end = min([f.start for f in db.children(db[feature.id],featuretype='exon')])
                    last_exon = [x for x in db.children(db[feature.id],featuretype='exon') if x.start == max_end ][0]
                    last_exon.start = last_exon.start - extend_dictionary[feature.id]
                    if verbose > 2:
                        print("Gene end chage: %s --> %s; %s" % (str(last_exon.start),str(str(last_exon.start-extend_dictionary[feature.id])),str(extend_dictionary[feature.id])))
                
                if infmt == 'gff':
                    if not 'Parent' in last_exon.attributes:
                        raise(ImportError('Missing "Parent" feature for the exon: %s' % last_exon.id))
                    else:
                        mrna_id = last_exon['Parent'][0]
                    if verbose > 2:
                        print("mRNA with the most downstream exon: %s" % mrna_id)
                else:
                    if not 'transcript_id' in last_exon.attributes:
                        raise(ImportError('Missing "transcript_id" feature for the exon: %s' % last_exon.id))
                    else:
                        mrna_id = last_exon['transcript_id'][0]
                    if verbose > 2:
                        print("mRNA with the most downstream exon: %s" % mrna_id)    

                if extension_mode == 'new_transcript':
                    if infmt == 'gff':
                        ###### Create fictional mRNA ########
                        # change last_exon id 
                        gene_id = feature.id
                        if outfmt == 'gff':
                            last_exon['ID'] = [tag + '~lastexon_'+last_exon['ID'][0]]
                            # now, change the parent
                            last_exon['Parent'] = [tag + '~'+mrna_id]
                            # duplicate last exon, change the metadata:
                            last_exon.source = last_exon.source + '~' + tag
                            last_exon['three_prime_ext'] = str(extend_dictionary[feature.id])
                        elif outfmt == 'gtf':
                            # need to add gtf attributes: 
                            last_exon['gene_id'] = gene_id
                            last_exon['transcript_id'] = [tag + '~'+mrna_id]
                            last_exon.source = last_exon.source + '~' + tag
                            last_exon['three_prime_ext'] = str(extend_dictionary[feature.id])   
                                   
                            #raise(NotImplementedError('gff-> gtf'))
                        # CAVE: mrna ranges should be also updated
                        mrna = db[mrna_id]
                        if mrna.featuretype == 'gene' and verbose > 2:
                            print('Gene picked instead of mRNA!')
                            # the feature can not be a parent of itself.
                            mrna['Parent'] = mrna_id
                            mrna.featuretype = 'transcript' 
                        # CAVE: in case of the same ids, the gene feature is actually picked
                        
                        # Add transcript ID if missing 
                        # Write down original gene and mRNA:
                        if outfmt == 'gff':
                            fout.write(str(feature)+'\n')
                            fout.write(str(mrna)+'\n')
                            for child_exon in db.children(mrna_id,featuretype = 'exon'):
                                fout.write(str(child_exon)+'\n')

                        elif outfmt == 'gtf':
                            feature['gene_id'] = feature['ID']
                            #fout.write(replace_gff_gtf(str(feature)))
                            fout.write(str_gtf(feature)+'\n')
                            
                            # ADD GENEID TO MRNA 
                            mrna['gene_id'] = mrna['Parent']
                            mrna['transcript_id'] = mrna['ID']
                            fout.write(str_gtf(mrna)+'\n')
                            for exon in db.children(mrna_id,featuretype = ['exon']):
                                exon['gene_id'] = mrna['gene_id']
                                exon['transcript_id'] = mrna['transcript_id']
                                fout.write(str_gtf(exon)+'\n')
                     
                        else:
                            print('Dont know how to write an output')

                        # now add the fake mRNA:
                        if mrna.strand == "+":
                            mrna.end = last_exon.end
                        elif mrna.strand == "-":
                            mrna.start = last_exon.start
                        # create fictional mrnaid
                        mrna.id = tag + '~' + mrna.id
                        mrna.source = mrna.source + '~' + tag
                        mrna['ID'][0] = mrna.id
                        # add extension length 
                        mrna['three_prime_ext'] = str(extend_dictionary[feature.id])

                        ####### Write down the gene and extended mrna ######
                        if verbose > 2:
                            print('adding %s: [%s/%s]' % (mrna.id,cnt,len(extend_dictionary)))
                        if outfmt == infmt:
                            fout.write(str(mrna)+'\n')
                        elif infmt == 'gff' and outfmt == 'gtf':
                            fout.write(str_gtf(mrna)+'\n')
                            if mrna.strand == '-':
                                fout.write(str_gtf(last_exon) + '\n')
                            for exon in db.children(mrna_id,featuretype = 'exon'):
                                if not exon.id == last_exon.id:
                                    exon['gene_id'] = mrna['gene_id']
                                    exon['transcript_id'] = mrna['transcript_id']
                                    exon.source = exon.source + '~' + tag
                                    fout.write(str_gtf(exon)+'\n')
                            if mrna.strand == '+':
                                fout.write(str_gtf(last_exon) + '\n')
                        cnt += 1
                    
                    elif infmt == 'gtf':
                        ###### Create fictional mRNA ########
                        new_mrna_id = tag + '~' + mrna_id
                        gene_id = feature.id

                        # now, change the parent transcript
                        last_exon['transcript_id'] =  new_mrna_id
                        # duplicate last exon, change the metadata:
                        last_exon.source = last_exon.source + '~' + tag
                        last_exon['three_prime_ext'] = str(extend_dictionary[gene_id])
                        # CAVE: mrna ranges should be also updated

                        mrna = db[mrna_id]
                        mrna.featuretype = 'transcript'
                        # if missing transcript_id (because they are the same as genes, add manually)
                        if not 'transcript_id' in [x for x in mrna.attributes]:
                            mrna['transcript_id'] = mrna['gene_id']

                        # Write down original gene and mRNA:
                        if outfmt == 'gtf':
                            fout.write(str(feature)+'\n')
                            fout.write(str(mrna)+'\n')
                            for child_exon in db.children(mrna_id,featuretype = ['exon']):
                                if not str(child_exon) in written_exons:
                                    fout.write(str(child_exon)+'\n')
                                    written_exons = written_exons + [str(child_exon)]

                        if mrna.strand == "+":
                            mrna.end = last_exon.end
                        elif mrna.strand == "-":
                            mrna.start = last_exon.start

                        
                        # add information on gene extension 
                        mrna['three_prime_ext'] = str(extend_dictionary[feature.id])

                        # ADD FAKE ID 
                        mrna['transcript_id'] = new_mrna_id

                        if not 'gene_id' in mrna.attributes:
                            mrna['gene_id'] = mrna_id
                        mrna.source = mrna.source + "~" + tag
                        #  write down this updated mRNA:
                        fout.write(str(mrna)+'\n')

                        ####### Write down the gene and extended mrna ######
                        if verbose > 2:
                            print('adding %s: [%s/%s]' % (mrna.id,cnt,len(extend_dictionary)))
                        cnt += 1
                        for child_exon in db.children(mrna_id,featuretype = 'exon'):
                            if not child_exon.id == last_exon.id and not str(child_exon) in written_exons:
                                child_exon.source = child_exon.source + '~' + tag
                                child_exon['transcript_id'] = new_mrna_id
                                fout.write(str(child_exon)+'\n')
                                written_exons = written_exons + [str(child_exon)]
                        fout.write(str(last_exon) + '\n')
                        written_exons  = written_exons + [str(last_exon)]
                elif extension_mode == 'new_exon':
                    raise(NotImplementedError())
                    ###### Only add a fictional exon ########
                elif extension_mode == 'replace_transcript':
                    raise(NotImplementedError('Implement replacing old mRNA with a new one!'))
                else:
                    raise(ValueError('Unknown extension type!'))
                
                # Write remaining transcripts per gene:
                for transcript in db.children(feature.id,featuretype = 'transcript'): 
                    if infmt == 'gtf':
                        transcript_id = transcript['transcript_id'][0]
                    elif infmt == 'gff':
                        transcript_id = transcript.id
                    
                    if transcript_id != mrna_id:
                        if infmt == outfmt:
                            fout.write(str(transcript) + '\n')
                        elif infmt == 'gff' and outfmt == 'gtf':
                            transcript['gene_id'] = feature.id
                            transcript['transcript_id'] = transcript.id
                            fout.write(str_gtf(transcript) + '\n')
                        else:
                            raise(NotImplementedError())
                        for exon in db.children(transcript.id):
                            if exon.featuretype in ['exon'] and not str(exon) in written_exons:
                            # write down the exon
                                if infmt == outfmt:
                                    fout.write(str(exon) + '\n')
                                elif infmt == 'gff' and outfmt == 'gtf':
                                    exon['gene_id'] = feature.id
                                    exon['transcript_id'] = transcript.id
                                    fout.write(str_gtf(exon) + '\n')
                                    written_exons = written_exons + [written_exons]

                                else:
                                    raise(NotImplementedError())

            
            ######### Write the gene to the file as is ############
            # CAVE: what to do here?
            elif not feature.id in extend_dictionary.keys(): # write the gene and all children as they are in the file:
                written_exons = []
                #if verbose:
                #    print("%s shoudn't be extended - omitting..." % feature.id)
                if outfmt == infmt:  
                    fout.write(str(feature) + '\n') # genes and transcripts are children of themselves
                    for transcript in db.children(feature.id,featuretype = ['transcript','mRNA']):
                        fout.write(str(transcript) + '\n') 
                        for exon in db.children(transcript.id):
                            if exon.featuretype in ['exon'] and not str(exon) in written_exons:
                            # write down the exon
                                fout.write(str(exon) + '\n')  
                                written_exons = written_exons + [str(exon)]
                    # If exons are found as children of a gene, write them down - to accommondate the cases of the same IDs in genes and transcripts
                    for exon in db.children(feature.id):
                        if exon.featuretype in ['exon'] and not str(exon) in written_exons:
                        # write down the exon
                            fout.write(str(exon) + '\n') 
                            written_exons = written_exons + [str(exon)]
                elif outfmt == 'gtf' and infmt == 'gff':
                    feature['gene_id'] = feature['ID']
                    # remove unnecessary atributes
                    fout.write(str_gtf(feature) + '\n') # genes and transcripts are children of themselves
                    for transcript in db.children(feature.id):
                        if transcript.featuretype in ['transcript','mRNA']:
                            transcript['gene_id'] = feature['gene_id']
                            transcript['transcript_id'] = transcript.id
                            fout.write(str_gtf(transcript) + '\n')
                            for exon in db.children(transcript.id):
                                exon['gene_id'] = transcript['gene_id']
                                exon['transcript_id'] = transcript['transcript_id']

                                if exon.featuretype in ['exon'] and not str(exon) in written_exons:
                                # write down the exon
                                    fout.write(str_gtf(exon) + '\n')  
                                    written_exons = written_exons + [str(exon)]
                    # If exons are found as children of a gene, write them down - to accommondate the cases of the same IDs in genes and transcripts
                    # CAVE - one needs to prevent duplicated exons from being written down!
                    for exon in db.children(feature.id):
                        exon['gene_id'] = feature.id
                        exon['transcript_id'] = feature.id
                        if exon.featuretype in ['exon'] and not str(exon) in written_exons:
                        # write down the exon
                            fout.write(str_gtf(exon) + '\n') 
                            written_exons = written_exons + [str(exon)]
                else:
                    raise(NotImplementedError())           
            else:
                continue 

# Gene extension function 
def extend_genes(genefile,peaksfile,outfile,maxdist,temp_dir,verbose,extension_type,infmt=None,outfmt=None,tag = None,clip_mode = None):
    # files with extensions:
    extension_table = temp_dir + '/extensions.tsv'
    #extension_bed = temp_dir + '/extensions.bed'

    if not outfmt in ['gtf','gff']:
        print("Unknown output format: %s" % outfile)
        quit()

    if verbose > 1:
        print('======== Gene extension: Loading annotations ===')
    #genes = check_ext_read_file(genefile,featuretype = 'gene')
    if infmt == 'bed':
        genes = parse_bed(genefile)
    elif infmt == 'gff':
        genes = parse_gff(genefile,'gene')
    else:
        genes = parse_gtf(genefile,'gene')


    if verbose > 1:
        print('\t%s genes loaded' % len(genes))
    peaks = parse_bed(peaksfile)

    if verbose > 1:
        print('\t%s peaks loaded' % len(peaks))
    
    if len(genes) == 0:
        print('No genes loaded! Please, make sure your .gff/gtf. file contains "gene" features!')
        quit()
    #peaks_d = {peak.id:peak for peak in peaks}
    #genes_d = {gene.id:gene for gene in genes}

    ####################### Assign peaks to genes, determine how much to extend #################################
    peaks_path = '%s/%s' % (temp_dir,'_peaks_tmp')
    genes_path = '%s/%s' % (temp_dir,'_genes_tmp')

    with open(peaks_path,'w') as file:
        for peak in peaks:
            file.write('\t'.join([peak.chrom,str(peak.start),str(peak.end),peak.id,'0',peak.strand]) + "\n")
    file.close()

    with open(genes_path,'w') as file:
        for gene in genes:
            file.write('\t'.join([gene.chrom,str(gene.start),str(gene.end),gene.id,'0',gene.strand]) + "\n")
    file.close()

    if verbose > 1:
        print("Written temporary files:\n\t%s\n\t%s" % (peaks_path,genes_path) )

    # use bedtools to identify the closest genes for every peak 
    os.system('bedtools sort -i %s/_peaks_tmp > %s/_peaks_tmp_sorted' % (temp_dir,temp_dir))
    os.system('bedtools sort -i %s/_genes_tmp > %s/_genes_tmp_sorted' % (temp_dir,temp_dir))
    #cmd = "bedtools closest -id -s -D a -a %s/_peaks_tmp_sorted -b %s/_genes_tmp_sorted  | cut -f 4,10,13  | awk '$3>=-%s'" % (temp_dir,temp_dir,maxdist)
    # awk 'BEGIN{OFS=@\\t@}{if($NF==0){if($6==@+@){$NF=-($3-$9)}else{$NF=-($8-$2)}};print $0}'
    cmd = "bedtools closest -id -s -D a -a %s/_peaks_tmp_sorted -b %s/_genes_tmp_sorted | awk '$NF!=-1' | awk '($6 == @+@ && $8<=$2)||($6==@-@ && $9 >= $3)' | awk 'BEGIN{OFS=@\\t@}{if($6==@+@){$NF=-($3-$9)}else{$NF=-($8-$2)};print $0}' | cut -f 4,10,13 | awk '$3>=-%s'" % (temp_dir,temp_dir,str(maxdist))
    cmd = cmd.replace('@','"')
    if verbose > 1:
        print('Running:\n\t\t%s > [output]' % cmd)
    out = os.popen(cmd).read()
    # write to a file:
    cmd = cmd + " > %s/_genes_peaks_closest" % (temp_dir)
    if verbose > 1:
        print('Running:\n\t\t%s' % cmd)
    os.system(cmd)
    cmd = "cp %s/_genes_peaks_closest %s/genes_peaks_closest.tsv" % (temp_dir,temp_dir)
    if verbose > 1:
        print('Running:\n\t\t%s' % cmd)
    os.system(cmd)

    # Assign genes to the most downstream peaks below extension threshold:
    peaks2genes = {x.split('\t')[0]:[x.split('\t')[1],int(x.split('\t')[2])] for x in out.split('\n')[:-1]}
    genes2peaks = {gene:[(k,v[1]) for k,v in peaks2genes.items() if gene in v] for gene in set([v[0] for v in peaks2genes.values()])}
    if verbose > 1:
        print('======== Asssigning genes to peaks =============')
    if verbose > 0:
        print('\t%s peaks assigned to %s genes.' % (len(peaks2genes),len(genes2peaks)))
        print('\tAverage n peaks per gene: %s' % round(np.mean([len(v) for v in genes2peaks.values()]),1))

    if verbose > 1:
        print('======== Gene extension statistics =============')
    # Select the most downstream peak per gene if not more than threshold 
    # TODO: simplify
    extend = {k:[x for x in v if abs(x[1])<maxdist] for k,v in genes2peaks.items()} 
    extend = {k:v for k,v in extend.items() if v}
    extend = {k:[x for x in v if x[1] == min([x[1] for x in v ])][0] for k,v in extend.items()}
    extend_dictionary = {k:abs(v[1]) for k,v in extend.items()}
    # this dictionary will be the input to the extension function
    #################################### Clip gene extensions overlapping other genes ####################################################
    # Check whether the extension intervenes into another gene
    if verbose > 1:
        print("Checking extensions ...")
    upd_counter = 0
    # for each gene extension, check if it's extension interferes with any of the genes on the same chromosome and strand.
    # Drop extensions for the genes already overlapping other genes. 
    for gene in genes:
        if gene.id in extend_dictionary.keys():
            # Screen for gene overlaps - remove the genes from extension if overlapping other genes regardless of the strand 
            if gene.strand == '+':
                overlapped = [x for x in genes if x.start <= gene.end and x.end  > gene.end and x.chrom == gene.chrom]
            elif gene.strand == '-':
                overlapped = [x for x in genes if x.end >= gene.start and x.end < gene.end and x.chrom == gene.chrom]
            if len(overlapped)>0:
                if verbose > 2:
                    print('%s overlaps the genes: %s' % (gene.id,','.join([x.id for x in overlapped])))
                    del extend_dictionary[gene.id]      
            else:
                # get the list of genes to consider for extension clipping:
                if clip_mode == 'sense':
                    genes_consider = [x for x in genes if x.chrom == gene.chrom and x.strand == gene.strand]
                elif clip_mode == 'both':
                    genes_consider = [x for x in genes if x.chrom == gene.chrom]
                else:
                    raise(ValueError('Unknown extension format'))
                #if verbose > 2:
                #    print('%s: %s genes on the same strand' % (gene.id,len(genes_same_strand)))
                # How to treat the genes that are also in the extension dictionary? 
                # find genes that would have been overlapped if any:
                if gene.strand == '+':
                    overlapped = [x for x in genes_consider if x.start <= gene.end+extend_dictionary[gene.id] and x.start > gene.end]
                    if len(overlapped)>0:
                        # if multiple overlaps, select the one 
                        new_ext = min([x.start for x in overlapped]) - gene.end -1 
                        if verbose > 2:
                            print('Extension check: found downstream gene overlap: %s --> %s' % (gene.id,','.join([x.id for x in overlapped])))
                            print('Updated extension: %s -> %s' % (extend_dictionary[gene.id],new_ext))
                        extend_dictionary[gene.id] = new_ext
                        upd_counter += 1
                else:
                    overlapped = [x for x in genes_consider if x.end >= gene.start-extend_dictionary[gene.id] and x.end < gene.start]
                    if len(overlapped)>0:
                        new_ext = gene.start-max([x.end for x in overlapped])-1
                        if verbose > 2:
                            print('Exension check: Found downstream gene overlap: %s --> %s' % (gene.id,','.join([x.id for x in overlapped])))
                            print('Updated extension: %s -> %s' % (extend_dictionary[gene.id],new_ext))
                        extend_dictionary[gene.id] = new_ext
                        upd_counter +=1
    if verbose > 1:
        print('Updated %s extensions (to prevent extensions into downstream genes).' % upd_counter)
    ########################################################################################
    exts = [abs(v[1]) for v in extend.values()]
    if verbose > 0:
        print('\tMaximal allowed extension: %s' % maxdist)
        print('\tMean: %s\n\tMedian: %s\n\tMax: %s' % (round(np.mean(exts),1),round(np.median(exts),1),round(np.max(exts),1)))
    with open(extension_table,'w') as file:
        for k,v in extend.items():
            file.write(k + '\t' + "\t".join([v[0],str(abs(v[1]))]) + '\n')
    file.close()

    if verbose > 1:
        print('\tGene extensions have been written to %s' % (extension_table))
    if verbose > 1:
        print('======== Writing output files ==================')

    ############### Finally, extend the genes ###########################
    if (outfmt == 'gtf' or outfmt == 'gff') and infmt == 'bed':
        raise(NotImplementedError('Can not output .%s file for input in .%s format!' % (outfmt,infmt)))
    elif outfmt == 'bed':
        raise(NotImplementedError())
    elif outfmt == 'gff':
        if verbose > 1:
            print('\tOutptut format - gff') 
        db = gffutils_import_gxf(genefile)
        extend_gff(db,extend_dictionary,outfile,extension_mode = extension_type,tag = tag,verbose = verbose,infmt = infmt,outfmt = outfmt)
    elif outfmt == 'gtf':
        if verbose > 1:
            print('\tOutptut format - gtf') 
        db = gffutils_import_gxf(genefile)
        extend_gff(db,extend_dictionary,outfile,extension_mode = extension_type,tag = tag,verbose = verbose, infmt = infmt,outfmt = outfmt)
    else:
        raise(ValueError('Unknown output format!'))
    if verbose > 1:
        print('\tExtended genes written: %s' % outfile)

# Report functions 
############################# Orphan peaks ############################################
# Orphan peaks 
def add_orphan_peaks(infile = None,peaksbed = None,fmt = None,tmp_outfile = None,tag = None,verbose = False):
    tag = 'GE_orphan'
    """This function takes orphan peaks and appends them to the input file"""
    if not fmt:
        raise(ValueError('Please, specify file format!'))
    if fmt in ['gff','gtf']:
        regs = parse_bed(peaksbed)
        if fmt == 'gff':
            with open(infile,'a') as file:
                for reg in regs:
                    gid = 'g.'+reg.id
                    tid = 't.'+reg.id
                    file.write('\t'.join([reg.chrom,tag,'gene',str(reg.start),str(reg.end),'.',reg.strand,'.','ID=%s' % gid])+'\n')
                    file.write('\t'.join([reg.chrom,tag,'transcript',str(reg.start),str(reg.end),'.',reg.strand,'.','ID=%s; Parent=%s' % (gid,tid)])+'\n')
                    file.write('\t'.join([reg.chrom,tag,'exon',str(reg.start),str(reg.end),'.',reg.strand,'.','ID=%s; Parent=%s' % (gid, tid)])+'\n')
        elif fmt == 'gtf':
            # write to the output file:
            with open(infile,'a') as file:
                for reg in regs:
                    gid = 'g.'+reg.id
                    tid = 't.'+reg.id
                    file.write('\t'.join([reg.chrom,tag,'gene',str(reg.start),str(reg.end),'.',reg.strand,'.','gene_id "%s"' % gid])+'\n')
                    file.write('\t'.join([reg.chrom,tag,'transcript',str(reg.start),str(reg.end),'.',reg.strand,'.','gene_id "%s"; transcript_id "%s"' % (gid,tid)])+'\n')
                    file.write('\t'.join([reg.chrom,tag,'exon',str(reg.start),str(reg.end),'.',reg.strand,'.','gene_id "%s"; transcript_id "%s"' % (gid, tid)])+'\n')
        if verbose > -1:
            print('%s orphan peaks added' % (str(len(regs))))



def merge_orphan_distance(orphan_bed = None,orphan_merged_bed = None,genic_bed = None,tempdir = None,maxsize = None,maxdist = None,verbose = False):
    """This function merges orphan peaks by distance"""
    pref = 'Orphan merging: '
    #maxsize = 100000 # maximum size of a peak cluster 
    #maxdist = 10000 # maximum distance of merging 
    # merge orphan peaks by distance
    cmd = "bedtools merge -s -i %s -c 4,5,6 -o distinct,max,distinct -d %s | awk 'NF==6' | grep  , | awk '$3-$2<%s' > %s/_orphan_merged.bed" % (orphan_bed,maxdist,maxsize,tempdir)
    if verbose > 1:
        print('Maximum distance: %s; Maximum cluster size: %s' % (maxdist,maxsize))
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    # Outersect the clusters with genes:
    cmd = "bedtools intersect -a %s/_orphan_merged.bed -b %s -v > %s/_orphan_merged_nov.bed; mv %s/_orphan_merged_nov.bed %s/_orphan_merged.bed" % (tempdir,genic_bed,tempdir,tempdir,tempdir)
    
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    # rename orphan clusters:
    cmd = "awk '{print  $4@\\tpeak.cl@NR}' %s/_orphan_merged.bed > %s/_peak_to_cluster" % (tempdir,tempdir)
    cmd = cmd.replace('@','"')
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    
    
    # select the ones to remove 
    cmd = "cut -f 4 %s/_orphan_merged.bed  | sed 's/,/\\n/g' | sort > %s/_orphan_toremove.txt" % (tempdir,tempdir)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    cmd = "grep -w -v -f %s/_orphan_toremove.txt %s | cut -f 1-6 > %s/_orphan_singleton.bed" % (tempdir,orphan_bed,tempdir)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    #cmd = "sed -i 's/_peak_//g'  %s/_orphan_merged.bed; sed -i 's/,/./g' %s/_orphan_merged.bed" % (tempdir,tempdir)
    #if verbose > 1:
    #    print('Running:\n\t%s' % cmd)
    #ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    cmd = "awk 'BEGIN{OFS=@\\t@}FNR==NR { p2c[$1]=$2; next }{print $1,$2,$3,p2c[$4],$5,$6}' %s/_peak_to_cluster %s/_orphan_merged.bed > %s/orphan_clusters.bed" % (tempdir,tempdir,tempdir)
    cmd = cmd.replace('@','"')
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    cmd = "cat %s/orphan_clusters.bed %s/_orphan_singleton.bed | awk 'NF==6' > %s" % (tempdir,tempdir,orphan_merged_bed)
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)


# get median length of the gene:
def get_median_gene_length(inputfile = None,fmt = None):
    if fmt == 'gtf':
        regs = parse_gtf(inputfile,featuretype = 'gene')
    elif fmt == 'gff':
        regs = parse_gff(inputfile,featuretype = 'gene')
    elif fmt == 'bed':
        regs = parse_bed(inputfile)
    else:
        print("Unknown input format!")
    med = np.median([x.end - x.start for x in regs])
    return(med)


def subsample_bam(inputbam = None,outputbam = None,nreads = None,verbose = True,threads = '1'):
    #cmd = "samtools idxstats -@ %s %s | cut -f 3 |  awk -v ct=%s 'BEGIN{total=0}{total+=$1}END{print ct/total}' | sed 's/,/./g'" % (str(threads),str(inputbam),str(nreads))
    cmd = "samtools view -@ %s -c -F 256 %s | awk -v ct=%s 'BEGIN{OFMT=toreplace1}{if(ct<=$1){print ct/$1}else{print toreplace2}}' | sed 's/,/./g'" % (str(threads),str(inputbam),str(nreads))
    cmd = cmd.replace('toreplace1','"%f"')
    cmd = cmd.replace('toreplace2','"1.0"')
    if verbose > 1:
        print('Running:\n\t%s' % cmd)
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    frac = ps.communicate()[0].decode("utf-8").rstrip()
    # subsample the bam using samtools view    
    cmd = "samtools view -@ %s -F 256 -h -b -s %s %s -o %s" % (str(threads),str(frac),inputbam,outputbam)
    if verbose> 1:
        print('Subsampling %s to %s reads => %s\n%s' % (inputbam,nreads,outputbam,cmd),flush = False)
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if verbose > 1:
        print(ps.communicate())
    else:
        ps.communicate()


###################### Mapping estimate #####################
# Estimate intergenic mapping for an new annotation  
# create a chromosome size file - should be another function
def get_chrsizes(tempdir = None,bamfile = None,outfile = None,verbose = False):
    cmd = "samtools idxstats %s | cut -f 1-2 | awk '$2!=0' > %s" % (bamfile,outfile)
    cmd = cmd.replace('__','"')
    if verbose > 1 :
        print('Running:\n\t%s' % cmd)
    ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

def get_genic_beds(genomeanno = None,genomechr = None,genicbed = None,intergenicbed = None,verbose = False,infmt = None):
    """Splits the genome into bed files with genic and intergenic regions"""

    # Genic regions
    if infmt in ['gtf','gff']:
        cmd = "awk -F '\\t|;' 'BEGIN{OFS=@\\t@} $3~/gene/ {print $1,$4,$5,@reg@NR,0,$7}' %s | bedtools sort -i - -g %s | bedtools merge -i - > %s" % (genomeanno,genomechr,genicbed)
        cmd = cmd.replace('@','"')
        if verbose > 1 :
            print('Running:\n\t%s' % cmd)
        ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    else:
        raise(NotImplementedError())
    ## Intergenic regions
    #cmd = "complementBed -i %s -g %s > %s" % (genicbed,genomechr,intergenicbed)
    #if verbose > 1 :
    #    print('Running:\n\t%s' % cmd)
    #ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

def estimate_mapping(bamfile=None,genicbed=None,intergenicbed=None,threads=1,verbose = False):
    if bamfile is None or genicbed is None or intergenicbed is None:
        print('Missing arguments for estimate mapping')
        quit()
    else:
        # Genic reads
        cmd = "samtools view -@ %s -c --region %s %s" % (threads,genicbed,bamfile)
        if verbose > 2 :
            print('Running:\n\t%s' % cmd)
        ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        Ngen = ps.communicate()[0].decode("utf-8").rstrip()
        
        # Intergenic reads 
        #cmd = "samtools view -@ %s -c --region %s %s" % (threads,intergenicbed,bamfile)
        #if verbose > 1 :
        #    print('Running:\n\t%s' % cmd)
        #ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        #Nigen = ps.communicate()[0].decode("utf-8").rstrip()

        # Total reads 
        cmd = "samtools view -@ %s -c  %s" % (threads,bamfile)
        if verbose > 2 :
            print('Running:\n\t%s\n' % cmd)
        ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        Ntot = ps.communicate()[0].decode("utf-8").rstrip()
        
        # Mapped reads 
        cmd = "samtools view -F 4 -@ %s -c  %s" % (threads,bamfile)
        if verbose > 2 :
            print('Running:\n\t%s\n' % cmd)
        ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        Nmap = ps.communicate()[0].decode("utf-8").rstrip()

        Ntot = int(Ntot)
        Nmap = int(Nmap)
        Ngen = int(Ngen)
        #Nigen = int(Nigen)
        Nigen = int(Nmap - Ngen)
        
def count_reads(bamfile=None,bed = None,flags = '',threads=1,verbose = False):
    # Genic reads
    if bed:
        cmd = "samtools view %s -@ %s -c --region %s %s" % (flags,threads,bed,bamfile)
    else:
        cmd = "samtools view %s -@ %s -c  %s" % (flags,threads,bamfile)
    if verbose > 2 :
        print('Running:\n\t%s' % cmd)
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    Nread = ps.communicate()[0].decode("utf-8").rstrip()
    return(int(Nread))


# add missing genes to the annotation file 
def get_featuretypes(infile = None):
    """To quickly check whether the file contains genes"""
    with open(infile) as file:
        return(set(line.split('\t')[2] for line in file if not '#' in line))

def add_gene_features(infile = None,outfile = None, infmt = None,verbose = False):
    """Add gene features based on transcripts"""
    # record the maximal span per gene 
    if verbose > 1:
        print('Loading the database ...',flush = '')
    db = gffutils.create_db(infile, ':memory:',disable_infer_genes=True,disable_infer_transcripts=True, merge_strategy = 'create_unique',transform =gffutils_transform_func,keep_order = True)
    if verbose > 1:
        print('done.')   
    # create a transcript2gene and gene2transcript dicts:
    t2g = {}
    for feature in db.features_of_type('transcript'):
        # CAVE: this is where parsing may fail
        if infmt == 'gtf':
            geneid = feature['gene_id'][0]
        elif infmt == 'gff':
            geneid = feature['Parent'][0]
        t2g.update({feature.id:geneid})
    g2t = {}
    for t,g in t2g.items():
        g2t.update({g:[t for t,v in t2g.items() if v == g]})
    with open(outfile,'w') as ofile:
        for i,feature in enumerate(db.features_of_type("transcript")):
            if verbose > 1 and i % 1000 == 0:
                print('%s/%s transcrpts done.' % (i,len(t2g)))
            if verbose > 2:
                print(feature.id)   
            # check if there is a gene_id 
            if infmt == 'gtf':
                gene = feature
                gene.featuretype = 'gene'
                geneid = t2g[gene.id]
                ofile.write(str_gtf(gene,attributes = ['gene_id']) + '\n') 
                #transcripts = [feature for feature in db.features_of_type('transcript') if geneid in feature['gene_id']]
                transcripts = [db[id] for id in g2t[geneid]]
                if verbose > 2:
                    print('%s transcripts found.' % len(transcripts))
                for transcript in transcripts:
                    if verbose > 2:
                        print(transcript.id)   
                    ofile.write(str(transcript) + '\n')
                    for child in db.children(db[transcript.id]):
                        ofile.write(str(child) + '\n')
            if infmt == 'gff':
                """Just copy as parent"""
                if 'Parent' in feature.attributes:
                    #print('Found Parent! -> setting as a gene ID')
                    gene = feature
                    gene.featuretype = 'gene'
                    gene['ID'] = gene['Parent'][0]
                    geneid = gene['Parent'][0]
                    ofile.write(str_gff(gene,attributes = ['ID']) + ';\n')
                    # write the gene:
                    #geneid = gene['Parent'][0]
                    #o = "\t".join([gene.chrom,gene.source,str(gene.start),str(gene.end),'.',gene.strand,'.','gene_id "' + geneid + '"'])
                    #ofile.write(o + '\n')         
                    transcripts = transcripts = [db[id] for id in g2t[geneid]]
                    #print('%s transcripts found.' % len(transcripts))
                    for transcript in transcripts:
                        ofile.write(str(transcript) + '\n')
                        for child in db.children(db[transcript.id]):
                            ofile.write(str(child) + '\n')
                else:
                    print('No parent found for transcript %s, where are the genes?' % feature.id)

def add_transcript_features(infile = None,outfile = None, infmt = None,verbose = False):
    # guess transcript features from the exons 
    raise(NotImplementedError())

    
def mRNA2transcript(infile = None,outfile = None,verbose = False):
    """In .gff/.gtf file, change the features of type 'mRNA' into 'transcript' """
    cmd = "awk -F '\\t' 'BEGIN{OFS=@\\t@}{if($3==@mRNA@){$3=@transcript@};print $0}' %s > %s" % (infile,outfile)
    cmd = cmd.replace('@','"')
    if verbose > 1:
        print('Running: %s' % cmd)
    os.system(cmd)


######################### 5' clipping ########################################
def gffutils_transform_func(x):
    # Specific function for gffuitls import
    if '' in x.attributes.keys():
        x.attributes.pop('')
    for k in x.attributes.keys():
        v = x.attributes[k]
        if len(v)>0:
            v = v[0]
        if '"' in v:
            v=v.replace('"','')
            x.attributes[k] = v
    return(x)
    

def check_overlap(a,b):
    # check 5' overlap for 2 features
    if a.strand != b.strand or a.chrom != b.chrom:
        return(False)
    else:
        if a.end < b.start:
            return(False)
        elif a.start > b.end:
            return(False)
        else:
            if a.start < b.start and a.end < b.end:
                if a.strand == '+':
                    return('3')
                else:
                    return('5')
            elif a.start <= b.start and a.end > b.end:
                return('full_a')
            elif a.start >= b.start and a.end <= b.end:
                return('full_b')
            elif a.start > b.start and a.end > b.end:
                if a.strand == '+':
                    return('5')
                else:
                    return('3')
            else:
                return(False)


def clip5_process_gene(gene,genes,db,verbose = False,tag = '_5clip'):
    logstr = []
    overlapped_genes = [x for x in genes if check_overlap(gene, x) == "5"]
    if len(overlapped_genes) > 0:
        overlapped_genes = [x for x in genes if check_overlap(gene, x) == "5"]
        ovgene = db[overlapped_genes[0].id]
        if verbose > 2:
            print("5' clipping: found 5' overlap %s -> %s" % (ovgene.id, gene.id))
        # the gene should be clipped
        old_start = gene.start
        old_end = gene.end
        new_start = gene.start
        new_end = gene.end
        if gene.strand == "+":
            new_start = ovgene.end + 1
            gene.start = new_start
        else:
            new_end = ovgene.start - 1
            gene.end = new_end
        gene_range = [gene.start, gene.end]
        gene.source = gene.source + tag
        outstr = str(gene) + '\n'
        # gene - overlapping_gene - gene_strand - overlapping_gene_strand - old_start - new_start - old_end - new_end 
        logstr = "\t".join([gene.id,ovgene.id,gene.strand,ovgene.strand,str(old_start),str(new_start),str(old_end),str(new_end)]) + '\n'
        for child in db.children(db[gene.id]):
            if not child.featuretype == 'gene':
                # remove child features outside of the new range:
                if (not child.end < gene_range[0]) or (not child.start > gene_range[1]): 
                    overlap_type = check_overlap(child,ovgene)
                    if not overlap_type == 'full_b':   
                        if gene.strand == '+' and overlap_type == '5':
                            # modify start
                            child.start = new_start
                        elif gene.strand == '-' and overlap_type == '5':
                                child.end = new_end
                        child.source = child.source + tag
                        outstr = outstr + str(child)+'\n'
                    else:
                        if verbose > 2:
                            print('Child feature is fully contained within a downstream gene - OMITTING:\n%s\n%s' % (str(child),str(ovgene)))
                        pass
                else:
                    if verbose > 2:
                        print('Child is outside of the new gene range\n%s\n%s' % (str(child),str(gene)))
                    pass
        return(outstr,logstr)
    else:
        outstr = str(gene) + '\n'
        for child in db.children(db[gene.id]):
            if not child.featuretype == 'gene':
                outstr = outstr + str(child) + '\n'
        return(outstr,logstr)


def clip5_worker_process(genes, infile, i, results,logs,verbose,tag):
    db = gffutils.create_db(
        infile,
        ":memory:",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        merge_strategy="create_unique",
        transform=gffutils_transform_func
    )
    all_genes = [x for x in db.features_of_type("gene")]
    
    #if verbose > 2:
        #print('process %s: total number of genes %s' % (i,len(all_genes)))
        #print('process %s: number of genes in chunk %s' % (i,len(genes)))
    result = []
    log = []
    for gene in genes:
        result_out,log_out = clip5_process_gene(gene, all_genes, db,verbose,tag)
        result = result + [result_out]
        if len(log_out)>0:
            log = log + [log_out]
    # update the results vector
    results[i] = result
    # update logs vector 
    logs[i] = log


def clip_5_overlaps(infile = None,outfile = None,threads = 1,verbose = False,tag = '_5clip'):
    logfile = outfile + '5clip.log'
    import math
    # Load the database
    if verbose > 1:
        print('Loading gene database.')
    db = gffutils.create_db(
        infile,
        ":memory:",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        merge_strategy="create_unique",
        transform=gffutils_transform_func
    )
    genes = [x for x in db.features_of_type("gene")]
    print("%s genes loaded." % len(genes))
    # Split the genes into chunks for each worker process
    chunk_size = math.floor(len(genes)/threads)
    chunks = [genes[i*(chunk_size):(i+1)*chunk_size] for i in range(0,threads-1)]
    chunks = chunks + [genes[sum([len(x) for x in chunks]):]]  
    # Create a list to hold the results from each worker process
    manager = multiprocessing.Manager()
    results = manager.list([[] for _ in range(threads)])
    logs = manager.list([[] for _ in range(threads)])

    # Start the worker processes
    processes = []
    for i in range(threads):
        p = multiprocessing.Process(target=clip5_worker_process, args=(chunks[i], infile, i, results,logs,verbose,tag))
        p.start()
        processes.append(p)

    
    # Wait for the worker processes to finish
    for p in processes:
        p.join()

    # Combine the results from all the worker processes
    final_results = results

    # Write the final results to the output file
    with open(outfile, "w") as outf:
        for result in final_results:
            outf.write(''.join(result))
    with open(logfile, "w") as outf:
        outf.write("\t".join(['gene_id','ov_gene_id','gene_strand','ov_gene_strand','old_start','new_start','old_end','new_end']) + '\n')
        for log in logs:
            if len(log)>0:
                outf.write(''.join(log))