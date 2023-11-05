import gffutils
infile = 'tmp_test1/sample_genome.fixed.gff'
outfile = 'boo.gff'
verbose = 3
infmt = 'gff'
outfmt = 'gff'

# GXF helper 
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
##########################
db = gffutils_import_gxf(infile)
# Create a dictionary to store the longest transcript for each gene
g2tid = {}
g2t = {}
# Iterate over all genes in the database
cnt = 1
for gene in db.features_of_type('gene'):
    longest_transcript = None
    transcripts = [x for x in db.children(gene, featuretype='transcript')]
    if len(transcripts)>0:
        lengths = {x.id:x.end-x.start for x in transcripts}
        t2model= {x.id:x for x in transcripts}
        g2tid[gene.id] = max(lengths, key=lengths.get)
        g2t[gene.id] = t2model[max(lengths, key=lengths.get)]
    cnt += 1 
if verbose:
    print('%s genes - %s transcripts' % (cnt,len(g2t)))

with open(outfile,'w') as ofile:
        for i,gene in enumerate(db.features_of_type("gene")):
            if gene.id in g2tid:
                # check if there is a gene_id 
                if infmt == 'gtf':
                    ofile.write(str_gtf(gene,attributes = ['gene_id']) + '\n') 
                    transcript = g2t[gene.id]
                    transcript_id = transcript.id
                    ofile.write(str(transcript) + '\n')
                    for child in db.children(transcript_id):
                        ofile.write(str(child)+ '\n')

                if infmt == 'gff':
                    ofile.write(str_gff(gene,attributes = ['ID']) + ';\n')
                    # write the gene:
                    #o = "\t".join([gene.chrom,gene.source,str(gene.start),str(gene.end),'.',gene.strand,'.','gene_id "' + gene.id + '"'])
                    #ofile.write(o + '\n')         
                    transcript = g2t[gene.id]
                    transcript_id = transcript.id
                    ofile.write(str(transcript) + '\n')
                    for child in db.children(transcript_id):
                        ofile.write(str(child)+ '\n')


