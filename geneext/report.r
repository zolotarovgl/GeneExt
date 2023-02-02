# what to inlcude in the report?
# Number of genes in the annotation 
# Number of peaks 
# Number the genes with peaks assigned 
# Gene extension distribution   
#############################################################
library(stringr)

#maxdist = 10000
#quant = .25
#closest_gene_file = 'tmp/_genes_peaks_closest'
#allpeaks_cov_file = 'tmp/allpeaks_coverage.bed'
#allpeaks_noov_file = 'tmp/allpeaks_noov.bed'
#extension_file = 'tmp/extensions.tsv'


# script.R maxdist quant closest_gene_file allpeaks_cov_file allpeaks_noov_file extension_file
args <- commandArgs(trailingOnly = TRUE)


#args = c('10000','.25','tmp/genes_peaks_closest','tmp/allpeaks_coverage.bed','tmp/allpeaks_noov.bed','tmp/extensions.tsv',1,'dummy.gtf')
maxdist = as.integer(args[1])
quant = as.numeric(args[2])
closest_gene_file = args[3]
allpeaks_cov_file = args[4]
allpeaks_noov_file = args[5]
extension_file = args[6]
verbosity = as.integer(args[7])
pref = args[8]
pref = unlist(str_split(pref,'\\.'))[length(unlist(str_split(pref,'\\.')))-1]
# read command call:
tempdir = 'tmp'
callcmd = readLines(sprintf('%s/_callcmd',tempdir))

if(verbosity > 1){
    print(args)
}

pdf(sprintf('%s.pdf',pref),width = 8.27, height = 11.69 )
par(mfrow = c(1,1))

############### read coverage and peak length ###############
cov = read.table(allpeaks_cov_file,fill = T,nrows = -1)
covn = read.table(allpeaks_noov_file,fill = T,nrows = -1)
al = cov$V7
nov = covn$V7
ov = cov$V7[!cov$V4 %in% covn$V4] 

# closest upstream genes per peak 
clo = read.table(closest_gene_file)
clo$V3 = -clo$V3
clo = clo[clo$V2 !='.',]



############## Numerical summary ##################
library(stringr)
n_genes = 20000
n_peaks = length(al)
n_peaks_in_genes = length(ov)
n_peaks_noov = length(nov)

e = read.table(extension_file)
n_ext = length(unique(e$V1))

# after considering threshold:
f = clo[clo$V3<=maxdist,]
peaks_below_threshold = length(unique(f$V1))
genes_below_threshold = length(unique(f$V2))
mean_peaks_per_gene_below_threshold = round(mean(sapply(split(f$V1,f$V2),length)),1)
#########################################################
plot.new()
txt = sprintf('%s\nNumber of peaks: %s\nNumber of peaks in genes: %s\nNumber of peaks outside of genes: %s\n-----After filtering (max %s bp)-----:\nN genes: %s\nN peaks: %s\nMean N peaks per gene: %s\nN genes extended: %s',
callcmd,n_peaks,n_peaks_in_genes,n_peaks_noov,maxdist,genes_below_threshold,peaks_below_threshold,mean_peaks_per_gene_below_threshold,n_ext)
text(x=0, y=0.9, txt,pos = 4,cex = 0.5)  # first 2 numbers are xy-coordinates within [0, 1]


######## Distance distributions ############
#d = tapply(d$V2,d$V3,max)
# plot distances of 3 closest peaks per gene 
par(mfrow = c(3,2))
thr = 50000 # cut distribution here
dd = lapply(split(clo$V3,clo$V2),FUN = function(x) sort(x,decreasing = F)[1])
dd = unlist(dd)
dd = dd[dd<=thr]
plot(density(dd/1000),main = paste0('Distance to the nearest non-genic peak\nMedian: ',round(median(dd/1000),2),' Kb; Mean: ',round(mean(dd/1000),2),' Kb'),
xlab = sprintf('Distance to the nearest non-genic peak, Kb ( cut at %s Kb)',thr/1000))
abline(v = maxdist/1000,lty = 2,lwd =2,col = 'red')

########### Peak coverage plots ################
plot(density(log10(al)),main = 'Peak coverage',xlab = 'log10 peak coverage',ylim = c(0,1))
lines(density(log10(nov)),col = 'red')
lines(density(log10(ov)),col = 'blue')
# does it make sense? the peaks not overlapping the genes will get lower coverage 
legend(x = 3,y = 0.8,fill = c('black','red','blue'),legend = c('All peaks','Not overlapping genes','Overlapping genes'),cex = 1)
cutoff = quantile(log10(ov),quant)
abline(v = cutoff,lty = 2,col = 'red')
# What about peak width? 

al = cov$V3-cov$V2
nov = covn$V3-covn$V2
ov = (cov$V3-cov$V2)[!cov$V4 %in% covn$V4] 
plot(density(log10(al)),main = 'Peak length',xlab = 'log10 peak length',ylim = c(0,1))
lines(density(log10(nov)),col = 'red')
lines(density(log10(ov)),col = 'blue')
# does it make sense? the peaks not overlapping the genes will get lower coverage 
legend(x = 3,y = 0.8,fill = c('black','red','blue'),legend = c('All peaks','Not overlapping genes','Overlapping genes'),cex = 1)
cutoff = quantile(log10(ov),.25)
abline(v = cutoff,lty = 2,col = 'red')

al = cov$V7/(cov$V3-cov$V2)
nov = covn$V7/(covn$V3-covn$V2)
ov = cov$V7[!cov$V4 %in% covn$V4]/(cov$V3[!cov$V4 %in% covn$V4] - cov$V2[!cov$V4 %in% covn$V4] ) 
plot(density(log10(al)),main = 'Peak coverage / Peak width',xlab = 'log10 peak coverage / width',ylim = c(0,1))
lines(density(log10(nov)),col = 'red')
lines(density(log10(ov)),col = 'blue')
# does it make sense? the peaks not overlapping the genes will get lower coverage 
legend(x = 2,y = 0.8,fill = c('black','red','blue'),legend = c('All peaks','Not overlapping genes','Overlapping genes'),cex = 1)
cutoff = quantile(log10(ov),quant)
abline(v = cutoff,lty = 2,col = 'red')

# Gene extension density plot
plot(density(e$V3),main = 'Gene extension length',xlab = 'Distance from TES, bp')

garbage <- dev.off()
#####################################################

