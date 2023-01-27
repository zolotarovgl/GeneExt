# what to inlcude in the report?
# Number of genes in the annotation 
# Number of peaks 
# Number the genes with peaks assigned 
# Gene extension distribution   
#############################################################


#maxdist = 10000
#quant = .25
#closest_gene_file = 'tmp/_genes_peaks_closest'
#allpeaks_cov_file = 'tmp/allpeaks_coverage.bed'
#allpeaks_noov_file = 'tmp/allpeaks_noov.bed'
#extension_file = 'tmp/extensions.tsv'


# script.R maxdist quant closest_gene_file allpeaks_cov_file allpeaks_noov_file extension_file
args <- commandArgs(trailingOnly = TRUE)


#args = c('10000','.25','tmp/genes_peaks_closest','allpeaks_coverage.bed','tmp/allpeaks_noov.bed','tmp/extensions.tsv',1)
maxdist = as.integer(args[1])
quant = as.numeric(args[2])
closest_gene_file = args[3]
allpeaks_cov_file = args[4]
allpeaks_noov_file = args[5]
extension_file = args[6]
verbosity = as.integer(args[7])

if(verbosity > 0){
    print(args)
}

pdf('report.pdf',width = 8.27, height = 11.69 )
par(mfrow = c(3,2))

############### read coverage and peak length ###############
cov = read.table(allpeaks_cov_file,fill = T,nrows = -1)
covn = read.table(allpeaks_noov_file,fill = T,nrows = -1)
al = cov$V7
nov = covn$V7
ov = cov$V7[!cov$V4 %in% covn$V4] 


############## Numerical summary ##################
library(stringr)
n_genes = 20000
n_peaks = length(al)
n_peaks_in_genes = length(ov)
n_peaks_noov = length(nov)

e = read.table(extension_file)
n_ext = length(unique(e$V1))

plot.new()
txt = sprintf('Number of peaks: %s\nNumber of peaks in genes: %s\nNumber of peaks outside of genes: %s\nN genes extended: %s',
n_peaks,n_peaks_in_genes,n_peaks_noov,n_ext)
text(x=0.1, y=0.9, txt,pos = 4,cex = 1)  # first 2 numbers are xy-coordinates within [0, 1]


# Gene extension density plot
plot(density(e$V3),main = 'Gene extension',xlab = 'Distance from TES, bp')
######## Distance distributions ############
d = read.table(closest_gene_file)
#d = tapply(d$V2,d$V3,max)
# plot distances of 3 closest peaks per gene 
d$V3 = -d$V3
d = d[d$V2 !='.',]
dd = lapply(split(d$V3,d$V2),FUN = function(x) sort(x,decreasing = F)[1])
dd = unlist(dd)
plot(density(dd/1000),main = paste0('Distance to the nearest non-genic peak\nMedian: ',round(median(dd/1000),2),' Kb; Mean: ',round(mean(dd/1000),2),' Kb'),
xlab = 'Distance to the nearest non-genic peak, Kb')
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
garbage <- dev.off()
#####################################################

# Report numerical values:
par(mfrow = c(3,2))
