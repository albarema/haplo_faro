#--------------------------------------------------------------------------------
# Libraries
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(RColorBrewer))
quiet(library(optparse))


#--------------------------------------------------------------------------------
option_list = list(
  make_option(c("-p", "--panel"), type="character", default=NULL, help="panel name"),
  make_option(c("-b", "--best"), type="character", default=NULL, help="best runs"),
  make_option(c("-w", "--windowS"), type="character", default=NULL, help="window size"),
  make_option(c("-v", "--vcf"), type="character", default="output.txt", help="vcf order samples"),
  make_option(c("-m", "--meta"), type="character", default="output.txt", help="meta file name"),
  make_option(c("-k", "--ks"), type="character", default=3, help="K best"),
  make_option(c("-o", "--outfile"), type="character", default="output.txt", help="Output supervised file name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##----------------------------------------------------------------------------------------------
panel <- opt$panel
windowSize = opt$windowS  
pathIds = opt$vcf
sInfo = opt$meta
ks = opt$ks
allbest = opt$best
outfile = opt$outfile

#panelvcf ="MESOAtSc" 
#infile="EUR/haplonet.admixture.k2.s979.q"
#sInfo = "MESOAtSc.sampleInfo_clusterInfo.txt" #"neo.impute.1000g.sampleInfo_clusterInfo.txt"
#pathIds = paste(panelvcf,"_vcf_order.txt", sep="") #paste("/home/ssd-alba/samples/",panelvcf,"_samples_vcf_order.txt", sep="")   
#windowSize="wl"
#panel = "ancient.AtSc" # "ALL"
#--------------------------------------------------------------------------------

readIn = function(infile, samids,ks, samInfo){
  ids = read_tsv(samids, col_names=F,show_col_types = FALSE)   
  q <- read.table(infile) %>% mutate(K = ks)
  qf = cbind(ids, q) %>% as_tibble()
  allq = merge(samInfo, qf, by=1) %>% as_tibble()
  colnames(allq) = c("sampleId","popId","clusterIBDFine","region","color","fill","shape","clusterAlias", "ageAverage","groupLabel","country","site","training",paste("Q", 1:ks, sep=""), "K")
  ftab = gather(allq,"Q", "value", -c(sampleId, popId, clusterIBDFine, region, color, fill, shape, clusterAlias, ageAverage, groupLabel, K, country, site, training))
  return(ftab)
}
#--------------------------------------------------------------------------------

ld = read_tsv(allbest, show_col_types = FALSE)# contains the k and the best seed for each (all best file)

samsinfo = read_tsv(sInfo,show_col_types = FALSE)
samsinfo[samsinfo$popId == "Faroese_modern",]$sampleId <- paste("SUBJ-FOFG-",samsinfo[samsinfo$popId == "Faroese_modern",]$sampleId,sep="")

filepa = paste(panel,"/haplonet-unsuperv-",windowSize,".admixture.k", ld[ld$k == ks,]$k,".s", ld[ld$k == ks,]$seed,".q", sep="")


ftab = readIn(filepa, pathIds, ks, samsinfo)

# it is a faroese viking 
ftab[ftab$sampleId == "VK24",]$training = 0 
# 
ftab1 = ftab %>% filter(training != 0,K ==ks, Q == "Q1", value >=.95) %>% mutate(supervised = 1)
ftab2 =  ftab %>% filter(training != 0,K ==ks, Q == "Q2", value >=.95) %>% mutate(supervised = 2)
ftab3 =  ftab %>% filter(training != 0,K ==ks, Q == "Q3", value >=.95) %>% mutate(supervised = 3)

fnew = rbind(ftab1, ftab2, ftab3)

samsinfo$supervised <- fnew[match(samsinfo$sampleId, fnew$sampleId),]$supervised

samsinfo[is.na(samsinfo$supervised),]$supervised <- 0 

write_out = samsinfo %>% select(supervised) 

write_tsv(write_out, outfile, col_names=F)

