#--------------------------------------------------------------------------------
# Libraries
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library("optparse"))
quiet(library(tidyverse))
#--------------------------------------------------------------------------------
option_list = list(
    make_option(c("-m", "--meta"), type="character", default=NULL, help="meta mesoneo"),
    make_option(c("-s", "--samples"), type="character", default=NULL, help="samples vcf order"),
    make_option(c("-o", "--out"), type="character", default="output.txt", help="Output metadata all name"),
    make_option(c("-t", "--training"), type="character", default="output.txt", help="Output unsupervised file name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##----------------------------------------------------------------------------------------------
sinfo <- opt$meta
samIds = opt$samples  
outfile1 = opt$out
outfile2 = opt$training

# ------------- test
# samIds="samples/MESOAtSc_vcf_order.txt"
# sinfo="/home/alba/mesoneo/neo.impute.1000g.sampleInfo_clusterInfo.txt"
# outfile2="samples/MESOAtSc_supervised_ks.txt"
# outfile3="samples/MESOAtSc_unsupervised.txt"


readIn = function(samids,sInfo){
    ids = read_tsv(samids, col_names=F,show_col_types = FALSE)   
    sampleInfo = read_tsv(sInfo,show_col_types = FALSE)
    
    ids$sampleId = ids$X1
    ids$popId = sampleInfo$popId[match(ids$sampleId, sampleInfo$sampleId)]
    ids$clusterIBDFine = sampleInfo$clusterIBDFine[match(ids$sampleId, sampleInfo$sampleId)]
    ids$region = sampleInfo$region[match(ids$sampleId, sampleInfo$sampleId)]
    ids$color = sampleInfo$color[match(ids$sampleId, sampleInfo$sampleId)]
    ids$fill = sampleInfo$fill[match(ids$sampleId, sampleInfo$sampleId)]
    ids$shape = sampleInfo$shape[match(ids$sampleId, sampleInfo$sampleId)]
    ids$clusterAlias = sampleInfo$clusterAlias[match(ids$sampleId, sampleInfo$sampleId)]
    ids$ageAverage = sampleInfo$ageAverage[match(ids$sampleId, sampleInfo$sampleId)]
    ids$groupLabel = sampleInfo$groupLabel[match(ids$sampleId, sampleInfo$sampleId)]
    ids$country = sampleInfo$country[match(ids$sampleId, sampleInfo$sampleId)]
    ids$site = sampleInfo$site[match(ids$sampleId, sampleInfo$sampleId)]

    
    ids[is.na(ids$popId),]$popId <- "Faroese_modern"
    ids[ids$popId == "Faroese_modern",]$sampleId = apply(ids[ids$popId == "Faroese_modern",],1, function(x) strsplit(x, "-")[[1]][3])
    ids[ids$popId =="Faroese_modern",]$color <- "black"
    ids[ids$popId =="Faroese_modern",]$fill <- "black"
    ids[ids$popId =="Faroese_modern",]$shape <- 3
    ids[ids$popId =="Faroese_modern",]$region <- "NorthernEurope"
    ids[ids$popId =="Faroese_modern",]$groupLabel <- "Faroese_modern"
    ids[ids$popId =="Faroese_modern",]$clusterAlias <- "Faroese_modern"
    ids[ids$popId =="Faroese_modern",]$clusterIBDFine <- "Faroe_modern"
    ids[ids$popId =="Faroese_modern",]$ageAverage <- 0
    ids[ids$popId =="Faroese_modern",]$country <- "Faroe_Islands"
    ids[ids$popId =="Faroese_modern",]$site <- "Faroes"
    return(ids[,-1])

}

df = readIn(samIds,sinfo)
                                                         
df[df$groupLabel == "Faroe_Historical",]$clusterAlias <- "Faroese"
                                                         
df = df %>% mutate(training = ifelse(clusterAlias %in% c("Faroese","Faroese_modern"), 0, 1))


print("printing meta and training files...")
                                                         
write_tsv(df, outfile1)
write_tsv(df[,which(names(df) == "training")], outfile2, col_names=F)

                                                         
#df[df$training == 1,]$supervised <- dense_rank(df$clusterAlias[df$training == 1])
#write_tsv(df[,which(names(df) == "supervised")], outfile3, col_names=F)
