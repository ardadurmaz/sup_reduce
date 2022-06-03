library(data.table)


setwd('~/wgs_ecoli/')

var_ft <- as.data.frame(fread('data/split_gatk_ecoli_calls.table', header=TRUE, sep='\t', stringsAsFactors = FALSE))
gt_idx <- c(2,3,4,5,which(grepl(pattern='\\.GT', colnames(var_ft))))
mut_mat <- sapply(which(grepl(pattern='\\.GT', colnames(var_ft))), function(i){
  ifelse(var_ft[,i] == var_ft$REF , 0, 1)
})
colnames(mut_mat) <- colnames(var_ft)[which(grepl(pattern='\\.GT', colnames(var_ft)))]
rownames(mut_mat) <- paste(var_ft$CHROM, var_ft$POS, var_ft$REF, var_ft$ALT, sep='_')

vaf_mat <- sapply(which(grepl(pattern='\\.AD', colnames(var_ft))), function(i){
  local_ad <- sapply(stringi::stri_split_fixed(pattern = ',', str=var_ft[,i]), function(s){
    temp <- as.numeric(s)
    return(temp[2]/sum(temp))
  })
  return(local_ad)
})
vaf_mat[is.na(vaf_mat)] <- 0
colnames(vaf_mat) <- colnames(var_ft)[which(grepl(pattern='\\.AD', colnames(var_ft)))]
rownames(vaf_mat) <- paste(var_ft$CHROM, var_ft$POS, var_ft$REF, var_ft$ALT, sep='_')


# Map Column Names to Samples #
name_map <- setNames(paste(paste0('X', c(10, 11, 12, 12, 12, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9)), paste0('p', c(10, 10, 10, 2, 8, 10, 10, 5, 8, 10, 10, 10, 10, 10, 10, 10)), sep=''),
                     paste0('S', 2:17))

colnames(vaf_mat) <- name_map[gsub(pattern='\\.AD', replacement = '', colnames(vaf_mat))]
colnames(mut_mat) <- name_map[gsub(pattern='\\.GT', replacement = '', colnames(mut_mat))]

mic_ft <- read.table('data/mic_ft.txt', header=TRUE, sep='\t')
mic_ft <- mic_ft[,c(1, which(grepl(pattern='_mle', colnames(mic_ft))))]
row_ids <- mic_ft$X
mic_ft <- data.matrix(mic_ft[,-1])
rownames(mic_ft) <- row_ids

common_ids <- intersect(colnames(vaf_mat), paste0(rownames(mic_ft), 'p10'))
vaf_mat <- vaf_mat[,match(common_ids, table=colnames(vaf_mat))]
mut_mat <- mut_mat[,match(common_ids, table=colnames(mut_mat))]
mic_mat <- t(mic_ft[match(common_ids, table=paste0(rownames(mic_ft), 'p10')),])
row_ids <- sort(union(rownames(vaf_mat), rownames(mut_mat)))
vaf_mat <- vaf_mat[match(row_ids, table=rownames(vaf_mat)),]
mut_mat <- mut_mat[match(row_ids, table=rownames(mut_mat)),]

require(dplyr)
gff_ft <- rtracklayer::readGFF('~/resources/GCF_000019425.1_ASM1942v1_genomic.gtf') %>%
  filter(source == 'RefSeq')

gene_idx <- sapply(1:length(var_ft$POS), function(x){
  temp <- which(gff_ft$seqid == var_ft$CHROM[x] & gff_ft$start <= var_ft$POS[x] & gff_ft$end >= var_ft$POS[x])
  if(length(temp) == 0){
    return(NA)
  }else{
    return(temp[1])
  }
  })

var_ft$gene <- gff_ft$gene[gene_idx]

# DO CCA #
## Do PCA PreProcessing ##
require(FactoMineR)
mca_res <- FactoMineR::MCA(data.frame(ifelse(t(mut_mat) == 1, 'Mut', 'NonMut')), ncp = 5, graph = FALSE)
factoextra::fviz_screeplot(mca_res, addlabels = TRUE, ylim = c(0, 45))
mut_lat <- mca_res$ind$coord
svd_res <- svd(scale(t(mic_mat)), nu = 3, nv = 3)
mic_lat <- scale(t(mic_mat)) %*% svd_res$v

fwrite(mut_lat, file='data/MutationLatentCoordinates.tsv', col.names = FALSE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')
fwrite(mic_lat, file='data/MICLatentCoordinates.tsv', col.names = FALSE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')
fwrite(scale(t(mic_mat)), file='data/MICMat.tsv', col.names = FALSE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')

fwrite(data.table('id' = rownames(mic_lat)), file='data/ObsAnnot.tsv', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')
fwrite(data.table('id' = rownames(mic_mat)), file='data/FeatAnnot.tsv', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = '\t')
