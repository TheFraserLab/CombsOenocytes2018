library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts <- read.table(args[1], header=TRUE)
samples <- read.table('DESeqDesign.txt', header=TRUE)

# Female Fat Bodies

fbfemale <- DESeqDataSetFromMatrix(countData=counts[ , 1:4], colData=samples[1:4, ], design = ~Replicate + AlignsTo)
fbfemale <- DESeq(fbfemale)
res <- results(fbfemale, contrast=c('AlignsTo', 'ref', 'alt'))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/fbfemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

# Male Fat Bodies

fbmale <- DESeqDataSetFromMatrix(countData=counts[ , 5:8], colData=samples[5:8, ], design = ~Replicate + AlignsTo)
fbmale <- DESeq(fbmale)
res <- results(fbmale, contrast=c('AlignsTo', 'ref', 'alt'))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/fbmale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

# Female Oenocytes

oefemale <- DESeqDataSetFromMatrix(countData=counts[ , 9:12], colData=samples[9:12, ], design = ~Replicate + AlignsTo)
oefemale <- DESeq(oefemale)
res <- results(oefemale, contrast=c('AlignsTo', 'ref', 'alt'))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/oefemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

# Male oenocytes

oemale <- DESeqDataSetFromMatrix(countData=counts[ , 13:16], colData=samples[13:16, ], design = ~Replicate + AlignsTo)
oemale <- DESeq(oemale)
res <- results(oemale, contrast=c('AlignsTo', 'ref', 'alt'))
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/oemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

