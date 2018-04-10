library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts <- read.table(args[1], header=TRUE)
samples <- read.table('DESeqDesign.txt', header=TRUE)

fbfemale <- DESeqDataSetFromMatrix(countData=counts[ , 1:4], colData=samples[1:4, ], design = ~AlignsTo)
fbfemale <- DESeq(fbfemale)
res <- results(fbfemale)
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/fbfemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

fbmale <- DESeqDataSetFromMatrix(countData=counts[ , 5:8], colData=samples[5:8, ], design = ~AlignsTo)
fbmale <- DESeq(fbmale)
res <- results(fbmale)
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/fbmale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

oefemale <- DESeqDataSetFromMatrix(countData=counts[ , 9:12], colData=samples[9:12, ], design = ~AlignsTo)
oefemale <- DESeq(oefemale)
res <- results(oefemale)
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/oefemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

oemale <- DESeqDataSetFromMatrix(countData=counts[ , 13:16], colData=samples[13:16, ], design = ~AlignsTo)
oemale <- DESeq(oemale)
res <- results(oemale)
res$gene <- rownames(res)
res <- res[, c(7,1,2,3,4,5,6)]
write.table(res, paste(dirname(args[1]), '/combined/oemale_deseq.tsv', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

