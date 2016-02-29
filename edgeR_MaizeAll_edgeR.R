rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_GF/HtSeq_ALLsamples/")
getwd()
library("edgeR")

####################################################################################################################################
# EDGER (using htseq-counts previously Uniquely mapped reads only)
####################################################################################################################################

# Import (to select prot. coding genes only) 
#gene_anno = read.table("TAIR10.gene.list.txt",sep="\t")

# Import the 14 htseq-count tables to start Diff.Expr.Analysis
B73_wt_06h.B1 <- read.table("B73_wt_06h.B1_counts.txt",header=F )
B73_wt_06h.B2 <- read.table("B73_wt_06h.B2_counts.txt",header=F )
B73_mt_06h.B1<- read.table("B73_mt_06h.B1_counts.txt",header=F )
B73_mt_06h.B2 <- read.table("B73_mt_06h.B2_counts.txt",header=F )
B73_wt_22h.B1 <- read.table("B73_wt_22h.B1_counts.txt",header=F )
B73_wt_22h.B2 <- read.table("B73_wt_22h.B2_counts.txt",header=F )
B73_mt_22h.B1<- read.table("B73_mt_22h.B1_counts.txt",header=F )
B73_mt_22h.B2 <- read.table("B73_mt_22h.B2_counts.txt",header=F )
B73_wt_48h.B1 <- read.table("B73_wt_48h.B1_counts.txt",header=F )
B73_wt_48h.B2 <- read.table("B73_wt_48h.B2_counts.txt",header=F )
B73_mt_48h.B1<- read.table("B73_mt_48h.B1_counts.txt",header=F )
B73_mt_48h.B2 <- read.table("B73_mt_48h.B2_counts.txt",header=F )

Mo17_wt_06h.B1 <- read.table("Mo17_wt_06h.B1_counts.txt",header=F )
Mo17_wt_06h.B2 <- read.table("Mo17_wt_06h.B2_counts.txt",header=F )
Mo17_mt_06h.B1<- read.table("Mo17_mt_06h.B1_counts.txt",header=F )
Mo17_mt_06h.B2 <- read.table("Mo17_mt_06h.B2_counts.txt",header=F )
Mo17_wt_22h.B1 <- read.table("Mo17_wt_22h.B1_counts.txt",header=F )
Mo17_wt_22h.B2 <- read.table("Mo17_wt_22h.B2_counts.txt",header=F )
Mo17_mt_22h.B1<- read.table("Mo17_mt_22h.B1_counts.txt",header=F )
Mo17_mt_22h.B2 <- read.table("Mo17_mt_22h.B2_counts.txt",header=F )
Mo17_wt_48h.B1 <- read.table("Mo17_wt_48h.B1_counts.txt",header=F )
Mo17_wt_48h.B2 <- read.table("Mo17_wt_48h.B2_counts.txt",header=F )
Mo17_mt_48h.B1<- read.table("Mo17_mt_48h.B1_counts.txt",header=F )
Mo17_mt_48h.B2 <- read.table("Mo17_mt_48h.B2_counts.txt",header=F )


# Merge single htseq-count files into one:
files <- c(
        "B73_wt_06h.B1_counts.txt","B73_wt_06h.B2_counts.txt",
        "B73_mt_06h.B1_counts.txt","B73_mt_06h.B2_counts.txt",
        "B73_wt_22h.B1_counts.txt","B73_wt_22h.B2_counts.txt",
        "B73_mt_22h.B1_counts.txt","B73_mt_22h.B2_counts.txt",
        "B73_wt_48h.B1_counts.txt","B73_wt_48h.B2_counts.txt",
        "B73_mt_48h.B1_counts.txt","B73_mt_48h.B1_counts.txt",
        "Mo17_wt_06h.B1_counts.txt","Mo17_wt_06h.B2_counts.txt",
        "Mo17_mt_06h.B1_counts.txt","Mo17_mt_06h.B2_counts.txt",
        "Mo17_wt_22h.B1_counts.txt","Mo17_wt_22h.B2_counts.txt",
        "Mo17_mt_22h.B1_counts.txt","Mo17_mt_22h.B2_counts.txt",
        "Mo17_wt_48h.B1_counts.txt","Mo17_wt_48h.B2_counts.txt",
        "Mo17_mt_48h.B1_counts.txt","Mo17_mt_48h.B2_counts.txt"
        )
c <- 1
for (filename in files) {
        if(c == 1){ # if it is the first file just read file
                MaizeAll = read.table(filename,sep="\t")
        }
        else{ # else merge the other files
                tmp = read.table(filename,sep="\t")
                names(tmp) = c(c,c+1)
                MaizeAll = merge(MaizeAll,tmp,by=c(1))
        }
        c = c+1
}
names(MaizeAll) = c("GeneID",
                    "B73_wt_06h.B1","B73_wt_06h.B2","B73_mt_06h.B1","B73_mt_06h.B2",
                    "B73_wt_22h.B1","B73_wt_22h.B2","B73_mt_22h.B1","B73_mt_22h.B2",
                    "B73_wt_48h.B1","B73_wt_48h.B2","B73_mt_48h.B1","B73_48_22h.B2",
                    "Mo17_wt_06h.B1","Mo17_wt_06h.B2","Mo17_mt_06h.B1","Mo17_mt_06h.B2",
                    "Mo17_wt_22h.B1","Mo17_wt_22h.B2","Mo17_mt_22h.B1","Mo17_mt_22h.B2",
                    "Mo17_wt_48h.B1","Mo17_wt_48h.B2","Mo17_mt_48h.B1","Mo17_mt_48h.B2")
                            
                    
                    
str(MaizeAll)
names(MaizeAll)
###get protein coding genes###
#htseqAllCount_BL <- merge(htseqAllCount_BL,gene_anno,by=c(1))
#htseqAllCount_BL <- htseqAllCount_BL[htseqAllCount_BL$V2=="protein_coding_gene",]
#htseqAllCount_BL = htseqAllCount_BL[,-32]
####

# make table of only 6 numeric values for downstream analysis
cm <- MaizeAll[,-1]
rownames(cm) <- MaizeAll[,1]

# build DGEList
group <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
y <- DGEList(counts = cm, group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12))
str(y)
dim(y)


# paramenter for filtering low expressed genes
min.cpm <- 2
n.min.cpm <- 2
keep <- rowSums(cpm(y)>min.cpm) >= n.min.cpm
table(keep)
y <- y[keep,]
dim(y)

y$samples$lib.size <- colSums(y$counts)


# TMM normalization
y <- calcNormFactors(y, method="TMM")
y$samples


# prepare for edgeR glm
PRT  <- factor(c("B73_wt_06h", "B73_wt_06h", "B73_mt_06h","B73_mt_06h","B73_wt_22h","B73_wt_22h","B73_mt_22h","B73_mt_22h",
                 "B73_wt_48h","B73_wt_48h","B73_mt_48h","B73_mt_48h","Mo17_wt_06h", "Mo17_wt_06h", "Mo17_mt_06h","Mo17_mt_06h",
                 "Mo17_wt_22h","Mo17_wt_22h","Mo17_mt_22h","Mo17_mt_22h","Mo17_wt_48h","Mo17_wt_48h","Mo17_mt_48h","Mo17_mt_48h"), 
               
               levels= c("B73_wt_06h","B73_mt_06h","B73_wt_22h","B73_mt_22h","B73_wt_48h","B73_mt_48h",
                         "Mo17_wt_06h","Mo17_mt_06h","Mo17_wt_22h","Mo17_mt_22h","Mo17_wt_48h","Mo17_mt_48h"))
sample.names <- c("B73_wt_06h_1","B73_wt_06h_2","B73_mt_06h_1","B73_mt_06h_2",
                  "B73_wt_22h_1","B73_wt_22h_2","B73_mt_22h_1","B73_mt_22h_2",
                  "B73_wt_48h_1","B73_wt_48h_2","B73_mt_48h_1","B73_48_22h_2",
                  "Mo17_wt_06h_1","Mo17_wt_06h_2","Mo17_mt_06h_1","Mo17_mt_06h_2",
                  "Mo17_wt_22h_1","Mo17_wt_22h_2","Mo17_mt_22h_1","Mo17_mt_22h_2",
                  "Mo17_wt_48h_1","Mo17_wt_48h_2","Mo17_mt_48h_1","Mo17_48_22h_2")

targets <- as.data.frame(cbind(sample.names,PRT))
design <- model.matrix(~0+PRT)


my.contrasts <- makeContrasts(
        c1 = (PRTB73_wt_06h-PRTB73_mt_06h), # comparison: B73_wt_06h / B73_mt_06h
        c2 = (PRTB73_wt_22h-PRTB73_mt_22h),  # comparison: B73_wt_22h/B73_mt_22h
        c3 = (PRTB73_wt_48h-PRTB73_mt_48h), # comparison B73_wt_48h/ B73_mt_48h
        c4 = (PRTMo17_wt_06h-PRTMo17_mt_06h), # comparison Mo17_wt_06h / Mo17_mt_06h
        c5 = (PRTMo17_wt_22h-PRTMo17_mt_22h), # comparison Mo17_wt_22h/ Mo17_mt_22h
        c6 = (PRTMo17_wt_48h-PRTMo17_mt_48h), # comparison Mo17_wt_48h/ Mo17_mt_48h
        c7 =(PRTB73_wt_06h-PRTMo17_wt_06h), #  comparison: B73_wt_06h / Mo17_wt_06h 
        c8 =(PRTB73_wt_06h-PRTMo17_mt_06h), #  comparison: B73_wt_06h / Mo17_mt_06h 
        c9 =(PRTB73_wt_22h-PRTMo17_wt_22h), #  comparison: B73_wt_22h / Mo17_wt_22h
        c10 =(PRTB73_wt_22h-PRTMo17_mt_22h), #  comparison: B73_wt_22h / Mo17_mt_22h
        c11 =(PRTB73_wt_48h-PRTMo17_wt_48h),# comparison B73_wt_48h / Mo17_wt_48h
        c12=(PRTB73_wt_48h-PRTMo17_mt_48h),# #comparison B73_wt_48h/Mo17_mt_48h
        c13=(PRTB73_mt_06h-PRTMo17_mt_06h),# comparison B73.MT.06h_vs_Mo17.MT.06h  
        c14=(PRTB73_wt_22h-PRTMo17_mt_22h),# comparison B73.MT.22h_vs_Mo17.MT.22h
        c15=(PRTB73_mt_48h-PRTMo17_mt_48h),levels=design #B73.MT.48h_vs_Mo17.MT.48h
)
interesting.contrasts <- c("c1", "c2", "c3", "c4", "c5", "c6","c7", "c8",
                           "c9", "c10", "c11", "c12", "c13","c14","c15")


# variance estimate
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

# QC plots

#plotMDS(y) # MDS plot
#plotMeanVar(y) # Mean vs Variance
#plotBCV(y) #


# EdgeR analysis & report
fit <- glmFit(y,design)
fit

# fdr threshold
fdr.t <- 1
project.name <- "RNAseq_GF_Maize"
cat("experiment: ", project.name, "\n")
cat("thresholds: ", min.cpm, " ", n.min.cpm,"\n")
for (my.contrast in interesting.contrasts) {
        lrt <- glmLRT(fit, contrast=my.contrasts[,my.contrast])
        etable <- topTags(lrt, n=nrow(lrt$table), adjust.method="BH")
        etable <- etable$table[etable$table$FDR<fdr.t,]
        etable <- etable[ order(etable$FDR), ]  
        cat(my.contrast," ", dim(etable)[1], " genes\n")
        # write result (c1,...,c6) in workign dir. 
        write.table( etable[,], file=paste(project.name,my.contrast, min.cpm,n.min.cpm,sep="."), row.names=TRUE)
}

#### ENDS HERE ####


################################################
#RPKM CALCULTING (Qiwen TAIR10.gene.length) #### 
################################################

# calculate gene lenght 
len = read.table("TAIR10.gene.length",sep="\t")
names(len) = c("ID","length")
htseqAllCount_BL = merge(htseqAllCount_BL,len,by=c(1))
x = htseqAllCount_BL[,c(2:31)] # should this be 8 instead of 7 Qiwen?!!!!!!!!!!!!!############!!!!!!
rownames(x) = htseqAllCount_BL[,1]

rpkm = rpkm(x,htseqAllCount_BL[,32],normalized.lib.sizes=TRUE)
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

str(rpkm)

### Merge Table ### 

# merge p-values/FDR table with RPKM table
etable$name = rownames(etable)
edgeR_table = merge(rpkm,etable,by=c("name"))
str(table)


#final table
edgeR_table 

#write.table(edgeR_table, file="YFPs_edgeR_DEGs_table")


#### pheatmap of ALL conditions and reps 

edgeR_table_hm =edgeR_table[,1:31]
test=edgeR_table_hm[2:31]
test=as.matrix(test)
test
pheatmap(test, scale="row", color = colorRampPalette(c("navy", "black", "firebrick3"))(50))


scale="row", color = colorRampPalette(c("navy","lightgreen", "black", "firebrick3"))(50))


test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)
pheatmap(test, display_numbers = TRUE)
pheatmap(test, display_numbers = TRUE, number_format = "%.1e")
pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
                                                                            "1e-4", "1e-3", "1e-2", "1e-1", "1"))
pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")

