library("circlize")
library("ComplexHeatmap")
library("ggplot2")
library('DESeq2')
library('rjson')

#setwd('/Users/Hsieh/Documents/GitHub/QuantTE/')

json_data <- fromJSON(file='test.json')
outdir <- json_data$output_dir
label_list <- json_data$read_label

merge_xprs <- data.frame()

for(label in label_list){
    xprs_file <- read.table(paste0(outdir, '/kallisto/', label, '/merge_abundance.tsv'), sep='\t', header=T, stringsAsFactors=F)
    name <- xprs_file$target_id
    rownames(xprs_file) <- name
    colnames(xprs_file) <- c('target_id', 'avg_length', paste0(label, '_est_counts'), paste0(label, '_est_tpm'))
    if(nrow(merge_xprs) == 0){
        merge_xprs <- xprs_file
        next
    } else {
        merge_xprs <- cbind(merge_xprs, xprs_file[merge_xprs$target_id, c(3, 4)])
    }
}

count_data <- merge_xprs[, seq(3, ncol(merge_xprs), 2)]
count_data <- round(count_data)
colnames(count_data) <- label_list
rownames(count_data) <- merge_xprs$target_id
col_data <- data.frame(condition=label_list)
rownames(col_data) <- label_list

DESeq_data <- DESeqDataSetFromMatrix(countData=count_data, colData=col_data, design=~condition)
DESeq_model <- DESeq(DESeq_data)
DESeq_result <- results(DESeq_model, cooksCutoff=F)

merge_xprs$DESeq2_P = DESeq_result$pvalue
merge_xprs$DESeq2_P_ADJ = DESeq_result$padj

merge_xprs[grepl("NR_\\d+", merge_xprs$target_id, perl=TRUE), "category"] = "ncRNA"
merge_xprs[grepl("NM_\\d+", merge_xprs$target_id, perl=TRUE), "category"] = "mRNA"
merge_xprs[grepl("LINE", merge_xprs$target_id, perl=TRUE), "category"] = "LINE"
merge_xprs[grepl("LTR", merge_xprs$target_id, perl=TRUE), "category"] = "LTR"
merge_xprs[grepl("DNA", merge_xprs$target_id, perl=TRUE), "category"] = "DNA"
merge_xprs[grepl("Satellite", merge_xprs$target_id, perl=TRUE), "category"] = "Satellite"
merge_xprs[grepl("SINE", merge_xprs$target_id, perl=TRUE), "category"] = "SINE"
merge_xprs$category <- factor(merge_xprs$category)

write.table(merge_xprs, 'Kallisto_DESeq2_Report.tsv', sep='\t', row.names=F, quote=F)

deg_condition <- (merge_xprs[, 3] + 1)/(merge_xprs[, 5] + 1) > 2 | (merge_xprs[, 5] + 1)/(merge_xprs[, 3] + 1) > 2
deg_condition <- deg_condition & merge_xprs[, 4] > 0.1 & merge_xprs[, 6] > 0.1
deg <- merge_xprs[deg_condition, ]

scaled_xprs <- log2(as.matrix(deg[, c('A_est_tpm', 'B_est_tpm')]) + 1)
colnames(scaled_xprs) <- c('A', 'B')

gene_condition <- deg$category == 'mRNA' | deg$category == 'ncRNA'
gene_deg <- deg[gene_condition, ]
gene_scaled_xprs <- scaled_xprs[gene_condition, ]
te_deg <- deg[!gene_condition, ]
te_scaled_xprs <- scaled_xprs[!gene_condition, ]

category_cols <- c('firebrick1', 'dodgerblue', 'darkorchid1', 'seagreen1', 'sienna1', 'gold', 'gray31')
names(category_cols) <- levels(merge_xprs$category)

jpeg('DEG_All_Heatmap.jpeg', height=2000, width=2000, res=250)
print(Heatmap(scaled_xprs, name="log2(TPM + 1)", km=3 ,col=colorRamp2(c(0, max(scaled_xprs)), c("#6D7D91", "#f0650e")),
              show_row_names=F, show_column_names=T) +
      Heatmap(deg$DESeq2_P, name="DESeq2 P", col=colorRamp2(c(0, 1), c("#6D7D91", "#f0650e")), show_row_names=F, width=unit(5, "mm")) +
      Heatmap(deg$category, name="Category", col=category_cols, width=unit(5, "mm")))
dev.off()

jpeg('DEG_mRNA_ncRNA_Heatmap.jpeg', height=2000, width=2000, res=250)
print(Heatmap(gene_scaled_xprs, name="log2(TPM + 1)", km=3 ,col=colorRamp2(c(0, max(scaled_xprs)), c("#6D7D91", "#f0650e")),
              show_row_names=F, show_column_names=T)  +
      Heatmap(gene_deg$DESeq2_P, name="DESeq2 P", col=colorRamp2(c(0, 1), c("#6D7D91", "#f0650e")), show_row_names=F, width=unit(5, "mm")) +
      Heatmap(gene_deg$category, name="Category", col=category_cols, width=unit(5, "mm")))
dev.off()

jpeg('DEG_TE_Heatmap.jpeg', height=2000, width=2000, res=250)
print(Heatmap(te_scaled_xprs, name="log2(TPM + 1)", km=3 ,col=colorRamp2(c(0, max(scaled_xprs)), c("#6D7D91", "#f0650e")),
              show_row_names=F, show_column_names = T)  +
      Heatmap(te_deg$DESeq2_P, name="DESeq2 P", col=colorRamp2(c(0, 1), c("#6D7D91", "#f0650e")), show_row_names=F, width=unit(5, "mm")) +
      Heatmap(te_deg$category, name="Category", col=category_cols, width=unit(5, "mm")))
dev.off()

for(target in c('DNA', 'LINE', 'LTR', 'Satellite', 'SINE')){
    te <- subset(deg, deg$category == target)
    if(nrow(te) == 0){
        next
    }
    cat_data <- data.frame(Name=te$target_id, TPM=te$A_est_tpm, condition='A')
    cat_data <- rbind(cat_data, data.frame(Name=te$target_id, TPM=te$B_est_tpm, condition='B'))
    jpeg(paste0(target, '_Barplot.jpeg'), height=2000, width=2000, res=250)
    f <- ggplot(data=cat_data, aes(x=Name, y=TPM, fill=condition)) +
         geom_bar(stat="identity", position=position_dodge()) + 
         scale_fill_brewer(palette="Paired") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(f)
    dev.off()
}
