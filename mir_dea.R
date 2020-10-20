####miRNA####

####binding miRNA reads data####
demo <- read.csv("data/miRNAs_expressed_all_samples_07_09_2020_t_16_02_45.csv", sep = "\t")
mir <- data.frame(demo[, c(1, 3)])
file_list <- list.files(path="E:/work/Project/mideep/data")
for (i in 1:length(file_list)){
    temp <- read.csv(paste0("data/", file_list[i]), sep = "\t")
    mir <- cbind(mir, temp[, 2])
}
colnames(mir) <- c("miRNA", "precursor", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19",
                   "m1", "m2", "m4", "m5", "m6", "m7", "m8", "m9")
mir_count <- mir[!duplicated(mir$miRNA), -2]
mir_count <- mir_count[rowSums(mir_count[, -1]) > 15, ]
#write.csv(mir_count, "mir_count.csv", row.names = FALSE)


####different expression analysis####
library(DESeq2)
mir_count <- read.csv("mir_count.csv")
mir_count <- data.frame(mir_count[,-1], row.names = mir_count[,1])
group <- read.csv("group.csv")


dem <- function(data, design){
    de <- DESeqDataSetFromMatrix(countData = data ,colData = group,design = design)
    #de <- de[rowSums(counts(de) >= 10) >= 3,]# prefilter
    de <- DESeq(de)
    return(de)
}
get_de <- function(t, c){
    x <- results(t, contrast = c)
    x <- subset(as.data.frame(x), padj < 0.1 & abs(log2FoldChange) > 2)
    return(x)
}

#obesity or not
de_f <- dem(mir_count, ~fat)
dem_f <- results(de_f)
dem_f <- subset(as.data.frame(dem_f), padj < 0.1 & abs(log2FoldChange) > 2)

#lung cancer or not
de_c <- dem(mir_count, ~cancer)
dem_c <- results(de_c)
dem_c <- subset(as.data.frame(dem_c), padj < 0.1 & abs(log2FoldChange) > 2)

#control, lung cancer, obesity, lung cancer & obesity
de_n <- dem(mir_count, ~name)

dem_l_c <- get_de(de_n, c("name", "l", "c"))
dem_o_c <- get_de(de_n, c("name", "o", "c"))
dem_ol_c <- get_de(de_n, c("name", "ol", "c"))
dem_l_o <- get_de(de_n, c("name", "l", "o"))
dem_l_ol <- get_de(de_n, c("name", "l", "ol"))
dem_o_ol <- get_de(de_n, c("name", "o", "ol"))

#lung cancer, obesity interaction
de_fc <- dem(mir_count, ~fat*cancer)

res_c_i <- results(de_fc, name = "cancer")
res_f_i <- results(de_fc, name = "fat")
res_b_i <- results(de_fc, name = "fat.cancer")

dem_i_l <- subset(as.data.frame(res_c_i), padj < 0.1 & abs(log2FoldChange) > 2)
dem_i_o <- subset(as.data.frame(res_f_i), padj < 0.1 & abs(log2FoldChange) > 2)
dem_i_i <- subset(as.data.frame(res_b_i), padj < 0.1 & abs(log2FoldChange) > 2)


#venn
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(row.names(dem_c), row.names(dem_l_c), row.names(dem_l_o)),
    category.names = c("cancer" , "cancer/\n control " , "cancer/\n obesity"),
    filename = "venn_diagramm.png",
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 27, 0),
    cat.dist = c(0.3, 0.2, 0.1),
    cat.fontfamily = "sans",
    rotation = 1
)


#dem intersect

insec_c_il <- intersect(row.names(dem_i_l), row.names(dem_c))#7
insec_c_lc <- intersect(row.names(dem_l_c), row.names(dem_c))#7
insec_c_lo <- intersect(row.names(dem_l_o), row.names(dem_c))#19
insec_il_lc <- intersect(row.names(dem_i_l), row.names(dem_l_c))#23
insec_il_lo <- intersect(row.names(dem_i_l), row.names(dem_l_o))#3
insec_lc_lo <- intersect(row.names(dem_l_c), row.names(dem_l_o))#3

dem_all <- union(union(row.names(dem_l_c), row.names(dem_c)), row.names(dem_l_o))

#write.csv(dem_all, "dem.csv", row.names = FALSE)


####target predict####

#targetscan
library(targetscan.Hs.eg.db)
library(isomiRs)

dem_target <- mirna2targetscan(dem_all)
#write.csv(dem_target, "targetscan.csv", row.names = FALSE)


#####
#miRNAtap
library(miRNAtap)
library(org.Hs.eg.db)#human database
dem_all <- readLines("dem.csv")

#multiple variable
# for (n in 1:length(dem_all)){
#     name <- paste0("m", n)
#     target <- data.frame(getPredictedTargets(dem_all[n], species ="hsa", method ="geom", min_src = 2))
#     assign(name, target)#assignment
# }
# 1: In getPredictedTargets(dem_all[n], species = "hsa", method = "geom",  :no targets found for mirna miR-1255b-5
# 2: In getPredictedTargets(dem_all[n], species = "hsa", method = "geom",  :no targets found for mirna miR-10401-3

# mrna_list <- list()
# for (n in 1:length(dem_all)){
#     target <- data.frame(getPredictedTargets(dem_all[n], species ="hsa", method ="geom", min_src = 2))
#     mrna_list[[n]] <- row.names(target)
# }
# 
# m_len <- sapply(mrna_list, length)
# seq_max <- seq_len(max(m_len))
# mrna_down <- data.frame(sapply(mrna_list, "[", i = seq_max))
# colnames(mrna_down) <- dem_all

#write.csv(mrna_down, "mrna_down.csv", row.names = FALSE)


#mRNA count
mrna_down <- read.csv("mrna_down.csv")
gene_name <- c("TAFA5", "SLC2A12", "MFGE8", "VGLL3", "VAMP8", "ATXN7L3", "TFDP3", "RBFOX3")
gene_id <- c(25817, 154091, 4240, 389136, 8673, 56970, 51270, 146713)
#gene_id <- mapIds(org.Hs.eg.db, gene_name, "ENTREZID", "SYMBOL")

#count related miRNA
down_all <- unique(unlist(mrna_down))
down_count <- c()
for (n in 1:length(down_all)){
    x <- length(which(mrna_down == down_all[n]))
    down_count <- c(down_count, x)
}
down_count_all <- data.frame(cbind(down_all, down_count))

down_count_all <- down_count_all[order(down_count_all$down_count, decreasing = TRUE), ]
library(org.Hs.eg.db)
symbol <- mapIds(org.Hs.eg.db, keys = as.character(down_count_all$down_all), column = "SYMBOL", keytype = "ENTREZID")

#write.csv(data.frame(unlist(down_count_all)), "down_count_all.csv", row.names = FALSE)

#related miRNA list
for (n in 1:length(gene_id)){
    k <- which(mrna_down == gene_id[n], arr.ind = T)
    if (length(k) != 0){
        mirna <- sapply(k[, 2], function(x) colnames(mrna_down)[x])
        a <- cbind(mirna, k[, 1])
        assign(paste0("mi_", gene_id[n]), a)
    }else{
        next
    }
}
#####

#kegg####
library(clusterProfiler)
library(org.Hs.eg.db)
mRNA_enrich <- enrichKEGG(gene = mrna_down[, 1], organism = "hsa", keyType = "kegg", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)

a <- data.frame(mRNA_enrich)




####mRNA####
library(limma)
mrna_raw <- read.csv("mRNA_expression.csv", stringsAsFactors = FALSE)
group_m <- read.csv("group_m.csv", stringsAsFactors = FALSE)
mrna_count <- data.frame(mrna_raw[,20:37], row.names = mrna_raw$GeneSymbol)


#control, lung cancer, obesity, lung cancer & obesity
design_g <- model.matrix(~0+group_m$name)
colnames(design_g) <- c("c","l","o","ol")
#linear model
fit_g <- lmFit(mrna_count, design_g)
contrast.matrix_g <- makeContrasts(o-c, l-c, ol-c, o-l, o-ol, l-ol, levels = design_g)
fit_g2 <- contrasts.fit(fit_g, contrast.matrix_g)
fit_g2 <- eBayes(fit_g2)

deg_o_c <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 1, sort.by = "logFC", lfc = 2, p.value = 0.05)
deg_l_c <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 2, sort.by = "logFC", lfc = 2, p.value = 0.05)
deg_ol_c <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 3, sort.by = "logFC", lfc = 2, p.value = 0.05)
deg_o_l <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 4, sort.by = "logFC", lfc = 2, p.value = 0.05)
deg_o_ol <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 5, sort.by = "logFC", lfc = 2, p.value = 0.05)
deg_l_ol <- topTable(fit_g2, adjust = "fdr", number=Inf, coef = 6, sort.by = "logFC", lfc = 2, p.value = 0.05)


#interaction design different expression analysis
design_i <- model.matrix(~group_m$fat*group_m$cancer)
colnames(design_i) <- c("intercept","fat","cancer","interaction")
#linear model
fit_i <- lmFit(mrna_count, design_i)
fit_i <- eBayes(fit_i)

deg_fat <- topTable(fit_i, adjust = "fdr", number=Inf, coef = "fat", lfc = 2, p.value = 0.05)
deg_cancer <- topTable(fit_i, adjust = "fdr", number=Inf, coef = "cancer", lfc = 2, p.value = 0.05)
deg_interaction <- topTable(fit_i, adjust = "fdr", number=Inf, coef = "interaction", lfc = 2, p.value = 0.05)


#venn
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(row.names(deg_fat), row.names(deg_cancer), row.names(deg_interaction)),
    category.names = c("fat" , "cancer" , "interaction"),
    filename = "deg_i.png",
)

venn.diagram(
    x = list(row.names(deg_l_c), row.names(deg_o_c), row.names(deg_ol_c)),
    category.names = c("l_c", "o_c", "ol_c"),
    filename = "deg_g.png",
)


#intersection
insec_i <- Reduce(intersect, list(row.names(deg_fat), row.names(deg_cancer), row.names(deg_interaction)))
insec_g <- Reduce(intersect, list(row.names(deg_l_c), row.names(deg_o_c), row.names(deg_ol_c)))

insec_m <- union(insec_i, insec_g)

insec_m_id <- mrna_raw[mrna_raw$GeneSymbol %in% insec_m, ]$EntrezID
gene_id <- c(25817, 154091, 4240, 389136, 8673, 56970, 51270, 146713)
gene_id %in% insec_m_id


insec_m_symbol <- mrna_raw[mrna_raw$GeneSymbol %in% insec_m, c(58, 66)]
colnames(insec_m_symbol) <- c("symbol", "gene")


####miRNA-mRNA interaction####
dem_target <- read.csv("targetscan.csv", stringsAsFactors = FALSE)
dem_deg <- dem_target[dem_target$gene %in% insec_m_id, ]
dem_deg[dem_deg$PCT == "NULL", ]$PCT <- NA
dem_deg$PCT <- as.numeric(dem_deg$PCT)
deg_mi <- unique(dem_deg$gene)

gene_id %in% deg_mi


library(clusterProfiler)
library(org.Hs.eg.db)
deg_mi_enrich <- enrichKEGG(gene = insec_m_id, organism = "hsa", keyType = "kegg", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
c <- data.frame(deg_mi_enrich)
#write.csv(c, "path.csv", row.names = FALSE)

deg_mi_symbol <- mapIds(org.Hs.eg.db, keys = as.character(deg_mi), column = "SYMBOL", keytype = "ENTREZID")


####miRNA-mRNA regulation network####

#link
deg_link <- dem_deg[, c(2, 3, 4, 6)]
deg_link <- merge(deg_link, insec_m_symbol, by = "gene")

#node
deg_node <- data.frame(name = c(unique(deg_link$mir), unique(deg_link$symbol)),
                       group = c(rep("miRNA", 28), rep("gene", 44)),
                       size = c(rep(2, 28), rep(1, 44)))

#write.csv(deg_node, "deg_node.csv", row.names = FALSE)

#deg_link need to convert node name to number, start from 0
deg_link$mir <- match(deg_link$mir, deg_node$name)-1
deg_link$symbol <- match(deg_link$symbol, deg_node$name)-1


library(networkD3)
library(RColorBrewer)
brewer.pal(11,"Set1")
my_color <- 'd3.scaleOrdinal() .domain(["miRNA", "gene"]) 
.range(["#FF7F00", "#55BCCF"])'
forceNetwork(Links = deg_link, Nodes = deg_node,
             Source = "mir", Target = "symbol",
             NodeID = "name", Group = "group",
             fontSize = 20, fontFamily = "arial", 
             Nodesize = "size", colourScale = my_color, 
             clickAction = NULL, opacityNoHover = 1, 
             charge=-500, zoom = TRUE)
