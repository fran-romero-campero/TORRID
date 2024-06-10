
# Wheat TranscriptOmic Response to a micRoalgae bIostimulant under
# Drought stress (TORRID)
# Authors Marcos Ramos-González and Francisco J. Romero-Campero

# Load data, create count matrix and visualize samples PCA
.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")
library(DESeq2)

gene.expression.ta <- read.csv(file = "../TORRID/gene_count_matrix.csv")
head(gene.expression.ta)
rownames(gene.expression.ta) <- gene.expression.ta$gene_id
gene.expression.ta <- gene.expression.ta[,2:ncol(gene.expression.ta)]
head(gene.expression.ta)
nrow(gene.expression.ta)
gene.expression.ta <- gene.expression.ta[-which(apply(X = gene.expression.ta,MARGIN = 1,FUN = sum) == 0),]
nrow(gene.expression.ta) # 75198 of 99386 genes showed expression in at least one sample
rownames(gene.expression.ta) <- unlist(strsplit(x = rownames(gene.expression.ta),split=".v2.2"))

experimental.design <- read.csv(file="../TORRID/experimental_design.csv")

colnames(gene.expression.ta) <- paste(experimental.design$treatment, 1:3,sep="_")
head(gene.expression.ta)


dds.wheat <- DESeqDataSetFromMatrix(countData=gene.expression.ta, colData=experimental.design, design = ~ treatment)
dds.wheat <- DESeq(dds.wheat)
rld <- rlog(dds.wheat)
colData(rld)

png("pca.png", width = 3000, height = 2000, res = 300)
plotPCA(object = rld,intgroup=c("treatment"))
dev.off()


# Change PCA format
library(FactoMineR)
library(factoextra)

counted_wheat <- log2(counts(dds.wheat, normalized=TRUE) + 1)

pca.gene.expression <- data.frame(colnames(counted_wheat),
                                  t(counted_wheat))
colnames(pca.gene.expression)[1] <- "Sample"
pca.gene.expression[1:6, 1:6]

res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = T,quali.sup = 1,  )
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 4)
fviz_dend(res.hcpc,k=4,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 1400      # Augment the room for labels
)


fviz_pca_ind(res.pca, col.ind = c(rep("H20_180", 3), rep("LRM_180",3), rep("LRM_30",3), rep("H2O_30",3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


##### Differential gene expression analysis under drought stress (not treated)
res.h2o.30.180 <- results(dds.wheat,contrast=c("treatment","H2O_30","H2O_180"))
res.h2o.30.180$padj[is.na(res.h2o.30.180$padj)] <- 1

genes.h2o.30.180 <- rownames(res.h2o.30.180)

activated.genes.h2o.30.180 <- genes.h2o.30.180[res.h2o.30.180$log2FoldChange > 1 & res.h2o.30.180$padj < 0.05]
length(activated.genes.h2o.30.180)

repressed.genes.h2o.30.180 <- genes.h2o.30.180[res.h2o.30.180$log2FoldChange < -1 & res.h2o.30.180$padj < 0.05]
length(repressed.genes.h2o.30.180)

logfc.h2o.30.180 <- res.h2o.30.180$log2FoldChange
names(logfc.h2o.30.180) <- genes.h2o.30.180
logqval.h2o.30.180 <- -log10(res.h2o.30.180$padj)
names(logqval.h2o.30.180) <- genes.h2o.30.180

# Volcano plot
plot(logfc.h2o.30.180[logqval.h2o.30.180 > 0],logqval.h2o.30.180[logqval.h2o.30.180 > 0],pch=19,cex=0.7,xlim=c(-12,12),ylim=c(0,35),col="grey",xlab="log2FC",ylab="-log10(q-value)")
points(logfc.h2o.30.180[activated.genes.h2o.30.180],logqval.h2o.30.180[activated.genes.h2o.30.180],cex=0.7,col="red",pch=19)
points(logfc.h2o.30.180[repressed.genes.h2o.30.180],logqval.h2o.30.180[repressed.genes.h2o.30.180],cex=0.7,col="blue",pch=19)

# ggplot2 volcano plot
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

df <- as.data.frame(res.h2o.30.180)
head(df)

df$diffexpressed <- "NO"
df$diffexpressed[df$log2FoldChange > 1 & df$padj < 0.05] <- "UP"
df$diffexpressed[df$log2FoldChange < -1 & df$padj < 0.05] <- "DOWN"
df$names <- rownames(df)

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

volcano_plot_h2o <-  ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />

png("volcano_h2o.png", width = 3000, height = 2000, res = 300)
plot(volcano_plot_h2o)
dev.off()

# GO enrichment analysis under drought stress (not treated)
taestivium2atha <- read.table(file="../TORRID/taestivium_atha.csv",header=T,fill = T)

head(taestivium2atha)

activated.genes.h2o.30.180.atha <- subset(taestivium2atha,locusName %in% activated.genes.h2o.30.180)[[2]]

activated.genes.h2o.30.180.atha <- unique( activated.genes.h2o.30.180.atha[activated.genes.h2o.30.180.atha != ""])

length(activated.genes.h2o.30.180)
length(activated.genes.h2o.30.180.atha)


repressed.genes.h2o.30.180.atha <- subset(taestivium2atha,locusName %in% repressed.genes.h2o.30.180)[[2]]

repressed.genes.h2o.30.180.atha <- unique( repressed.genes.h2o.30.180.atha[repressed.genes.h2o.30.180.atha != ""])

length(repressed.genes.h2o.30.180)
length(repressed.genes.h2o.30.180.atha)


library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)
#install.packages("../TORRID/org.Taestivum.eg.db/",repos = NULL, type="source")
library(org.Taestivum.eg.db)

# Activated genes H20 under drought no treatment
activated.genes.h2o.30.180.atha.enrich.go <- enrichGO(gene          = activated.genes.h2o.30.180.atha,
                                                      OrgDb         = org.At.tair.db,
                                                      ont           = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.05,
                                                      readable      = FALSE,
                                                      keyType = "TAIR")

activated.genes.h2o.30.180.atha.enrich.go.df <- as.data.frame(activated.genes.h2o.30.180.atha.enrich.go)

activated.genes.h2o.30.180.enrich.go <- enrichGO(gene = activated.genes.h2o.30.180, 
                                                 OrgDb = org.Taestivum.eg.db,
                                                 universe = rownames(gene.expression.ta), 
                                                 ont           = "BP",
                                                 pAdjustMethod = "BH",
                                                 pvalueCutoff  = 0.05,
                                                 readable      = FALSE,
                                                 keyType = "GID")

activated.genes.h2o.30.180.enrich.go.df <- as.data.frame(activated.genes.h2o.30.180.enrich.go)
head(activated.genes.h2o.30.180.enrich.go.df)

write.table(x = activated.genes.h2o.30.180.enrich.go.df,file="activated_genes_h2o_30_180_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


intersect(activated.genes.h2o.30.180.atha.enrich.go.df$Description,
          activated.genes.h2o.30.180.enrich.go.df$Description)

setdiff(activated.genes.h2o.30.180.atha.enrich.go.df$Description,
        activated.genes.h2o.30.180.enrich.go.df$Description)

setdiff(activated.genes.h2o.30.180.enrich.go.df$Description,
        activated.genes.h2o.30.180.atha.enrich.go.df$Description)


treeplot(x = pairwise_termsim(activated.genes.h2o.30.180.enrich.go))

{
library(treemap)

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006766","vitamin metabolic process",1.256768546208505,2.4660664317727994,0.8207383156461785,0.03913661,"vitamin metabolic process"),
                     c("GO:0006563","L-serine metabolic process",0.2663767935051533,1.4768020745959938,0.8162945558299953,0.29970322,"vitamin metabolic process"),
                     c("GO:0006817","phosphate ion transport",0.2426906405041399,10.279245651666049,0.9402625518575654,-0,"phosphate ion transport"),
                     c("GO:0015698","inorganic anion transport",0.8722912628806088,5.190747108887932,0.9486819408786794,0.25520912,"phosphate ion transport"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,1.3825896521155858,0.93668431916998,0.37687047,"phosphate ion transport"),
                     c("GO:0034486","vacuolar transmembrane transport",0.025445977479964865,2.6776125250934606,0.9484483038087386,0.19211745,"phosphate ion transport"),
                     c("GO:0071836","nectar secretion",0.00012077226816333578,3.0448540126109958,0.97007760379404,0.13980114,"phosphate ion transport"),
                     c("GO:0072530","purine-containing compound transmembrane transport",0.132753370113172,2.5908533053098775,0.9257403406439652,0.24048243,"phosphate ion transport"),
                     c("GO:0098656","monoatomic anion transmembrane transport",0.4548751919662879,2.516341325301193,0.9371903523548951,0.29818431,"phosphate ion transport"),
                     c("GO:1901264","carbohydrate derivative transport",0.41917836013343907,1.3162680441762238,0.9361451451937864,0.4160251,"phosphate ion transport"),
                     c("GO:1901679","nucleotide transmembrane transport",0.10048745659222856,2.3731648134601735,0.8955766626667911,0.39035476,"phosphate ion transport"),
                     c("GO:1905011","transmembrane phosphate ion transport from cytosol to vacuole",0.0009292070428077058,3.4069290072598006,0.9451671847180927,0.53391375,"phosphate ion transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.4047522359383369,1.313443483625142,0.9377320640914398,0.32867354,"phosphate ion transport"),
                     c("GO:0009813","flavonoid biosynthetic process",0.06199314467029023,4.098813776766317,0.8902418214470647,0.04437928,"flavonoid biosynthetic process"),
                     c("GO:0009312","oligosaccharide biosynthetic process",0.15951551884210055,1.890728518352544,0.8817172375612645,0.13137012,"flavonoid biosynthetic process"),
                     c("GO:0046349","amino sugar biosynthetic process",0.16032395361674495,1.4768020745959938,0.8237217682202084,0.14106407,"flavonoid biosynthetic process"),
                     c("GO:0016036","cellular response to phosphate starvation",0.07367847780013462,58.437263285695195,0.9337176454799332,0,"cellular response to phosphate starvation"),
                     c("GO:0000160","phosphorelay signal transduction system",2.0739975137672992,2.325112263035866,0.8787074548780844,0.38056179,"cellular response to phosphate starvation"),
                     c("GO:0002213","defense response to insect",0.000404217387322185,2.0641658594883507,0.9595305939122576,0.40657397,"cellular response to phosphate starvation"),
                     c("GO:0009735","response to cytokinin",0.0338384177472334,3.1808353925886697,0.9456338276874637,0.6140009,"cellular response to phosphate starvation"),
                     c("GO:0009736","cytokinin-activated signaling pathway",0.03294125232659148,2.9403628806797477,0.8778864568872735,0.27705991,"cellular response to phosphate starvation"),
                     c("GO:0009861","jasmonic acid and ethylene-dependent systemic resistance",0.00010844856733034232,2.039323022807327,0.9611294916787141,0.41577308,"cellular response to phosphate starvation"),
                     c("GO:0010046","response to mycotoxin",0.0002735861584924545,1.9919137875125033,0.9631214731700426,0.14319862,"cellular response to phosphate starvation"),
                     c("GO:0010438","cellular response to sulfur starvation",0.006469942937321558,2.0665452561891056,0.9410425515767487,0.69574796,"cellular response to phosphate starvation"),
                     c("GO:0061760","antifungal innate immune response",0.0015700394861233652,5.2440325251435125,0.9567569120910712,0.39453375,"cellular response to phosphate starvation"),
                     c("GO:0071369","cellular response to ethylene stimulus",0.038408046016107374,1.3548443793256149,0.934474712996174,0.65359876,"cellular response to phosphate starvation"),
                     c("GO:0090548","response to nitrate starvation",0.00011337804766353971,1.7875358893054079,0.9589557997083563,0.53468156,"cellular response to phosphate starvation"),
                     c("GO:0098869","cellular oxidant detoxification",0.7848373522893541,1.5594953414353876,0.9306523803872746,0.54607613,"cellular response to phosphate starvation"),
                     c("GO:0032107","regulation of response to nutrient levels",0.07313130548314971,5.316677594415389,0.834885327322508,-0,"regulation of response to nutrient levels"),
                     c("GO:0002684","positive regulation of immune system process",0.4153949839777101,1.3291303527470133,0.8703173946899884,0.56744823,"regulation of response to nutrient levels"),
                     c("GO:0080036","regulation of cytokinin-activated signaling pathway",0.003295357602742447,5.27125913127071,0.8952620871243417,0.41115455,"regulation of response to nutrient levels"),
                     c("GO:1900150","regulation of defense response to fungus",0.006139667754997334,2.005496928896609,0.8542335199281061,0.63401875,"regulation of response to nutrient levels"),
                     c("GO:1900367","positive regulation of defense response to insect",0.00017006707149530955,3.932513543308072,0.8573698014535559,0.52312945,"regulation of response to nutrient levels"),
                     c("GO:1905034","regulation of antifungal innate immune response",0.00037217576515640203,5.2440325251435125,0.8583797470084782,0.54390328,"regulation of response to nutrient levels"),
                     c("GO:2000068","regulation of defense response to insect",0.00025140349699306625,2.7748044400795004,0.873690815285751,0.53329505,"regulation of response to nutrient levels"),
                     c("GO:0046475","glycerophospholipid catabolic process",0.05117540007908858,19.4248594488226,0.7557130715577124,0.00712519,"glycerophospholipid catabolic process"),
                     c("GO:0006011","UDP-glucose metabolic process",0.04799835000434287,1.5875816119081179,0.7399927139996335,0.35670137,"glycerophospholipid catabolic process"),
                     c("GO:0006308","DNA catabolic process",0.09500094498137987,2.7621777376072436,0.8423503840725731,0.36845124,"glycerophospholipid catabolic process"),
                     c("GO:0006643","membrane lipid metabolic process",1.10151702785462,5.294723910021435,0.8415719209119874,0.59578979,"glycerophospholipid catabolic process"),
                     c("GO:0009110","vitamin biosynthetic process",1.0904306265852592,2.737879857421467,0.7413751174240135,0.59380742,"glycerophospholipid catabolic process"),
                     c("GO:0009129","pyrimidine nucleoside monophosphate metabolic process",0.29525615403719013,1.6050767418332463,0.7047695221495438,0.49165031,"glycerophospholipid catabolic process"),
                     c("GO:0009208","pyrimidine ribonucleoside triphosphate metabolic process",0.16502421311444862,1.8888538753238655,0.7035503701411095,0.45622658,"glycerophospholipid catabolic process"),
                     c("GO:0009247","glycolipid biosynthetic process",0.7048811812848927,8.962191092901932,0.7379178379605761,0.56300419,"glycerophospholipid catabolic process"),
                     c("GO:0009395","phospholipid catabolic process",0.12197752610480252,11.145253678990128,0.7410733796097034,0.69336156,"glycerophospholipid catabolic process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,5.252544183992574,0.8410648590005543,0.57787118,"glycerophospholipid catabolic process"),
                     c("GO:0019374","galactolipid metabolic process",0.01169519209051078,2.7231276464427947,0.829535208466126,0.68689132,"glycerophospholipid catabolic process"),
                     c("GO:0019375","galactolipid biosynthetic process",0.008735039150425755,2.926776963467647,0.8052385709160328,0.67287367,"glycerophospholipid catabolic process"),
                     c("GO:0032262","pyrimidine nucleotide salvage",0.038208402062612876,3.340190881781737,0.6533258183523374,0.32609505,"glycerophospholipid catabolic process"),
                     c("GO:0032958","inositol phosphate biosynthetic process",0.03756017539879743,3.2897278900751576,0.7424661718220883,0.37421144,"glycerophospholipid catabolic process"),
                     c("GO:0042360","vitamin E metabolic process",0.004374913795712673,3.031012464547948,0.8060126988455398,0.42128918,"glycerophospholipid catabolic process"),
                     c("GO:0046434","organophosphate catabolic process",0.9710632603168847,5.877946127808392,0.7676924463189478,0.41248629,"glycerophospholipid catabolic process"),
                     c("GO:0046503","glycerolipid catabolic process",0.1769338375994535,12.969119549206347,0.8436397456527538,0.65008118,"glycerophospholipid catabolic process"),
                     c("GO:0046505","sulfolipid metabolic process",0.17912499160755974,18.891686402549443,0.8339906611943675,0.56911536,"glycerophospholipid catabolic process"),
                     c("GO:0046506","sulfolipid biosynthetic process",0.17912006212722653,18.891686402549443,0.797953515482551,0.40555077,"glycerophospholipid catabolic process"),
                     c("GO:1903509","liposaccharide metabolic process",1.0088230796691766,7.9151382786097635,0.7962299131599007,0.46872312,"glycerophospholipid catabolic process"),
                     c("GO:0051262","protein tetramerization",0.04135094577502621,13.067648638754903,0.9845647187311816,0.00702396,"protein tetramerization"),
                     c("GO:0051259","protein complex oligomerization",0.23363025565172313,6.932708406644491,0.9835170301707383,0.53264918,"protein tetramerization"),
                     c("GO:0055081","monoatomic anion homeostasis",0.027077635470253197,20.447897610183766,0.9786192784256096,-0,"monoatomic anion homeostasis"),
                     c("GO:0030002","intracellular monoatomic anion homeostasis",0.0007492810106460015,9.7562611664949,0.9809243966867315,0.4346327,"monoatomic anion homeostasis"),
                     c("GO:0030643","intracellular phosphate ion homeostasis",0.04212487418733819,9.7562611664949,0.977417104647664,0.58356474,"monoatomic anion homeostasis"),
                     c("GO:0055062","phosphate ion homeostasis",0.04733779963969443,20.447897610183766,0.9777604007788775,0.51004304,"monoatomic anion homeostasis"),
                     c("GO:0072527","pyrimidine-containing compound metabolic process",0.9998119896200919,1.875163160149058,0.8527763343306781,0.06850747,"pyrimidine-containing compound metabolic process"),
                     c("GO:0018105","peptidyl-serine phosphorylation",0.0005496370571515077,1.854113144556883,0.8738302374777484,0.38053845,"pyrimidine-containing compound metabolic process"),
                     c("GO:0018209","peptidyl-serine modification",0.0017499655182850694,1.8270052194926965,0.9260011402844723,0.12122825,"pyrimidine-containing compound metabolic process"),
                     c("GO:1904966","positive regulation of vitamin E biosynthetic process",1.9717921332789515E-05,5.024448309977004,0.9003498039612475,0.08476606,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:0030656","regulation of vitamin metabolic process",0.004552375087707779,4.112646499057077,0.9172569616557398,0.50934064,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:0031330","negative regulation of cellular catabolic process",0.1353783183905996,1.388684722594566,0.9128020898925289,0.18202491,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:0043068","positive regulation of programmed cell death",0.20637269414930823,1.5594953414353876,0.8931994418351102,0.28065938,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:1904965","regulation of vitamin E biosynthetic process",1.9717921332789515E-05,5.024448309977004,0.9256858299224145,0.64749804,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:2000904","regulation of starch metabolic process",0.002077775960442695,2.021561581984971,0.9119338315137988,0.12069711,"positive regulation of vitamin E biosynthetic process"),
                     c("GO:2001006","regulation of cellulose biosynthetic process",0.002356291599268347,1.6050767418332463,0.9054472423106539,0.69507598,"positive regulation of vitamin E biosynthetic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)
}

# Repressed genes H20 under drought no treatment
repressed.genes.h2o.30.180.atha.enrich.go <- enrichGO(gene          = repressed.genes.h2o.30.180.atha,
                                                      OrgDb         = org.At.tair.db,
                                                      ont           = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.05,
                                                      readable      = FALSE,
                                                      keyType = "TAIR")

repressed.genes.h2o.30.180.atha.enrich.go.df <- as.data.frame(repressed.genes.h2o.30.180.atha.enrich.go)

repressed.genes.h2o.30.180.enrich.go <- enrichGO(gene = repressed.genes.h2o.30.180, 
                                                 OrgDb = org.Taestivum.eg.db,
                                                 universe = rownames(gene.expression.ta), 
                                                 ont           = "BP",
                                                 pAdjustMethod = "BH",
                                                 pvalueCutoff  = 0.05,
                                                 readable      = FALSE,
                                                 keyType = "GID")

repressed.genes.h2o.30.180.enrich.go.df <- as.data.frame(repressed.genes.h2o.30.180.enrich.go)


treeplot(x = pairwise_termsim(repressed.genes.h2o.30.180.enrich.go))

write.table(x = repressed.genes.h2o.30.180.enrich.go.df,file="repressed_genes_h2o_30_180_enrich_go.tsv",quote = F,sep = "\t",row.names = F)

{
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000122","negative regulation of transcription by RNA polymerase II",0.552067290955774,10.262242265116685,0.9406486988719428,0,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,1.8517034140920416,0.9457194576771062,0.28773891,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0010439","regulation of glucosinolate biosynthetic process",0.0001602081108289148,3.201191097494853,0.9238335769256238,0.1844076,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0042762","regulation of sulfur metabolic process",0.0025781182142622285,2.210747334427634,0.9570538858182825,0.20631318,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0043255","regulation of carbohydrate biosynthetic process",0.060962883280651976,1.9206717397674704,0.9305754649967788,0.63135342,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0043455","regulation of secondary metabolic process",0.06959440334408058,1.6526534393900982,0.9494158308365452,0.25380113,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0006094","gluconeogenesis",0.2574593635823993,2.210747334427634,0.7823681369431171,-0,"gluconeogenesis"),
                     c("GO:0009065","glutamine family amino acid catabolic process",0.20902475456856845,1.3586326572249678,0.8309778561963788,0.5727513,"gluconeogenesis"),
                     c("GO:0016143","S-glycoside metabolic process",0.005523482713347663,1.4442842042233532,0.870035272093434,0.64112133,"gluconeogenesis"),
                     c("GO:0019757","glycosinolate metabolic process",0.005523482713347663,1.4442842042233532,0.8417208924187091,0.19588594,"gluconeogenesis"),
                     c("GO:0019761","glucosinolate biosynthetic process",0.0008355469164769557,2.28928181528941,0.7616875461426931,0.57161592,"gluconeogenesis"),
                     c("GO:0030187","melatonin biosynthetic process",0.013573324097458981,1.3053373185857426,0.8728210811308253,0.33394023,"gluconeogenesis"),
                     c("GO:0033014","tetrapyrrole biosynthetic process",0.699579525186539,1.4452746063672393,0.8404429366453027,0.1698995,"gluconeogenesis"),
                     c("GO:0043649","dicarboxylic acid catabolic process",0.12580773232369688,1.3586326572249678,0.8507023797140587,0.2683099,"gluconeogenesis"),
                     c("GO:0051555","flavonol biosynthetic process",0.050921531841928915,1.3586326572249678,0.8122175867350676,0.13410507,"gluconeogenesis"),
                     c("GO:1901659","glycosyl compound biosynthetic process",0.32514359329736586,2.2005561169217653,0.8446417789997231,0.15537221,"gluconeogenesis"),
                     c("GO:0010438","cellular response to sulfur starvation",0.006469942937321558,3.304901536398002,0.9479800441880912,-0,"guard cell development"),
                     c("GO:0009640","photomorphogenesis",0.010953305300364573,2.2286830349886397,0.86215649133931,0.57361878,"guard cell development"),
                     c("GO:0009647","skotomorphogenesis",0.00023415031582687548,3.201191097494853,0.8697664378511087,0.12680683,"guard cell development"),
                     c("GO:0010167","response to nitrate",0.017245786945691028,2.130698255709449,0.970932212859384,0.15658045,"guard cell development"),
                     c("GO:0010441","guard cell development",0.00021196765432748727,1.6526534393900982,0.8777564430315863,0.66001247,"guard cell development"),
                     c("GO:0010442","guard cell morphogenesis",0.00015034915016252004,1.656840990974971,0.8790848358960972,0.69967204,"guard cell development"),
                     c("GO:0016036","cellular response to phosphate starvation",0.07367847780013462,1.6571018610896444,0.9442974599744555,0.69574796,"guard cell development"),
                     c("GO:0060776","simple leaf morphogenesis",0.0001158427878301384,1.6620587484381693,0.901947743292177,0.26941437,"guard cell development"),
                     c("GO:0015714","phosphoenolpyruvate transport",0.0052844029171875894,3.201191097494853,0.8605010063158842,-0,"phosphoenolpyruvate transport"),
                     c("GO:0015711","organic anion transport",1.9925107391193801,1.4856911601342944,0.8502741303609308,0.48423164,"phosphoenolpyruvate transport"),
                     c("GO:0015713","phosphoglycerate transmembrane transport",0.003985484849390081,3.158224376545195,0.8524787227129823,0.55136091,"phosphoenolpyruvate transport"),
                     c("GO:0015718","monocarboxylic acid transport",0.39675661883789076,2.1236755288983984,0.8323957405251643,0.51639335,"phosphoenolpyruvate transport"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,1.9223921184104689,0.8683820407439538,0.29818639,"phosphoenolpyruvate transport"),
                     c("GO:0015849","organic acid transport",1.7778688417111328,1.4442842042233532,0.8516094524770416,0.57946783,"phosphoenolpyruvate transport"),
                     c("GO:0034486","vacuolar transmembrane transport",0.025445977479964865,1.5607284578452274,0.9105281654249293,0.26214864,"phosphoenolpyruvate transport"),
                     c("GO:0035435","phosphate ion transmembrane transport",0.14767983656209366,1.3586326572249678,0.8831014219193966,0.66354944,"phosphoenolpyruvate transport"),
                     c("GO:0042873","aldonate transmembrane transport",0.023853755332342113,3.158224376545195,0.8373788388869667,0.60559877,"phosphoenolpyruvate transport"),
                     c("GO:0098656","monoatomic anion transmembrane transport",0.4548751919662879,2.2090165262963883,0.8914168833599501,0.18022361,"phosphoenolpyruvate transport"),
                     c("GO:1905011","transmembrane phosphate ion transport from cytosol to vacuole",0.0009292070428077058,1.7034829842002448,0.9089569288198147,0.21104468,"phosphoenolpyruvate transport"),
                     c("GO:0055081","monoatomic anion homeostasis",0.027077635470253197,1.3586326572249678,0.9808457036716208,-0,"monoatomic anion homeostasis"),
                     c("GO:0055062","phosphate ion homeostasis",0.04733779963969443,1.3586326572249678,0.9808457036716208,0.51004304,"monoatomic anion homeostasis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)
}
  
## Differential gene expression analysis of LRM treated plants under Drought Stress
## i.e LRM 30 vs H20 30
res.lrm.h2o.30 <- results(dds.wheat,contrast=c("treatment","LRM_30","H2O_30"))
res.lrm.h2o.30$padj[is.na(res.lrm.h2o.30$padj)] <- 1

genes.lrm.h2o.30 <- rownames(res.lrm.h2o.30)

activated.genes.lrm.h2o.30 <- genes.lrm.h2o.30[res.lrm.h2o.30$log2FoldChange > 1 & res.lrm.h2o.30$padj < 0.05]
length(activated.genes.lrm.h2o.30)

repressed.genes.lrm.h2o.30 <- genes.lrm.h2o.30[res.lrm.h2o.30$log2FoldChange < -1 & res.lrm.h2o.30$padj < 0.05]
length(repressed.genes.lrm.h2o.30)

logfc.lrm.h2o.30 <- res.lrm.h2o.30$log2FoldChange
names(logfc.lrm.h2o.30) <- genes.lrm.h2o.30
logqval.lrm.h2o.30 <- -log10(res.lrm.h2o.30$padj)
names(logqval.lrm.h2o.30) <- genes.lrm.h2o.30

plot(logfc.lrm.h2o.30[logqval.lrm.h2o.30 > 0],logqval.lrm.h2o.30[logqval.lrm.h2o.30 > 0],pch=19,cex=0.7,xlim=c(-12,12),ylim=c(0,35),col="grey",xlab="log2FC",ylab="-log10(q-value)")
points(logfc.lrm.h2o.30[activated.genes.lrm.h2o.30],logqval.lrm.h2o.30[activated.genes.lrm.h2o.30],cex=0.7,col="red",pch=19)
points(logfc.lrm.h2o.30[repressed.genes.lrm.h2o.30],logqval.lrm.h2o.30[repressed.genes.lrm.h2o.30],cex=0.7,col="blue",pch=19)

# ggplot2 volcano plot
# Loading relevant libraries 
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

df_lrm_30 <- as.data.frame(res.lrm.h2o.30)
head(df_lrm_30)

df_lrm_30$diffexpressed <- "NO"
df_lrm_30$diffexpressed[df_lrm_30$log2FoldChange > 1 & df_lrm_30$padj < 0.05] <- "UP"
df_lrm_30$diffexpressed[df_lrm_30$log2FoldChange < -1 & df_lrm_30$padj < 0.05] <- "DOWN"
df_lrm_30$names <- rownames(df_lrm_30)

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

volcano_plot_lrm_30 <-  ggplot(data = df_lrm_30, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') + xlim(c(-10,10)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />

plot(volcano_plot_lrm_30)
png("volcano_lrm_30.png", width = 3000, height = 2000, res = 300)
plot(volcano_plot_lrm_30)
dev.off()

# GO enrichment activated genes (LRM 30 vs H20 30)
activated.genes.lrm.h2o.30.enrich.go <- enrichGO(gene = activated.genes.lrm.h2o.30, 
                                                 OrgDb = org.Taestivum.eg.db,
                                                 universe = rownames(gene.expression.ta), 
                                                 ont           = "BP",
                                                 pAdjustMethod = "BH",
                                                 pvalueCutoff  = 0.05,
                                                 readable      = FALSE,
                                                 keyType = "GID")

activated.genes.lrm.h2o.30.enrich.go.df <- as.data.frame(activated.genes.lrm.h2o.30.enrich.go)

write.table(x = activated.genes.lrm.h2o.30.enrich.go.df,file="activated_genes_lrm_h2o_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)

# Treemap
{
revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005982","starch metabolic process",0.040983699490203,0.9010512635831647,0.07962189,"starch metabolic process"),
                     c("GO:0005984","disaccharide metabolic process",0.2026632601985772,0.8857305302732708,0.47009629,"starch metabolic process"),
                     c("GO:0009311","oligosaccharide metabolic process",0.3704479822996164,0.8931154800567558,0.414014,"starch metabolic process"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,0.94250869059614,-0,"regulation of carbohydrate metabolic process"),
                     c("GO:0030656","regulation of vitamin metabolic process",0.004552375087707779,0.9348806505492289,0.19115702,"regulation of carbohydrate metabolic process"),
                     c("GO:0032107","regulation of response to nutrient levels",0.07313130548314971,0.9074508710231366,0.58176449,"regulation of carbohydrate metabolic process"),
                     c("GO:0032881","regulation of polysaccharide metabolic process",0.057874563851903815,0.8942708333752016,0.24142817,"regulation of carbohydrate metabolic process"),
                     c("GO:0080148","negative regulation of response to water deprivation",0.0001774612919951056,0.9337815898885051,0.42513778,"regulation of carbohydrate metabolic process"),
                     c("GO:1900367","positive regulation of defense response to insect",0.00017006707149530955,0.9174922904813134,0.42440006,"regulation of carbohydrate metabolic process"),
                     c("GO:1900425","negative regulation of defense response to bacterium",0.0013432833907962855,0.9155837914566843,0.52914681,"regulation of carbohydrate metabolic process"),
                     c("GO:1901141","regulation of lignin biosynthetic process",0.00025140349699306625,0.9419100948493303,0.64033671,"regulation of carbohydrate metabolic process"),
                     c("GO:1904965","regulation of vitamin E biosynthetic process",1.9717921332789515E-05,0.9387742155369267,0.64749804,"regulation of carbohydrate metabolic process"),
                     c("GO:1904966","positive regulation of vitamin E biosynthetic process",1.9717921332789515E-05,0.93207958286999,0.50934064,"regulation of carbohydrate metabolic process"),
                     c("GO:1905034","regulation of antifungal innate immune response",0.00037217576515640203,0.9160088444501951,0.48673677,"regulation of carbohydrate metabolic process"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,0.9282321045735893,0.10790585,"regulation of carbohydrate metabolic process"),
                     c("GO:2000762","regulation of phenylpropanoid metabolic process",0.005555524335513445,0.9347617055354593,0.19348436,"regulation of carbohydrate metabolic process"),
                     c("GO:0006767","water-soluble vitamin metabolic process",1.2065050999910578,0.8283392366983222,0.03842157,"water-soluble vitamin metabolic process"),
                     c("GO:0006563","L-serine metabolic process",0.2663767935051533,0.8193015915611126,0.2985208,"water-soluble vitamin metabolic process"),
                     c("GO:0006766","vitamin metabolic process",1.256768546208505,0.8498242814287484,0.35118373,"water-soluble vitamin metabolic process"),
                     c("GO:0009404","toxin metabolic process",0.08737257416575693,0.9228528451525968,0.05058592,"toxin metabolic process"),
                     c("GO:0009407","toxin catabolic process",0.009016019529418004,0.8529198286216787,0.59727878,"toxin metabolic process"),
                     c("GO:0009809","lignin biosynthetic process",0.024309732263162874,0.8270574608808813,0.6358361,"toxin metabolic process"),
                     c("GO:0009625","response to insect",0.000539778096485113,0.9643009689825971,-0,"response to insect"),
                     c("GO:0009749","response to glucose",0.058823488816044316,0.9094081045506305,0.51531967,"response to insect"),
                     c("GO:0016036","cellular response to phosphate starvation",0.07367847780013462,0.9491532934950722,0.37275592,"response to insect"),
                     c("GO:0043434","response to peptide hormone",0.19767216136121488,0.9039099314704905,0.69760009,"response to insect"),
                     c("GO:0046686","response to cadmium ion",0.02872408190154112,0.9378122200782615,0.14134599,"response to insect"),
                     c("GO:0060359","response to ammonium ion",0.0006901272466476329,0.9375770671137353,0.49996243,"response to insect"),
                     c("GO:0061760","antifungal innate immune response",0.0015700394861233652,0.9626447014747748,0.39537663,"response to insect"),
                     c("GO:0071369","cellular response to ethylene stimulus",0.038408046016107374,0.9164104648674045,0.35931117,"response to insect"),
                     c("GO:0009718","anthocyanin-containing compound biosynthetic process",0.0029700119007514203,0.9020842985252061,0.03652024,"anthocyanin-containing compound biosynthetic process"),
                     c("GO:0006011","UDP-glucose metabolic process",0.04799835000434287,0.7594673361208526,0.33285723,"anthocyanin-containing compound biosynthetic process"),
                     c("GO:0046349","amino sugar biosynthetic process",0.16032395361674495,0.8232211988786658,0.10773451,"anthocyanin-containing compound biosynthetic process"),
                     c("GO:0010230","alternative respiration",0.007818155808451042,0.9417187154417953,0.02717564,"alternative respiration"),
                     c("GO:0010431","seed maturation",0.002551006072429643,0.9908235224396884,0,"seed maturation"),
                     c("GO:0009704","de-etiolation",0.0020876349211090897,0.9609560445563063,0.59324692,"seed maturation"),
                     c("GO:0015977","carbon fixation",0.042674511244489705,0.9506694176650013,-0,"carbon fixation"),
                     c("GO:0017014","protein nitrosylation",9.612486649734888E-05,0.9333055858469922,0.03072367,"protein nitrosylation"),
                     c("GO:0006575","cellular modified amino acid metabolic process",1.1845763067288293,0.8899084982258868,0.10312018,"protein nitrosylation"),
                     c("GO:0006749","glutathione metabolic process",0.33273992249082307,0.863833746073326,0.18598634,"protein nitrosylation"),
                     c("GO:0006750","glutathione biosynthetic process",0.05906996283270419,0.8239761731621957,0.6671716,"protein nitrosylation"),
                     c("GO:0018119","peptidyl-cysteine S-nitrosylation",9.119538616415149E-05,0.9318198163087279,0.38938339,"protein nitrosylation"),
                     c("GO:0018198","peptidyl-cysteine modification",0.01781021244384213,0.9111242942037798,0.28022578,"protein nitrosylation"),
                     c("GO:0019184","nonribosomal peptide biosynthetic process",0.19618345830058925,0.8335428267346694,0.46794114,"protein nitrosylation"),
                     c("GO:0046505","sulfolipid metabolic process",0.17912499160755974,0.8380347877586894,0.59841245,"protein nitrosylation"),
                     c("GO:0046506","sulfolipid biosynthetic process",0.17912006212722653,0.8002245494126037,0.59841108,"protein nitrosylation"),
                     c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,0.8779515081200668,0.07391108,"monoterpene biosynthetic process"),
                     c("GO:0006643","membrane lipid metabolic process",1.10151702785462,0.8317855568036494,0.30215704,"monoterpene biosynthetic process"),
                     c("GO:0006694","steroid biosynthetic process",0.2949431320360321,0.7495016244081962,0.49532952,"monoterpene biosynthetic process"),
                     c("GO:0009247","glycolipid biosynthetic process",0.7048811812848927,0.7316725200417646,0.59578979,"monoterpene biosynthetic process"),
                     c("GO:0009395","phospholipid catabolic process",0.12197752610480252,0.7615761089655728,0.69336156,"monoterpene biosynthetic process"),
                     c("GO:0016042","lipid catabolic process",1.4093137798594644,0.8239087852360161,0.58366296,"monoterpene biosynthetic process"),
                     c("GO:0043692","monoterpene metabolic process",8.380116566435544E-05,0.9020535856788482,0.50228334,"monoterpene biosynthetic process"),
                     c("GO:0046475","glycerophospholipid catabolic process",0.05117540007908858,0.7747625911416651,0.65008118,"monoterpene biosynthetic process"),
                     c("GO:0046503","glycerolipid catabolic process",0.1769338375994535,0.830597232201768,0.5198576,"monoterpene biosynthetic process"),
                     c("GO:1903509","liposaccharide metabolic process",1.0088230796691766,0.7846528490840348,0.56223164,"monoterpene biosynthetic process"),
                     c("GO:0051259","protein complex oligomerization",0.23363025565172313,0.9796462491856014,-0,"protein complex oligomerization"),
                     c("GO:0051262","protein tetramerization",0.04135094577502621,0.9810813081466744,0.53264918,"protein complex oligomerization"),
                     c("GO:0055081","monoatomic anion homeostasis",0.027077635470253197,0.9776304492940234,-0,"monoatomic anion homeostasis"),
                     c("GO:0030002","intracellular monoatomic anion homeostasis",0.0007492810106460015,0.980041630011044,0.4346327,"monoatomic anion homeostasis"),
                     c("GO:0030643","intracellular phosphate ion homeostasis",0.04212487418733819,0.9763730144016386,0.50632071,"monoatomic anion homeostasis"),
                     c("GO:0055062","phosphate ion homeostasis",0.04733779963969443,0.9767320873423171,0.58356474,"monoatomic anion homeostasis"),
                     c("GO:0071836","nectar secretion",0.00012077226816333578,0.978056720368052,-0,"nectar secretion"),
                     c("GO:0006817","phosphate ion transport",0.2426906405041399,0.9670929941021582,0.13980114,"nectar secretion"),
                     c("GO:0015713","phosphoglycerate transmembrane transport",0.003985484849390081,0.9416117155462926,0.55136091,"nectar secretion"),
                     c("GO:0015714","phosphoenolpyruvate transport",0.0052844029171875894,0.9464009226641851,0.25038389,"nectar secretion"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,0.9475198029600548,0.32578074,"nectar secretion"),
                     c("GO:0035627","ceramide transport",0.019542924780961007,0.9495592565263469,0.18863509,"nectar secretion"),
                     c("GO:0042873","aldonate transmembrane transport",0.023853755332342113,0.9361795865988178,0.60559877,"nectar secretion"),
                     c("GO:0090481","pyrimidine nucleotide-sugar transmembrane transport",0.05787702859207042,0.9307611543532037,0.35289052,"nectar secretion"),
                     c("GO:0120009","intermembrane lipid transfer",0.11033409355779031,0.9320938834326131,0.60157124,"nectar secretion"),
                     c("GO:1901264","carbohydrate derivative transport",0.41917836013343907,0.9471295251252203,0.4160251,"nectar secretion"),
                     c("GO:1901334","lactone metabolic process",0.08600464337329466,0.8884992236195826,0.05918817,"lactone metabolic process"),
                     c("GO:0008655","pyrimidine-containing compound salvage",0.07123099081470212,0.7470027091689898,0.31123743,"lactone metabolic process"),
                     c("GO:0009173","pyrimidine ribonucleoside monophosphate metabolic process",0.21775979371899418,0.7371605672875675,0.47971038,"lactone metabolic process"),
                     c("GO:0009208","pyrimidine ribonucleoside triphosphate metabolic process",0.16502421311444862,0.7426340610988151,0.45622658,"lactone metabolic process"),
                     c("GO:0019336","phenol-containing compound catabolic process",0.04817581129633798,0.843080780469754,0.15525694,"lactone metabolic process"),
                     c("GO:0032958","inositol phosphate biosynthetic process",0.03756017539879743,0.7770140263632731,0.67851077,"lactone metabolic process"),
                     c("GO:0043647","inositol phosphate metabolic process",0.05878158823321214,0.8034227702582816,0.50083657,"lactone metabolic process"),
                     c("GO:0046244","salicylic acid catabolic process",2.21826614993882E-05,0.8689139087482665,0.51820171,"lactone metabolic process"),
                     c("GO:0046434","organophosphate catabolic process",0.9710632603168847,0.802895435030072,0.35595493,"lactone metabolic process"),
                     c("GO:1901336","lactone biosynthetic process",0.05878651771354534,0.8134133279504976,0.16682536,"lactone metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

library(treemap)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  position.legend = "none"
)

}

# GO enrichment repressed genes (LRM 30 vs H20 30)
repressed.genes.lrm.h2o.30.enrich.go <- enrichGO(gene = repressed.genes.lrm.h2o.30, 
                                                 OrgDb = org.Taestivum.eg.db,
                                                 universe = rownames(gene.expression.ta), 
                                                 ont           = "BP",
                                                 pAdjustMethod = "BH",
                                                 pvalueCutoff  = 0.05,
                                                 readable      = FALSE,
                                                 keyType = "GID")

repressed.genes.lrm.h2o.30.enrich.go.df <- as.data.frame(repressed.genes.lrm.h2o.30.enrich.go)

write.table(x = repressed.genes.lrm.h2o.30.enrich.go.df,file="repressed_genes_lrm_h2o_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)

{
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000103","sulfate assimilation",0.10156701278519879,1.6222524673914884,0.9192600516952587,0.06990154,"sulfate assimilation"),
                     c("GO:0019762","glucosinolate catabolic process",0.002250307772104603,1.4654779848950146,0.8245535545984743,0.40918029,"sulfate assimilation"),
                     c("GO:0000122","negative regulation of transcription by RNA polymerase II",0.552067290955774,14.727254539059688,0.8636977237090934,-0,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0000957","mitochondrial RNA catabolic process",0.004922086112697582,3.7593699017577746,0.8403657414912994,0.66263955,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0000958","mitochondrial mRNA catabolic process",0.004781595923201457,3.7593699017577746,0.7451378185285054,0.46886586,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0006109","regulation of carbohydrate metabolic process",0.1417866428237562,6.235778772857743,0.9109584525548557,0.28773891,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0008360","regulation of cell shape",0.6483080002409529,1.4654779848950146,0.9190655248579425,0.59613965,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0009895","negative regulation of catabolic process",0.16801640767669945,3.51270300897809,0.8681468226816202,0.68477829,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0010322","regulation of isopentenyl diphosphate biosynthetic process, methylerythritol 4-phosphate pathway",0.0001355607091629279,7.888848341032003,0.8442230175143017,0.18081773,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0010906","regulation of glucose metabolic process",0.04743392450619177,5.08094493867497,0.8599618401929324,0.61888491,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0016556","mRNA modification",0.09789208519680015,1.3183008510674807,0.8358874090701082,0.46888121,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0019216","regulation of lipid metabolic process",0.10650635207906255,2.8904409432094935,0.9124880507464505,0.28090334,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0019747","regulation of isoprenoid metabolic process",0.002585512434762025,3.942316547809627,0.8897385864265288,0.68089704,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0033238","regulation of amine metabolic process",0.01681199267636966,1.5035002169942018,0.9218218100941408,0.24618673,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0045995","regulation of embryonic development",0.01473914619626016,2.9796029472354633,0.9447609628632405,0.13698677,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0050821","protein stabilization",0.10143638155636905,1.3037680057003005,0.9347535092094796,0.15044246,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0090352","regulation of nitrate assimilation",0.00011091330749694102,2.8994319347436623,0.8989267781905774,0.47591064,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1900367","positive regulation of defense response to insect",0.00017006707149530955,1.6347361473891968,0.9375881584629424,0.48673677,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1902326","positive regulation of chlorophyll biosynthetic process",0.00019717921332789513,4.553892109004702,0.9057701919318755,0.18649588,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1902369","negative regulation of RNA catabolic process",0.07210843831401124,6.026938095453567,0.8649427170803534,0.65315277,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1903314","regulation of nitrogen cycle metabolic process",0.00024400927649327026,2.8994319347436623,0.9363949670267521,0.18868876,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1903725","regulation of phospholipid metabolic process",0.0073326019956311,6.178021906641218,0.8734486485950375,0.59373684,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1904965","regulation of vitamin E biosynthetic process",1.9717921332789515E-05,2.3654122038121663,0.900556031352932,0.44751716,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:1905034","regulation of antifungal innate immune response",0.00037217576515640203,2.5055361963165494,0.9415686005489895,0.10851536,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:2000068","regulation of defense response to insect",0.00025140349699306625,1.5569497246385358,0.9453142611799386,0.49422838,"negative regulation of transcription by RNA polymerase II"),
                     c("GO:0006081","cellular aldehyde metabolic process",0.642313752155785,2.4316003561251254,0.9189159536251771,0.06162237,"cellular aldehyde metabolic process"),
                     c("GO:0006775","fat-soluble vitamin metabolic process",0.050263446217447064,7.368891268808366,0.8554361564660422,0.02907341,"fat-soluble vitamin metabolic process"),
                     c("GO:0006006","glucose metabolic process",0.4804074353520837,6.235778772857743,0.8154858742767039,0.23877753,"fat-soluble vitamin metabolic process"),
                     c("GO:0006020","inositol metabolic process",0.11004325421813167,1.538012172890235,0.8602323737404128,0.25417668,"fat-soluble vitamin metabolic process"),
                     c("GO:0006580","ethanolamine metabolic process",0.013943035122448785,1.5032319026109122,0.8485949062979806,0.53767207,"fat-soluble vitamin metabolic process"),
                     c("GO:0006766","vitamin metabolic process",1.256768546208505,3.14835551417746,0.8457384694714511,0.31789605,"fat-soluble vitamin metabolic process"),
                     c("GO:0010236","plastoquinone biosynthetic process",0.024430504531326204,6.322417088469348,0.7871181970197374,0.68109931,"fat-soluble vitamin metabolic process"),
                     c("GO:0019520","aldonic acid metabolic process",0.06686840071982243,2.9450993037092212,0.8400290087221827,0.43854599,"fat-soluble vitamin metabolic process"),
                     c("GO:0042372","phylloquinone biosynthetic process",0.033601802691239926,6.989454615885736,0.7619374016380759,0.55426757,"fat-soluble vitamin metabolic process"),
                     c("GO:0046177","D-gluconate catabolic process",0.041346016294693005,2.9450993037092212,0.8251261358106107,0.42248785,"fat-soluble vitamin metabolic process"),
                     c("GO:0009768","photosynthesis, light harvesting in photosystem I",0.017839789325841314,14.132822804546084,0.8918970910719564,0.0487292,"photosynthesis, light harvesting in photosystem I"),
                     c("GO:0009051","pentose-phosphate shunt, oxidative branch",0.054428857098998855,2.7545452823987726,0.7420119510787526,0.44486773,"photosynthesis, light harvesting in photosystem I"),
                     c("GO:0009112","nucleobase metabolic process",0.6825752827771746,1.6208434346010863,0.748668781519601,0.56688699,"photosynthesis, light harvesting in photosystem I"),
                     c("GO:0009123","nucleoside monophosphate metabolic process",1.15699586426459,1.6820613278853145,0.7216552231696831,0.65136519,"photosynthesis, light harvesting in photosystem I"),
                     c("GO:0009156","ribonucleoside monophosphate biosynthetic process",0.8452407395521883,2.3574741911489747,0.682859919300484,0.48935317,"photosynthesis, light harvesting in photosystem I"),
                     c("GO:0016036","cellular response to phosphate starvation",0.07367847780013462,3.51270300897809,0.9659366662425385,0.0067085,"cellular response to phosphate starvation"),
                     c("GO:0000304","response to singlet oxygen",0.002134464984274465,2.280550066908692,0.9659240626466802,0.54752824,"cellular response to phosphate starvation"),
                     c("GO:0009647","skotomorphogenesis",0.00023415031582687548,2.162069459950188,0.97532750992334,0.14209398,"cellular response to phosphate starvation"),
                     c("GO:0010167","response to nitrate",0.017245786945691028,2.9796029472354633,0.9692867822664131,0.18056805,"cellular response to phosphate starvation"),
                     c("GO:0010196","nonphotochemical quenching",0.0005693549784842973,1.3551572464477073,0.9788698051810433,0.47447992,"cellular response to phosphate starvation"),
                     c("GO:0010343","singlet oxygen-mediated programmed cell death",0.0019717921332789512,2.7817392539681842,0.9530119084256998,0.51379914,"cellular response to phosphate starvation"),
                     c("GO:0048598","embryonic morphogenesis",0.16375240718848372,1.5035002169942018,0.9958206574837686,0.37025908,"cellular response to phosphate starvation"),
                     c("GO:0061760","antifungal innate immune response",0.0015700394861233652,2.5055361963165494,0.9795523399182661,0.39453375,"cellular response to phosphate starvation"),
                     c("GO:0097468","programmed cell death in response to reactive oxygen species",0.002447486985432498,2.559328741250107,0.9575456075572536,0.61430271,"cellular response to phosphate starvation"),
                     c("GO:1990066","energy quenching",0.000579213939150692,1.3551572464477073,0.9788581489335669,0.49228548,"cellular response to phosphate starvation"),
                     c("GO:0019695","choline metabolic process",0.06108858502914852,3.3762819021020554,0.924730980254578,0.07671199,"choline metabolic process"),
                     c("GO:0019750","chloroplast localization",0.009112144395915353,3.560301425816415,0.9860018817870662,-0,"chloroplast localization"),
                     c("GO:0045036","protein targeting to chloroplast",0.031610292636628186,2.4264499883406616,0.9532492125217118,0.5362014,"chloroplast localization"),
                     c("GO:0051644","plastid localization",0.009161439199247327,3.560301425816415,0.9859994308177215,0.51936881,"chloroplast localization"),
                     c("GO:0072598","protein localization to chloroplast",0.06766944127396701,2.3323060521986805,0.965402431063719,0.15675683,"chloroplast localization"),
                     c("GO:0071941","nitrogen cycle metabolic process",0.2034347238707226,3.438427458462438,0.923688699279676,0.08380609,"nitrogen cycle metabolic process"),
                     c("GO:0080144","intracellular amino acid homeostasis",0.0009045596411417188,1.4654779848950146,1,-0,"intracellular amino acid homeostasis"),
                     c("GO:0097164","ammonium ion metabolic process",0.08456030563566783,3.3762819021020554,0.9283971050996583,0.0785081,"ammonium ion metabolic process"),
                     c("GO:1901259","chloroplast rRNA processing",0.02042037228027014,24.620675620889035,0.8059582042700874,0,"chloroplast rRNA processing"),
                     c("GO:0000373","Group II intron splicing",0.009405448475740597,2.24959109513692,0.8316362208724755,0.33478559,"chloroplast rRNA processing"),
                     c("GO:0000959","mitochondrial RNA metabolic process",0.13719976137371603,1.5035002169942018,0.8360334776870517,0.30536567,"chloroplast rRNA processing"),
                     c("GO:0006353","DNA-templated transcription termination",0.20992684946954357,2.4264499883406616,0.7714274054287262,0.4662656,"chloroplast rRNA processing"),
                     c("GO:0006414","translational elongation",0.46125640425761183,1.417425968504699,0.7856024217952887,0.45686782,"chloroplast rRNA processing"),
                     c("GO:0006426","glycyl-tRNA aminoacylation",0.045033267583924654,6.902693829041655,0.7860092020565309,0.37539321,"chloroplast rRNA processing"),
                     c("GO:0006438","valyl-tRNA aminoacylation",0.03575598559684719,2.0334187907846415,0.789099736371623,0.61061417,"chloroplast rRNA processing"),
                     c("GO:0006656","phosphatidylcholine biosynthetic process",0.04863918244765853,4.348359560846457,0.7887091601323529,0.26737914,"chloroplast rRNA processing"),
                     c("GO:0006809","nitric oxide biosynthetic process",0.015387372860075615,3.877718112383114,0.8363619108373285,0.25888136,"chloroplast rRNA processing"),
                     c("GO:0008655","pyrimidine-containing compound salvage",0.07123099081470212,2.4264499883406616,0.7903772821393771,0.69069611,"chloroplast rRNA processing"),
                     c("GO:0009308","amine metabolic process",0.759497358636553,2.216805924891084,0.8805030010558187,0.1916293,"chloroplast rRNA processing"),
                     c("GO:0009828","plant-type cell wall loosening",4.6830063165375095E-05,1.5035002169942018,0.9806192672736802,0.14323902,"chloroplast rRNA processing"),
                     c("GO:0010239","chloroplast mRNA processing",0.0012742706661315222,2.8994319347436623,0.8456724993147391,0.30093069,"chloroplast rRNA processing"),
                     c("GO:0016120","carotene biosynthetic process",0.045947686185732764,2.9450993037092212,0.8260858420570822,0.41469173,"chloroplast rRNA processing"),
                     c("GO:0016123","xanthophyll biosynthetic process",0.02847021366438146,2.507151579954429,0.7938047952698833,0.6584158,"chloroplast rRNA processing"),
                     c("GO:0016485","protein processing",0.5034847975319472,2.789830248581815,0.792478451588303,0.69809989,"chloroplast rRNA processing"),
                     c("GO:0016540","protein autoprocessing",0.02521675664447119,9.209600435060736,0.8306970141718272,0.21508806,"chloroplast rRNA processing"),
                     c("GO:0019682","glyceraldehyde-3-phosphate metabolic process",0.26501132745285766,3.877718112383114,0.8600171901390183,0.37144194,"chloroplast rRNA processing"),
                     c("GO:0019860","uracil metabolic process",0.03451375655288144,3.560301425816415,0.7786561213551222,0.67115554,"chloroplast rRNA processing"),
                     c("GO:0031425","chloroplast RNA processing",0.003108037350080947,4.956100392801173,0.8418294872906436,0.31514581,"chloroplast rRNA processing"),
                     c("GO:0032544","plastid translation",0.006610433126817685,14.132822804546084,0.8034157166166224,0.24388633,"chloroplast rRNA processing"),
                     c("GO:0033014","tetrapyrrole biosynthetic process",0.699579525186539,14.727254539059688,0.7501108960062173,0.17970027,"chloroplast rRNA processing"),
                     c("GO:0042793","plastid transcription",0.0037932351163953828,4.893050545474161,0.8218923700270581,0.30582428,"chloroplast rRNA processing"),
                     c("GO:0043038","amino acid activation",0.9712111447268806,3.1957388194256136,0.9030635655262761,0.44853605,"chloroplast rRNA processing"),
                     c("GO:0043100","pyrimidine nucleobase salvage",0.01789894308983968,3.560301425816415,0.7541633428297179,0.34073754,"chloroplast rRNA processing"),
                     c("GO:0043572","plastid fission",0.035997530133173854,3.10801014767406,0.9651885676476899,0.28579969,"chloroplast rRNA processing"),
                     c("GO:0046470","phosphatidylcholine metabolic process",0.07026974214972863,3.628752333463245,0.8167500185825053,0.60490362,"chloroplast rRNA processing"),
                     c("GO:0046490","isopentenyl diphosphate metabolic process",0.18732764688200018,2.8352902343064144,0.8363552489998127,0.56664441,"chloroplast rRNA processing"),
                     c("GO:0051083","'de novo' cotranslational protein folding",0.03352046626574218,2.4264499883406616,0.8192571519122461,0.63022307,"chloroplast rRNA processing"),
                     c("GO:0051290","protein heterotetramerization",0.0013851839736284633,1.5032319026109122,0.9765595991657237,0.25189216,"chloroplast rRNA processing"),
                     c("GO:0051604","protein maturation",1.8440323267433083,1.5585602350127292,0.7973503615045632,0.37152587,"chloroplast rRNA processing"),
                     c("GO:0061077","chaperone-mediated protein folding",0.2637543099678923,5.3027915299625,0.7982570504766758,0.55114814,"chloroplast rRNA processing"),
                     c("GO:0140053","mitochondrial gene expression",0.3272164397774754,1.7769498184419086,0.8400883402414306,0.25768409,"chloroplast rRNA processing"),
                     c("GO:1904821","chloroplast disassembly",0.00021196765432748727,3.2105079505142533,0.9589661874239819,0.55339525,"chloroplast rRNA processing"),
                     c("GO:1901668","regulation of superoxide dismutase activity",1.9717921332789515E-05,1.5032319026109122,0.9566512490307988,0.08814938,"regulation of superoxide dismutase activity"),
                     c("GO:2001057","reactive nitrogen species metabolic process",0.11012952012396263,5.705403823246963,0.9270425106750579,0.0712206,"reactive nitrogen species metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
library(treemap)
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)
}

## Comparison between DEGs in LRM treated plants under Drought Stress and DEGs between Drought and Full Irrigation
length(intersect(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180))
length(intersect(activated.genes.lrm.h2o.30,repressed.genes.h2o.30.180))

length(intersect(repressed.genes.lrm.h2o.30,repressed.genes.h2o.30.180))
length(intersect(repressed.genes.lrm.h2o.30,activated.genes.h2o.30.180))

{
library(VennDiagram)

png("venn_activated.png", width = 1800, height = 1800, res = 300)
grid.newpage()
draw.pairwise.venn(area1 = length(activated.genes.h2o.30.180),cat.pos = c(190,180),cat.dist = 0.05,
                   area2 = length(activated.genes.lrm.h2o.30),
                   cross.area = length(intersect(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)),
                   lwd = 3,category = c("Drought","LRM"),euler.d = T,
                   col = c("darkorange","darkgreen"),
                   fill = c("darkorange","darkgreen"),alpha = 0.5, cat.col = c("darkorange","darkgreen"),
                   cex = 2,
                   cat.cex = 2)
dev.off()

png("venn_repressed.png", width = 1800, height = 1800, res = 300)
grid.newpage()
draw.pairwise.venn(area1 = length(repressed.genes.h2o.30.180),cat.pos = c(190,180),cat.dist = 0.05,
                   area2 = length(repressed.genes.lrm.h2o.30),
                   cross.area = length(intersect(repressed.genes.h2o.30.180,repressed.genes.lrm.h2o.30)),
                   lwd = 3,category = c("Drought","LRM"),euler.d = T,
                   col = c("darkorange","darkgreen"),
                   fill = c("darkorange","darkgreen"),alpha = 0.5, cat.col = c("darkorange","darkgreen"),
                   cex = 2,
                   cat.cex = 2)
dev.off()

}

sum(logfc.lrm.h2o.30[setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)] > 0)/length(setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30))
sum(logfc.lrm.h2o.30[setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)] < 0)/length(setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30))


sum(logfc.h2o.30.180[setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)] > 0)/length(setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180))
sum(logfc.h2o.30.180[setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)] < 0)/length(setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180))


col2rgb("darkorange")
transparent.darkorange <- rgb(red = 255,green = 140,blue = 0,max=255,alpha=150)
col2rgb("darkgreen")
transparent.darkgreen <- rgb(red = 0,green = 100,blue = 0,max=255,alpha=150)


boxplot(logfc.lrm.h2o.30[setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)],
        logfc.h2o.30.180[setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)],
        logfc.lrm.h2o.30[intersect(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)],
        logfc.h2o.30.180[intersect(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)],
        outline=F,
        col=c(transparent.darkorange,transparent.darkgreen,colorRampPalette(c(transparent.darkorange,transparent.darkgreen))(5)[2],colorRampPalette(c(transparent.darkorange,transparent.darkgreen))(5)[2]))
lines(x = c(0,10),y=c(0,0),lty=3)

# Make it a raincloud plot
# Rainclouds

library(ggdist)
library(ggplot2)
library(MetBrewer)
library(ggrepel)
library(cowplot)
library(RColorBrewer)

# Prepare matrix
first.rain <- data.frame(group="1", value=as.numeric(logfc.lrm.h2o.30[setdiff(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)]))
second.rain <- data.frame(group="2", value=as.numeric(logfc.h2o.30.180[setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)]))
third.rain <- data.frame(group="3", value=as.numeric(logfc.lrm.h2o.30[intersect(activated.genes.h2o.30.180,activated.genes.lrm.h2o.30)]))
fourth.rain <- data.frame(group="4", value=as.numeric(logfc.h2o.30.180[intersect(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)]))

cloud.data <- rbind(first.rain, second.rain, third.rain, fourth.rain)
head(cloud.data)

colors.cloud <- c(transparent.darkorange,transparent.darkgreen,colorRampPalette(c(transparent.darkorange,transparent.darkgreen))(5)[2],colorRampPalette(c(transparent.darkorange,transparent.darkgreen))(5)[2])
 
rainplot <- ggplot(as.data.frame(cloud.data), aes(x = factor(group, levels=unique(group)), y = value, fill = group)) +
  scale_fill_manual(values=colors.cloud, breaks = unique(cloud.data$group))+
  scale_color_manual(values=colors.cloud, breaks = unique(cloud.data$group))+
  ggdist::stat_halfeye(
    adjust = .35, 
    width = .5, 
    .width = 0, 
    justification = -.45, 
    point_colour = NA,
  ) + 
  stat_boxplot(aes(),geom = 'errorbar',linetype=1, width=0.1) +
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    lwd=1
  ) +
  geom_point(
    size = 1.0, aes(color = group),
    alpha = .1,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(xlim = c(1.0, NA), clip = "off")+
  
  xlab(NULL) +
  ylab(NULL)+
  geom_hline(yintercept = 0, col = "black", linetype = 'dashed') +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()) +
  ylim(-2,5.5)

png("raincloud.png", width = 1800, height = 1800, res = 300)
plot(rainplot)
dev.off()


#ggsave2("raincloud.png",width = 15,height = 10)


# Functional enrichment of the comparison between LRM and drought activated
activated.specific.lrm.30 <- setdiff(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)
length(activated.specific.lrm.30)

activated.specific.lrm.30.enrich.go <- enrichGO(gene = activated.specific.lrm.30, 
                                                OrgDb = org.Taestivum.eg.db,
                                                universe = rownames(gene.expression.ta), 
                                                ont           = "BP",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                readable      = FALSE,
                                                keyType = "GID")


activated.specific.lrm.30.enrich.go.df <- as.data.frame(activated.specific.lrm.30.enrich.go)

write.table(x = activated.specific.lrm.30.enrich.go.df,file="activated_specific_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)

{
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0009404","toxin metabolic process",0.08737257416575693,2.8895743801233826,0.8847610930577874,0.05058592,"toxin metabolic process"),
                     c("GO:0009407","toxin catabolic process",0.009016019529418004,3.031670009044683,0.8031732921038727,0.59727878,"toxin metabolic process"),
                     c("GO:0009809","lignin biosynthetic process",0.024309732263162874,1.4266081496120182,0.8345624487791715,0.6358361,"toxin metabolic process"),
                     c("GO:0010230","alternative respiration",0.007818155808451042,3.299451841774647,0.9551550271918405,0.02006951,"alternative respiration"),
                     c("GO:0015977","carbon fixation",0.042674511244489705,10.182564354913458,0.9580037276457208,0.03072367,"carbon fixation"),
                     c("GO:0016126","sterol biosynthetic process",0.20334106374439187,2.3117090363670103,0.8366403134124812,0.04757864,"sterol biosynthetic process"),
                     c("GO:0019852","L-ascorbic acid metabolic process",0.0667599521524921,1.609498729889624,0.9177617991330277,0.16889498,"sterol biosynthetic process"),
                     c("GO:0043692","monoterpene metabolic process",8.380116566435544E-05,1.3112354007929257,0.9076672622802958,0.50228334,"sterol biosynthetic process"),
                     c("GO:0043693","monoterpene biosynthetic process",7.147746483136198E-05,1.3227323141469884,0.8903673763070111,0.31207673,"sterol biosynthetic process"),
                     c("GO:1901334","lactone metabolic process",0.08600464337329466,1.360345270214386,0.9183982035493109,0.17212211,"sterol biosynthetic process"),
                     c("GO:0017014","protein nitrosylation",9.612486649734888E-05,10.182564354913458,0.9231940750001886,0,"protein nitrosylation"),
                     c("GO:0006575","cellular modified amino acid metabolic process",1.1845763067288293,4.525634723010098,0.8978741904401215,0.10312018,"protein nitrosylation"),
                     c("GO:0006749","glutathione metabolic process",0.33273992249082307,4.391485015365568,0.8261414519434429,0.18598634,"protein nitrosylation"),
                     c("GO:0006750","glutathione biosynthetic process",0.05906996283270419,1.7481965418254763,0.8256312820902051,0.6671716,"protein nitrosylation"),
                     c("GO:0018119","peptidyl-cysteine S-nitrosylation",9.119538616415149E-05,10.182564354913458,0.9171536206053698,0.38938339,"protein nitrosylation"),
                     c("GO:0018198","peptidyl-cysteine modification",0.01781021244384213,9.749062189632719,0.8982302901308479,0.28022578,"protein nitrosylation"),
                     c("GO:0019184","nonribosomal peptide biosynthetic process",0.19618345830058925,1.7481965418254763,0.8518380518151302,0.46794114,"protein nitrosylation"),
                     c("GO:0032881","regulation of polysaccharide metabolic process",0.057874563851903815,1.609498729889624,0.9885306505560326,0.09280684,"regulation of polysaccharide metabolic process"),
                     c("GO:0046686","response to cadmium ion",0.02872408190154112,3.4274075163814257,0.8964304710529556,-0,"response to peptide hormone"),
                     c("GO:0009704","de-etiolation",0.0020876349211090897,1.4473261028550646,0.9456012876384249,0.15144465,"response to peptide hormone"),
                     c("GO:0016036","cellular response to phosphate starvation",0.07367847780013462,2.2279450313908886,0.9284847842291067,0.18656216,"response to peptide hormone"),
                     c("GO:0043434","response to peptide hormone",0.19767216136121488,3.235047897478214,0.8414678144607108,0.40024394,"response to peptide hormone"),
                     c("GO:0080148","regulation of response to water deprivation",0.0001774612919951056,1.8738673377965496,0.969764734064127,-0,"regulation of response to water deprivation"),
                     c("GO:2000070","regulation of response to water deprivation",0.0013112417686305027,1.5082656543980195,0.9691897413335445,0.42513778,"regulation of response to water deprivation"),
                     c("GO:0120009","intermembrane lipid transfer",0.11033409355779031,1.609498729889624,0.8742450164595049,-0,"intermembrane lipid transfer"),
                     c("GO:0015713","phosphoglycerate transmembrane transport",0.003985484849390081,1.3112354007929257,0.8620085246994375,0.55136091,"intermembrane lipid transfer"),
                     c("GO:0015714","phosphoenolpyruvate transport",0.0052844029171875894,1.5623773116672997,0.8649661224630416,0.39685991,"intermembrane lipid transfer"),
                     c("GO:0015748","organophosphate ester transport",0.37259230624455725,1.3381790375665124,0.8900417202448845,0.37124925,"intermembrane lipid transfer"),
                     c("GO:0015790","UDP-xylose transmembrane transport",0.007542104909791988,1.609498729889624,0.8756101431121445,0.28243313,"intermembrane lipid transfer"),
                     c("GO:0035627","ceramide transport",0.019542924780961007,1.609498729889624,0.8863723408377396,0.60157124,"intermembrane lipid transfer"),
                     c("GO:0042873","aldonate transmembrane transport",0.023853755332342113,1.3112354007929257,0.851186606389152,0.60559877,"intermembrane lipid transfer"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)
}

activated.common.lrm.30.drought <- intersect(activated.genes.lrm.h2o.30,activated.genes.h2o.30.180)
length(activated.common.lrm.30.drought)

activated.common.lrm.30.drought.enrich.go <- enrichGO(gene = activated.common.lrm.30.drought, 
                                                      OrgDb = org.Taestivum.eg.db,
                                                      universe = rownames(gene.expression.ta), 
                                                      ont           = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.05,
                                                      readable      = FALSE,
                                                      keyType = "GID")


activated.common.lrm.30.drought.enrich.go.df <- as.data.frame(activated.common.lrm.30.drought.enrich.go)

write.table(x = activated.common.lrm.30.drought.enrich.go.df,file="activated_common_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


activated.specific.drought.lrm.30 <- setdiff(activated.genes.h2o.30.180, activated.genes.lrm.h2o.30)
length(activated.specific.drought.lrm.30)

activated.specific.drought.lrm.30.enrich.go <- enrichGO(gene = activated.specific.drought.lrm.30, 
                                                OrgDb = org.Taestivum.eg.db,
                                                universe = rownames(gene.expression.ta), 
                                                ont           = "BP",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                readable      = FALSE,
                                                keyType = "GID")


activated.specific.drought.lrm.30.enrich.go.df <- as.data.frame(activated.specific.drought.lrm.30.enrich.go)

write.table(x = activated.specific.drought.lrm.30.enrich.go,file="activated_specific_drought_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)

# Functional enrichment of the comparison between LRM and drought repressed
repressed.specific.lrm.30 <- setdiff(repressed.genes.lrm.h2o.30,repressed.genes.h2o.30.180)
length(repressed.specific.lrm.30)

repressed.specific.lrm.30.enrich.go <- enrichGO(gene = repressed.specific.lrm.30, 
                                                OrgDb = org.Taestivum.eg.db,
                                                universe = rownames(gene.expression.ta), 
                                                ont           = "BP",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                readable      = FALSE,
                                                keyType = "GID")


repressed.specific.lrm.30.enrich.go.df <- as.data.frame(repressed.specific.lrm.30.enrich.go)
head(repressed.specific.lrm.30.enrich.go.df)

write.table(x = repressed.specific.lrm.30.enrich.go.df,file="repressed_specific_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


repressed.common.lrm.30.drought <- intersect(repressed.genes.lrm.h2o.30,repressed.genes.h2o.30.180)
length(repressed.common.lrm.30.drought)

repressed.common.lrm.30.drought.enrich.go <- enrichGO(gene = repressed.common.lrm.30.drought, 
                                                      OrgDb = org.Taestivum.eg.db,
                                                      universe = rownames(gene.expression.ta), 
                                                      ont           = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.05,
                                                      readable      = FALSE,
                                                      keyType = "GID")


repressed.common.lrm.30.drought.enrich.go.df <- as.data.frame(repressed.common.lrm.30.drought.enrich.go)
head(repressed.common.lrm.30.drought.enrich.go.df)

write.table(x = repressed.common.lrm.30.drought.enrich.go.df,file="repressed_common_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


repressed.specific.drought.lrm.30 <- setdiff(repressed.genes.h2o.30.180, repressed.genes.lrm.h2o.30)
length(repressed.specific.drought.lrm.30)

repressed.specific.drought.lrm.30.enrich.go <- enrichGO(gene = repressed.specific.drought.lrm.30, 
                                                        OrgDb = org.Taestivum.eg.db,
                                                        universe = rownames(gene.expression.ta), 
                                                        ont           = "BP",
                                                        pAdjustMethod = "BH",
                                                        pvalueCutoff  = 0.05,
                                                        readable      = FALSE,
                                                        keyType = "GID")


repressed.specific.drought.lrm.30.enrich.go.df <- as.data.frame(repressed.specific.drought.lrm.30.enrich.go)
head(repressed.specific.drought.lrm.30.enrich.go.df)

write.table(x = repressed.specific.drought.lrm.30.enrich.go,file="repressed_specific_drought_lrm_30_enrich_go.tsv",quote = F,sep = "\t",row.names = F)



# Barplots for genes in all conditions

# Create barplots in FPKM
gene.expression.fpkm <- read.table(file="../TORRID/taestivium_gene_expression.tsv",header = T,sep = "\t", dec = ",", as.is = T, row.names = 1)
head(gene.expression.fpkm)
write.table(gene.expression.fpkm, "taestivium_gene_expression.tsv", quote = F, row.names = T, col.names = T, sep = "\t")

library(ggplot2)

barplots.ggplot2 <- function(gene, dataset = gene.expression.fpkm)
{
  conds <- c("h2o_180","LRM_180", "h2o_30", "LRM_30")
  
  means <- vector(mode = "numeric",length = length(conds))
  sds <- vector(mode = "numeric",length = length(conds))
  
  gene.matrix <- matrix(NA,nrow=length(conds),ncol=3)
  rownames(gene.matrix) <- conds
  for(i in 1:length(conds))
  {
    gene.matrix[i,] <- unlist(dataset[gene,paste(conds[i],1:3,sep="_")])
  }
  
  means <- rowMeans(gene.matrix)
  sds <- apply(X = gene.matrix,MARGIN = 1,FUN = sd)
  
  data.df <- data.frame(irrigation = rep(c("180", "30"), each = 2), 
                        treatment = rep(c("H20", "LRM"), 2),
                        condition = toupper(conds),
                        means_dat = means, sds_dat = sds)
  
  data.df$irrigation <- as.factor(data.df$irrigation)
  data.df$treatment <- as.factor(data.df$treatment)
  data.df$condition <- factor(data.df$condition, levels = toupper(conds))
  
  my_palette <- c("lightsalmon", "lightgreen", "chocolate3", "green4")
  
  p <- ggplot(data.df, aes(x=irrigation, y=means_dat, fill=condition)) + 
    geom_bar(position="dodge", stat="identity", colour="black", width = .6) + theme_bw() +
    ggtitle(gene) +
    geom_errorbar(aes(ymin=means_dat-sds_dat, ymax=means_dat+sds), width=.2,
                  position=position_dodge(.6)) + labs(x="Irrigation (ml)", y = "FPKM")+
    scale_fill_manual(values=my_palette) + 
    theme(plot.title = element_text(size=18, hjust = 0.5),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size=14, colour="black"),
          axis.text.x = element_text(size = 12, colour="black"),
          axis.title.y = element_text(size = 15, colour="black"))
    
  png(paste0(gene, ".png"), width = 1800, height = 1800, res = 300)
  plot(p)
  dev.off()
}


barplots.ggplot2(gene = "Traes_1AL_00744EBC6", dataset = gene.expression.fpkm)


# Barplots without grouping
barplots.ggplot2.no.group <- function(gene, dataset = gene.expression.fpkm)
{
  conds <- c("h2o_180","LRM_180", "h2o_30", "LRM_30")
  
  means <- vector(mode = "numeric",length = length(conds))
  sds <- vector(mode = "numeric",length = length(conds))
  
  gene.matrix <- matrix(NA,nrow=length(conds),ncol=3)
  rownames(gene.matrix) <- conds
  for(i in 1:length(conds))
  {
    gene.matrix[i,] <- unlist(dataset[gene,paste(conds[i],1:3,sep="_")])
  }
  
  means <- rowMeans(gene.matrix)
  sds <- apply(X = gene.matrix,MARGIN = 1,FUN = sd)
  
  data.df <- data.frame(irrigation = rep(c("180", "30"), each = 2), 
                        treatment = rep(c("H20", "LRM"), 2),
                        condition = toupper(conds),
                        means_dat = means, sds_dat = sds)
  
  data.df$irrigation <- as.factor(data.df$irrigation)
  data.df$treatment <- as.factor(data.df$treatment)
  data.df$condition <- factor(data.df$condition, levels = toupper(conds))
  
  my_palette <- c("lightsalmon", "lightgreen", "chocolate3", "green4")
  
  p <- ggplot(data.df, aes(x=condition, y=means_dat, fill=condition)) + 
    geom_bar(position="dodge", stat="identity", colour="black") + theme_bw() +
    ggtitle(gene) + scale_fill_manual(values=my_palette) +
    geom_errorbar(aes(ymin=means_dat-sds_dat, ymax=means_dat+sds), width=.2,
                  position=position_dodge(.6)) + labs(x="Irrigation (ml)", y = "FPKM") +
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size=14, colour="black"),
        axis.text.x = element_text(size = 12, colour="black"),
        axis.title.y = element_text(size = 15, colour="black"))
    
  
  
  png(paste0(gene, ".png"), width = 1800, height = 1800, res = 300)
  plot(p)
  dev.off()
}

barplots.ggplot2.no.group(gene = "Traes_1AL_00744EBC6", dataset = gene.expression.fpkm)

# Barplots only h2o

barplots.ggplot2.h2o <- function(gene, dataset = gene.expression.fpkm)
{
  conds <- c("h2o_180", "h2o_30")
  
  means <- vector(mode = "numeric",length = length(conds))
  sds <- vector(mode = "numeric",length = length(conds))
  
  gene.matrix <- matrix(NA,nrow=length(conds),ncol=3)
  rownames(gene.matrix) <- conds
  for(i in 1:length(conds))
  {
    gene.matrix[i,] <- unlist(dataset[gene,paste(conds[i],1:3,sep="_")])
  }
  
  means <- rowMeans(gene.matrix)
  sds <- apply(X = gene.matrix,MARGIN = 1,FUN = sd)
  
  data.df <- data.frame(irrigation = c("180", "30"), 
                        condition = toupper(conds),
                        means_dat = means, sds_dat = sds)
  
  data.df$irrigation <- as.factor(data.df$irrigation)
  data.df$condition <- factor(data.df$condition, levels = toupper(conds))
  
  my_palette <- c("lightsalmon", "chocolate3")
  
  p <- ggplot(data.df, aes(x=irrigation, y=means_dat, fill=condition)) + 
    geom_bar(position="dodge", stat="identity", colour="black", width = .6) + theme_bw() +
    ggtitle(gene) +
    geom_errorbar(aes(ymin=means_dat-sds_dat, ymax=means_dat+sds), width=.2,
                  position=position_dodge(.6)) + labs(x="Irrigation (ml)", y = "FPKM")+
    scale_fill_manual(values=my_palette) + 
    theme(plot.title = element_text(size=18, hjust = 0.5),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size=14, colour="black"),
          axis.text.x = element_text(size = 12, colour="black"),
          axis.title.y = element_text(size = 15, colour="black"))
  
  png(paste0(gene, ".png"), width = 1800, height = 1800, res = 300)
  plot(p)
  dev.off()
}

barplots.ggplot2.h2o(gene = "Traes_1AL_00744EBC6", dataset = gene.expression.fpkm)


# Analysis for specific genes
{
library(DESeq2)

gene.expression.ta <- read.csv(file = "../TORRID/gene_count_matrix.csv")
head(gene.expression.ta)
rownames(gene.expression.ta) <- gene.expression.ta$gene_id
gene.expression.ta <- gene.expression.ta[,2:ncol(gene.expression.ta)]
head(gene.expression.ta)
nrow(gene.expression.ta)
gene.expression.ta <- gene.expression.ta[-which(apply(X = gene.expression.ta,MARGIN = 1,FUN = sum) == 0),]
nrow(gene.expression.ta) # 75198 of 99386 genes showed expression in at least one sample
rownames(gene.expression.ta) <- unlist(strsplit(x = rownames(gene.expression.ta),split=".v2.2"))
experimental.design <- read.csv(file="../TORRID/experimental_design.csv")

colnames(gene.expression.ta) <- paste(experimental.design$treatment, 1:3,sep="_")
head(gene.expression.ta)


dds.wheat <- DESeqDataSetFromMatrix(countData=gene.expression.ta, colData=experimental.design, design = ~ treatment)
dds.wheat <- DESeq(dds.wheat)
}

# Create contrast (30 vs 180)
res.h2o.30.180 <- results(dds.wheat,contrast=c("treatment","H2O_30","H2O_180"))
res.h2o.30.180$padj[is.na(res.h2o.30.180$padj)] <- 1

genes.h2o.30.180 <- rownames(res.h2o.30.180)

# Create contrast (LRM 30 vs H20 30)
res.lrm.30 <- results(dds.wheat,contrast=c("treatment","LRM_30","H2O_30"))
res.lrm.30$padj[is.na(res.lrm.30$padj)] <- 1

genes.lrm.30 <- rownames(res.lrm.30)

# Analysis of OST genes
ost <- read.table("result_2/ost_genes.csv", sep = "\t", header = T)
ost

sapply(ost$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))


ost.h2o.30.180 <- res.h2o.30.180[ost$Gene,]
as.data.frame(ost.h2o.30.180)
subset(ost.h2o.30.180, padj < 0.05) # No significant exp for OST genes

ost.res.lrm.30 <- res.lrm.30[ost$Gene,]
as.data.frame(ost.res.lrm.30)
subset(ost.res.lrm.30, padj < 0.05) # Significant for OST genes (up)
sig.genes.ost.lrm.30 <- as.data.frame(subset(ost.res.lrm.30, padj < 0.05))


# SCAB
scab <- read.table("result_2/scab_genes.csv", sep = "\t", header = T)
scab
sapply(scab$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

scab_genes_in_dataset <- scab$Gene[scab$Gene %in% rownames(res.h2o.30.180)]
scab.h2o.30.180 <- res.h2o.30.180[scab_genes_in_dataset,]
as.data.frame(scab.h2o.30.180)
subset(scab.h2o.30.180, padj < 0.05) # No significant exp for SCAB genes


scab.res.lrm.30 <- res.lrm.30[scab_genes_in_dataset,]
as.data.frame(scab.res.lrm.30)
subset(scab.res.lrm.30, padj < 0.05) # No significant exp for SCAB genes

# VPS34
vps34 <- read.table("result_2/vps34.csv", sep = "\t", header = T)
vps34
sapply(vps34$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

vps34_genes_in_dataset <- vps34$Gene[vps34$Gene %in% rownames(res.h2o.30.180)]
vps34.h2o.30.180 <- res.h2o.30.180[vps34_genes_in_dataset,]
as.data.frame(vps34.h2o.30.180)
subset(vps34.h2o.30.180, padj < 0.05) # No significant exp for vps34 genes

vps34.res.lrm.30 <- res.lrm.30[vps34_genes_in_dataset,]
as.data.frame(vps34.res.lrm.30)
subset(vps34.res.lrm.30, padj < 0.05) # No significant exp for VPS34 genes

# PEPCK
pepck <- read.table("result_2/pepck_genes.csv", sep = "\t", header = T)
pepck
sapply(pepck$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

pepck_genes_in_dataset <- pepck$Gene[pepck$Gene %in% rownames(res.h2o.30.180)]
pepck.h2o.30.180 <- res.h2o.30.180[pepck_genes_in_dataset,]
as.data.frame(pepck.h2o.30.180)
subset(pepck.h2o.30.180, padj < 0.05) # Significant exp for pepck genes (down)

pepck.res.lrm.30 <- res.lrm.30[pepck_genes_in_dataset,]
as.data.frame(pepck.res.lrm.30)
subset(pepck.res.lrm.30, padj < 0.05) # Significant exp for pepck genes (down)


# ME
me <- read.table("result_2/me_genes.csv", sep = "\t", header = T)

me_genes_in_dataset <- me$Gene[me$Gene %in% rownames(res.h2o.30.180)]
me.h2o.30.180 <- res.h2o.30.180[me_genes_in_dataset,]
as.data.frame(me.h2o.30.180)
subset(me.h2o.30.180, padj < 0.05) # Significant exp for me genes (down)

me.res.lrm.30 <- res.lrm.30[me_genes_in_dataset,]
as.data.frame(me.res.lrm.30)
subset(me.res.lrm.30, padj < 0.05)  # significant for NADP-ME2 (down)

# pp2c
pp2c <- read.table("result_2/pp2c_genes.csv", sep = "\t", header = T)
sapply(pp2c$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

pp2c_genes_in_dataset <- pp2c$Gene[pp2c$Gene %in% rownames(res.h2o.30.180)]
pp2c.h2o.30.180 <- res.h2o.30.180[pp2c_genes_in_dataset,]
as.data.frame(pp2c.h2o.30.180)
subset(pp2c.h2o.30.180, padj < 0.05) # Significant exp for pp2c genes (down)

pp2c.res.lrm.30 <- res.lrm.30[pp2c_genes_in_dataset,]
as.data.frame(pp2c.res.lrm.30)
sig.genes.pp2c.lrm.30 <- as.data.frame(subset(pp2c.res.lrm.30, padj < 0.05))  # significant for some components (down)

# PYR/PYL/RCAR
pyr <- read.table("result_2/pyr_pyl_rcar_genes.csv", sep = "\t", header = T)
pyr
sapply(pyr$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

pyr_genes_in_dataset <- pyr$Gene[pyr$Gene %in% rownames(res.h2o.30.180)]
pyr.h2o.30.180 <- res.h2o.30.180[pyr_genes_in_dataset,]
as.data.frame(pyr.h2o.30.180)
subset(pyr.h2o.30.180, padj < 0.05) # PYL4 (down)


pyr.res.lrm.30 <- res.lrm.30[pyr_genes_in_dataset,]
as.data.frame(pyr.res.lrm.30)
sig.genes.pyl.lrm.30 <- as.data.frame(subset(pyr.res.lrm.30, padj < 0.05)) # PYL7 related (up)


# ABA biosynthesis
aba <- read.table("result_2/aba.csv", sep = "\t", header = T)
aba
sapply(aba$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

aba_genes_in_dataset <- aba$Gene[aba$Gene %in% rownames(res.h2o.30.180)]
aba.h2o.30.180 <- res.h2o.30.180[aba_genes_in_dataset,]
as.data.frame(aba.h2o.30.180)
subset(aba.h2o.30.180, padj < 0.05) # No significant changes for ABA synthesis genes


aba.res.lrm.30 <- res.lrm.30[aba_genes_in_dataset,]
as.data.frame(aba.res.lrm.30)
sig.genes.aba.lrm.30 <- as.data.frame(subset(aba.res.lrm.30, padj < 0.05)) # BETA-CAROTENE DIOXYGENASE (up)

# MDH-NADP
mdh <- read.table("result_2/mdh.csv", sep = "\t", header = T)
mdh
sapply(mdh$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

mdh_genes_in_dataset <- mdh$Gene[mdh$Gene %in% rownames(res.h2o.30.180)]
mdh.h2o.30.180 <- res.h2o.30.180[mdh_genes_in_dataset,]
as.data.frame(mdh.h2o.30.180)
subset(mdh.h2o.30.180, padj < 0.05) # No significant changes for mdh genes


mdh.res.lrm.30 <- res.lrm.30[mdh_genes_in_dataset,]
as.data.frame(mdh.res.lrm.30)
sig.genes.mdh.lrm.30 <- as.data.frame(subset(mdh.res.lrm.30, padj < 0.05)) # No significant changes for mdh synthesis genes

# CA
ca <- read.table("result_2/ca_genes.csv", sep = "\t", header = T)
ca
sapply(ca$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

ca_genes_in_dataset <- ca$Gene[ca$Gene %in% rownames(res.h2o.30.180)]
ca.h2o.30.180 <- res.h2o.30.180[ca_genes_in_dataset,]
as.data.frame(ca.h2o.30.180)
subset(ca.h2o.30.180, padj < 0.05) # No significant changes for ca genes


ca.res.lrm.30 <- res.lrm.30[ca_genes_in_dataset,]
as.data.frame(ca.res.lrm.30)
sig.genes.ca.lrm.30 <- as.data.frame(subset(ca.res.lrm.30, padj < 0.05)) # CA1 (up)

# slac1
slac1 <- read.table("result_2/slac1_genes.csv", sep = "\t", header = T)
slac1
sapply(slac1$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

slac1_genes_in_dataset <- slac1$Gene[slac1$Gene %in% rownames(res.h2o.30.180)]
slac1.h2o.30.180 <- res.h2o.30.180[slac1_genes_in_dataset,]
as.data.frame(slac1.h2o.30.180)
subset(slac1.h2o.30.180, padj < 0.05) # No significant changes for slac1 genes


slac1.res.lrm.30 <- res.lrm.30[slac1_genes_in_dataset,]
as.data.frame(slac1.res.lrm.30)
sig.genes.slac1.lrm.30 <- as.data.frame(subset(slac1.res.lrm.30, padj < 0.05)) # No signifislacant changes for slac1 genes

# quac1
quac1 <- read.table("result_2/quac1_genes.csv", sep = "\t", header = T)
quac1
sapply(quac1$Gene, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

quac1_genes_in_dataset <- quac1$Gene[quac1$Gene %in% rownames(res.h2o.30.180)]
quac1.h2o.30.180 <- res.h2o.30.180[quac1_genes_in_dataset,]
as.data.frame(quac1.h2o.30.180)
subset(quac1.h2o.30.180, padj < 0.05) # No significant changes for quac1 genes


quac1.res.lrm.30 <- res.lrm.30[quac1_genes_in_dataset,]
as.data.frame(quac1.res.lrm.30)
sig.genes.quac1.lrm.30 <- as.data.frame(subset(quac1.res.lrm.30, padj < 0.05)) # ALMT (down)

# Analysis of Phosphate Starvation Response
psr <- read.table("enrichment_tables/h2o_30_180/activated_genes_h2o_30_180_enrich_go.tsv", sep="\t", header = T)
psr.genes <- strsplit(psr[1,"geneID"], split = "/")[[1]]
psr_genes_in_dataset <- psr.genes[psr.genes %in% rownames(res.h2o.30.180)]
sapply(psr_genes_in_dataset, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))
sapply(psr_genes_in_dataset, function(x) barplots.ggplot2.h2o(gene = x, dataset = gene.expression.fpkm))


psr.h2o.30.180 <- res.h2o.30.180[psr_genes_in_dataset,]
as.data.frame(psr.h2o.30.180)
subset(psr.h2o.30.180, padj < 0.05)

psr.res.lrm.30 <- res.lrm.30[psr_genes_in_dataset,]
as.data.frame(psr.res.lrm.30)
sig.genes.psr.lrm.30 <- as.data.frame(subset(psr.res.lrm.30, padj < 0.05))

# 
ht1 <- read.table("result_2/ht1_genes.csv", sep="\t", header = T)
ht1_genes_in_dataset <- ht1[ht1 %in% rownames(res.h2o.30.180)]
sapply(ht1_genes_in_dataset, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))
sapply(ht1_genes_in_dataset, function(x) barplots.ggplot2.h2o(gene = x, dataset = gene.expression.fpkm))


ht1.h2o.30.180 <- res.h2o.30.180[ht1_genes_in_dataset,]
as.data.frame(ht1.h2o.30.180)
subset(ht1.h2o.30.180, padj < 0.05)

ht1.res.lrm.30 <- res.lrm.30[ht1_genes_in_dataset,]
as.data.frame(ht1.res.lrm.30)
sig.genes.ht1.lrm.30 <- as.data.frame(subset(ht1.res.lrm.30, padj < 0.05))


# Heatmaps for Result 2
barplot.heatmap.gene <- function(gene.id,gene.expression,max.fc)
{
  expression.h2o.180 <- as.numeric(gene.expression[gene.id, paste("h2o_180", 1:3,sep="_")])
  expression.h2o.30 <- as.numeric(gene.expression[gene.id, paste("h2o_30", 1:3,sep="_")])
  expression.lrm.30 <- as.numeric(gene.expression[gene.id, paste("LRM_30", 1:3,sep="_")])
  
  means <- c(mean(expression.h2o.180),
             mean(expression.h2o.30),
             mean(expression.lrm.30))
  
  sds <- c(sd(expression.h2o.180),
           sd(expression.h2o.30),
           sd(expression.lrm.30))
  
  # color palette blue to red
  rwb <- colorRampPalette(colors = c("blue", "white", "red"))
  range.colors <- rwb(max.fc*200)
  
  log2.fcs <- 100*log2(c(means[2]/means[1],
                         means[3]/means[1]))
  
  log2.fcs <- log2.fcs + max.fc*100
  log2.fcs <- round(log2.fcs)
  log2.fcs[log2.fcs > max.fc*200] <- max.fc*200
  log2.fcs[log2.fcs < 1] <- 1
  
  png(filename = paste0("heatmap_", gene.id,".png"),
      width = 700,height=450)
  
  plot(1, 1, col = "white", axes=F,xlab="",ylab="",xlim=c(-1,21),ylim=c(-1,11))  
  polygon(x = c(0, 0, 10, 10),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][1],lwd=6)
  polygon(x = c(10, 10, 20, 20),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][2],lwd=6)
  dev.off()
}

barplot.heatmap.gene(gene.id = "Traes_6AL_91F3E548F", 
                     gene.expression = gene.expression.fpkm, max.fc = 8)

# Heatmaps DEGs OST genes
sapply(rownames(sig.genes.ost.lrm.30[-1,]), function(x) barplot.heatmap.gene(gene.id = x, 
                       gene.expression = gene.expression.fpkm, max.fc = 2))

# Heatmaps DEGs PP2C genes
sapply(rownames(sig.genes.pp2c.lrm.30), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                             gene.expression = gene.expression.fpkm, max.fc = 2))

# Heatmaps DEGs PYL genes
sapply(rownames(sig.genes.pyl.lrm.30), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                         gene.expression = gene.expression.fpkm, max.fc = 2))
# Heatmaps DEGs ABA genes
sapply(rownames(sig.genes.aba.lrm.30[-2,]), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                        gene.expression = gene.expression.fpkm, max.fc = 2))

# Heatmaps DEGs PSR genes
sapply(rownames(sig.genes.psr.lrm.30), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                        gene.expression = gene.expression.fpkm, max.fc = 2))

# Heatmaps DEGs CA genes
sapply(rownames(sig.genes.ca.lrm.30), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                        gene.expression = gene.expression.fpkm, max.fc = 2))

# Heatmaps DEGs QUAC1 genes
sapply(rownames(sig.genes.quac1.lrm.30), function(x) barplot.heatmap.gene(gene.id = x, 
                                                                          gene.expression = gene.expression.fpkm, max.fc = 2))

# Result 3
## Differential gene expression analysis of LRM treated plants under Full Irrigation
res.lrm.180 <- results(dds.wheat,contrast=c("treatment","LRM_180","H2O_180"))
res.lrm.180$padj[is.na(res.lrm.180$padj)] <- 1

genes.lrm.180 <- rownames(res.lrm.180)

activated.genes.lrm.180 <- genes.lrm.180[res.lrm.180$log2FoldChange > log2(1.5) & res.lrm.180$padj < 0.05]
length(activated.genes.lrm.180)

repressed.genes.lrm.180 <- genes.lrm.180[res.lrm.180$log2FoldChange < -log2(1.5) & res.lrm.180$padj < 0.05]
length(repressed.genes.lrm.180)

# Volcano plot
logfc.lrm.180 <- res.lrm.180$log2FoldChange
names(logfc.lrm.180) <- genes.lrm.180
logqval.lrm.180 <- -log10(res.lrm.180$padj)
names(logqval.lrm.180) <- genes.lrm.180

plot(logfc.lrm.180[logqval.lrm.180 > 0],logqval.lrm.180[logqval.lrm.180 > 0],pch=19,cex=0.7,xlim=c(-7,7),ylim=c(0,5),col="grey",xlab="log2FC",ylab="-log10(q-value)")
points(logfc.lrm.180[activated.genes.lrm.180],logqval.lrm.180[activated.genes.lrm.180],cex=0.7,col="red",pch=19)
points(logfc.lrm.180[repressed.genes.lrm.180],logqval.lrm.180[repressed.genes.lrm.180],cex=0.7,col="blue",pch=19)

# With ggplot2
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

df.lrm.180 <- as.data.frame(res.lrm.180)
head(df.lrm.180)

df.lrm.180$diffexpressed <- "NO"
df.lrm.180$diffexpressed[df.lrm.180$log2FoldChange > log2(1.5) & df.lrm.180$padj < 0.05] <- "UP"
df.lrm.180$diffexpressed[df.lrm.180$log2FoldChange < -log2(1.5) & df.lrm.180$padj < 0.05] <- "DOWN"
df.lrm.180$names <- rownames(df.lrm.180)

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

volcano_plot_lrm.180 <-  ggplot(data = df.lrm.180, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "gray", linetype = 'dashed') + xlim(c(-8,8)) + ylim(c(0,5)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />

png("volcano_lrm_180.png", width = 3000, height = 2000, res = 300)
plot(volcano_plot_lrm.180)
dev.off()

# Functional enrichment of activated and repressed genes by LRM under irrigation
library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)
library(org.Taestivum.eg.db)

activated.genes.lrm.180.enrich.go <- enrichGO(gene = activated.genes.lrm.180, 
                                              OrgDb = org.Taestivum.eg.db,
                                              universe = rownames(gene.expression.ta), 
                                              ont           = "BP",
                                              pAdjustMethod = "BH",
                                              pvalueCutoff  = 0.05,
                                              readable      = FALSE,
                                              keyType = "GID")

activated.genes.lrm.180.enrich.go.df <- as.data.frame(activated.genes.lrm.180.enrich.go)

write.table(x = activated.genes.lrm.180.enrich.go.df,file="activated_genes_lrm_180_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


repressed.genes.lrm.180.enrich.go <- enrichGO(gene = repressed.genes.lrm.180, 
                                              OrgDb = org.Taestivum.eg.db,
                                              universe = rownames(gene.expression.ta), 
                                              ont           = "BP",
                                              pAdjustMethod = "BH",
                                              pvalueCutoff  = 0.05,
                                              readable      = FALSE,
                                              keyType = "GID")

repressed.genes.lrm.180.enrich.go.df <- as.data.frame(repressed.genes.lrm.180.enrich.go)

write.table(x = repressed.genes.lrm.180.enrich.go.df,file="repressed_genes_lrm_180_enrich_go.tsv",quote = F,sep = "\t",row.names = F)


length(intersect(activated.genes.lrm.h2o.30,activated.genes.lrm.180))

length(intersect(activated.genes.lrm.h2o.30,repressed.genes.lrm.180))
length(intersect(repressed.genes.lrm.h2o.30,repressed.genes.lrm.180))


length(intersect(activated.genes.lrm.180,activated.genes.h2o.30.180))
length(intersect(activated.genes.lrm.180,repressed.genes.h2o.30.180))
length(intersect(repressed.genes.lrm.180,repressed.genes.h2o.30.180))

# Barplots of activated and repressed genes
sapply(activated.genes.lrm.180, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))
sapply(repressed.genes.lrm.180, function(x) barplots.ggplot2(gene = x, dataset = gene.expression.fpkm))

# Treemap activated genes
{
  library(treemap)
  revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
  revigo.data <- rbind(c("GO:0009641","shade avoidance",0.0003204162216578296,0.772475649157364,-0,"shade avoidance"),
                       c("GO:0048571","long-day photoperiodism",0.0008059700344777714,0.7338628319070892,0.48734826,"shade avoidance"),
                       c("GO:0048574","long-day photoperiodism, flowering",0.0005964671203168828,0.7351363827379533,0.6593882,"shade avoidance"),
                       c("GO:0042752","regulation of circadian rhythm",0.08089030752760239,0.8729214057225224,0,"regulation of circadian rhythm"),
                       c("GO:0042754","negative regulation of circadian rhythm",0.004387237496545666,0.8863096948272162,0.11197257,"regulation of circadian rhythm"),
                       c("GO:0044092","negative regulation of molecular function",0.0036749275883986455,0.7382379715745363,0.10497025,"regulation of circadian rhythm"),
                       c("GO:1902446","regulation of shade avoidance",4.6830063165375095E-05,0.8533132663605822,0.29413431,"regulation of circadian rhythm"),
                       c("GO:2000030","regulation of response to red or far red light",0.004793919624034451,0.8358687121716452,0.11255422,"regulation of circadian rhythm"));
  
  stuff <- data.frame(revigo.data);
  names(stuff) <- revigo.names;
  
  stuff$frequency <- as.numeric( as.character(stuff$frequency) );
  stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
  stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
  
  png("treemap_act_lrm_180.png", width = 3000, height = 2000, res = 300)
  
  treemap(
    stuff,
    index = c("representative","description"),
    vSize = "uniqueness",
    type = "categorical",
    vColor = "representative",
    title = " ",
    inflate.labels = TRUE,     
    lowerbound.cex.labels = 0,   
    position.legend = "none"
  )
  
  dev.off()
}

# Treemap repressed genes
{
  library(treemap)
  revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
  revigo.data <- rbind(c("GO:0006690","icosanoid metabolic process",0.15592932189969946,0.943056910277058,0.00762984,"icosanoid metabolic process"),
                       c("GO:0019372","lipoxygenase pathway",0.019244691220802565,0.9333520111053946,0.49026099,"icosanoid metabolic process"),
                       c("GO:0019464","glycine decarboxylation via glycine cleavage system",0.07523619358542499,0.9079863313281709,0.32747518,"icosanoid metabolic process"),
                       c("GO:0030497","fatty acid elongation",0.07301299795515297,0.9358591188826687,0.32675746,"icosanoid metabolic process"),
                       c("GO:2001293","malonyl-CoA metabolic process",0.059511151322525345,0.9355170874352527,0.69078938,"icosanoid metabolic process"),
                       c("GO:2001295","malonyl-CoA biosynthetic process",0.053519367977523935,0.9358630638304197,0.22065194,"icosanoid metabolic process"),
                       c("GO:0010265","SCF complex assembly",0.00911953861641515,0.9952623385930559,0,"SCF complex assembly"),
                       c("GO:0034198","cellular response to amino acid starvation",0.06431985938755938,0.9444087521047673,0.006329,"cellular response to amino acid starvation"),
                       c("GO:0010255","glucose mediated signaling pathway",0.0005668902383176985,0.7792404521771,0.68472684,"cellular response to amino acid starvation"),
                       c("GO:0031929","TOR signaling",0.05839955350738934,0.8244146422917595,0.28523874,"cellular response to amino acid starvation"),
                       c("GO:0045472","response to ether",0.004956592475029964,0.9405507726884759,0.16610387,"cellular response to amino acid starvation"),
                       c("GO:0071324","cellular response to disaccharide stimulus",0.0006359029629824618,0.918440165661677,0.68801059,"cellular response to amino acid starvation"),
                       c("GO:0071329","cellular response to sucrose stimulus",0.0006211145219828697,0.9118559757590926,0.40470041,"cellular response to amino acid starvation"),
                       c("GO:1901355","response to rapamycin",0.00038696420615599424,0.945894781134579,0.39614169,"cellular response to amino acid starvation"),
                       c("GO:0040019","positive regulation of embryonic development",0.003036559885249585,0.7947335347313873,0.09411582,"positive regulation of embryonic development"),
                       c("GO:0033674","positive regulation of kinase activity",0.0016538406517877205,0.7404739477066846,0.27088337,"positive regulation of embryonic development"),
                       c("GO:1903313","positive regulation of mRNA metabolic process",0.2376527116036122,0.6188147824302259,0.5010169,"positive regulation of embryonic development"),
                       c("GO:2000234","positive regulation of rRNA processing",0.006615362607150882,0.7182209416846864,0.34470745,"positive regulation of embryonic development"),
                       c("GO:0055062","phosphate ion homeostasis",0.04733779963969443,0.9514409837742751,-0,"phosphate ion homeostasis"),
                       c("GO:0007035","vacuolar acidification",0.07733615220736707,0.7808175178660604,0.56567559,"phosphate ion homeostasis"),
                       c("GO:0055081","monoatomic anion homeostasis",0.027077635470253197,0.9507812383595506,0.51004304,"phosphate ion homeostasis"),
                       c("GO:0060151","peroxisome localization",0.00026372719782605974,0.9237601749115794,-0,"peroxisome localization"),
                       c("GO:0006817","phosphate ion transport",0.2426906405041399,0.8723866387719781,0.25520912,"peroxisome localization"),
                       c("GO:0007034","vacuolar transport",0.4278024859763679,0.8569329935756456,0.40954151,"peroxisome localization"),
                       c("GO:0015698","inorganic anion transport",0.8722912628806088,0.897757874539155,0.24257536,"peroxisome localization"),
                       c("GO:0032507","maintenance of protein location in cell",0.09837517426945347,0.8349959661860205,0.43055422,"peroxisome localization"),
                       c("GO:0034067","protein localization to Golgi apparatus",0.052126789783395674,0.858992929208611,0.19212086,"peroxisome localization"),
                       c("GO:0034486","vacuolar transmembrane transport",0.025445977479964865,0.9121999895003551,0.12320799,"peroxisome localization"),
                       c("GO:0051235","maintenance of location",0.2512506831027372,0.909139211462361,0.18246122,"peroxisome localization"),
                       c("GO:0051645","Golgi localization",0.012962068536142506,0.907799207770462,0.44322855,"peroxisome localization"),
                       c("GO:0051646","mitochondrion localization",0.04394385243028803,0.9014984337679232,0.57915923,"peroxisome localization"),
                       c("GO:0072665","protein localization to vacuole",0.21067613048018954,0.8458289400802537,0.66552598,"peroxisome localization"),
                       c("GO:0072666","establishment of protein localization to vacuole",0.20075801604979646,0.8382765837110473,0.59802436,"peroxisome localization"),
                       c("GO:0098656","monoatomic anion transmembrane transport",0.4548751919662879,0.8934159204464921,0.27102916,"peroxisome localization"),
                       c("GO:1905011","transmembrane phosphate ion transport from cytosol to vacuole",0.0009292070428077058,0.899367944334989,0.53391375,"peroxisome localization"),
                       c("GO:1902659","regulation of glucose mediated signaling pathway",0.004510474504875602,0.8556767355662085,-0,"regulation of glucose mediated signaling pathway"),
                       c("GO:0010565","regulation of cellular ketone metabolic process",0.06463288138871742,0.7346156201459705,0.22304529,"regulation of glucose mediated signaling pathway"),
                       c("GO:0010929","positive regulation of auxin mediated signaling pathway",0.00023415031582687548,0.7949223547611656,0.37948485,"regulation of glucose mediated signaling pathway"),
                       c("GO:0016241","regulation of macroautophagy",0.057734073662407695,0.6753851584128568,0.22128492,"regulation of glucose mediated signaling pathway"),
                       c("GO:0017148","negative regulation of translation",0.21779183534115995,0.6616381519770471,0.64688619,"regulation of glucose mediated signaling pathway"),
                       c("GO:0019216","regulation of lipid metabolic process",0.10650635207906255,0.7851157533281655,0.20726298,"regulation of glucose mediated signaling pathway"),
                       c("GO:0034249","negative regulation of amide metabolic process",0.23996217313971518,0.7053793003212597,0.68743154,"regulation of glucose mediated signaling pathway"),
                       c("GO:0046890","regulation of lipid biosynthetic process",0.06254771120777494,0.7264899303892461,0.23788923,"regulation of glucose mediated signaling pathway"),
                       c("GO:0050687","negative regulation of defense response to virus",0.0015798984467897596,0.8250227541409683,0.34242069,"regulation of glucose mediated signaling pathway"),
                       c("GO:0051248","negative regulation of protein metabolic process",0.3529088912741001,0.6973885750413564,0.31114933,"regulation of glucose mediated signaling pathway"),
                       c("GO:1900459","positive regulation of brassinosteroid mediated signaling pathway",0.0002021086936610925,0.7962258452701457,0.33965984,"regulation of glucose mediated signaling pathway"),
                       c("GO:1902661","positive regulation of glucose mediated signaling pathway",0.00035738732415680993,0.7910828546808045,0.38796548,"regulation of glucose mediated signaling pathway"),
                       c("GO:1903311","regulation of mRNA metabolic process",0.4779772015478174,0.7506124286150798,0.27758205,"regulation of glucose mediated signaling pathway"),
                       c("GO:2000232","regulation of rRNA processing",0.0076480887369557325,0.8074333137516976,0.10565847,"regulation of glucose mediated signaling pathway"));
  
  stuff <- data.frame(revigo.data);
  names(stuff) <- revigo.names;
  
  stuff$frequency <- as.numeric( as.character(stuff$frequency) );
  stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
  stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
  
  png("treemap_rep_lrm_180.png", width = 3000, height = 2000, res = 300)
  
  treemap(
    stuff,
    index = c("representative","description"),
    vSize = "uniqueness",
    type = "categorical",
    vColor = "representative",
    title = "Revigo TreeMap",
    inflate.labels = FALSE,      
    lowerbound.cex.labels = 0,  
    position.legend = "none"
  )
  
  dev.off()
  
  
}