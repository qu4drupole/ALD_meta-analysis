#######################################################
# PART 1
# AGGREGATION AND PROCESSING OF HUMAN ALD MICROARRAY DATA
# 
# The purpose of this script is to download and clean
# various microarray datasets for the purpose of downstream ALD modeling
#
#
#######################################################

library(limma)
library(readxl)
library(dplyr)
library(GEOquery)
library(sva)
library(gplots)
library(org.Hs.eg.db)


datapath <- "D:/Dropbox (SBG)/David-Smith/Jefferson-Desktop/Data/old_rat_mircoarray/"

# Reading in the old etoh rat data
metaf <- read_excel(paste0(datapath,"old_rat_meta.xls"))
names(metaf)[3] <- "trmt"
metaf$diet <- character(nrow(metaf))
metaf$diet[which(metaf$trmt=="chow")] <- "chow"
metaf$diet[which(metaf$trmt!="chow")] <- "liquid"
metaf$diet <- as.factor(metaf$diet)
metaf$trmt2 <- character(nrow(metaf))
metaf$trmt2[which(metaf$trmt %in% c("chow","CHO"))] <- "ctrl"
metaf$trmt2[which(metaf$trmt == "etoh")] <- "etoh"
metaf$trmt2 <- as.factor(metaf$trmt2)

load(paste0(datapath,"Old-Young_combined_llm.Rdata"))
llm.data <- c.old.young.llm

colnames(llm.data)[sapply(colnames(llm.data), function(x) grepl("Avg",x), USE.NAMES = F)] <- metaf$ID
llm.data <- data.frame(llm.data) %>% rename_at(vars(contains("ctrl")), ~paste0("y.",.)) %>%
  rename_if(!grepl("y.",names(.)), ~paste0("o.",.))
llm.data <- llm.data[!is.na(llm.data$Gene.symbol),]
rownames(llm.data) <- toupper(llm.data$Gene.symbol)
llm.data <- subset(llm.data, select = -c(Gene.symbol))

##################################
# Need to have some standarized codes for ALD
# `liverOutcome`:
# -Healthy
# -NAFLD (Steatosis)
# -NASH
# -ASH
# -Hep (Hepatitis)
# -Cir (Cirrhosis)
# -HCC
##################################

#human NASH
#GSE48452
#platform: Affymetrix Human Gene 1.1 ST Array 
#Human liver biopsy of different phases from control to NASH
gset <- getGEO("GSE48452", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6247", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'Healthy'
gset$liverOutcome[gset$`group:ch1`=='Nash'] <- 'NASH'
gset$liverOutcome[gset$`group:ch1`=='Steatosis'] <- 'NAFLD'

#Annotate and format
ex <- exprs(gset)
names(fData(gset))
ex <- data.frame(ex)
a <- sapply(fData(gset)$gene_assignment, function(x) trimws(unlist(strsplit(x,"//"))[2]), USE.NAMES = F)
ex$Gene.symbol <- a
ex <- ex[!is.na(ex$Gene.symbol),]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

GSE48452 <- list('ex'=ex, 'meta'=pData(gset))
GSE48452$meta$GSE <- 'GSE48452'

#Human HCC
#GSE62232
gset <- getGEO("GSE62232", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'trash'
gset$liverOutcome[gset$source_name_ch1=='Non tumor liver'] <- 'Healthy'
gset$liverOutcome[grepl('Al',gset$`etiology:ch1`)] <- 'HCC'
gset <- gset[gset$liverOutcome != 'trash',]

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
a <- sapply(fData(gset)$Gene.symbol, function(x) trimws(unlist(strsplit(x,"//"))[1]), USE.NAMES = F)
ex$Gene.symbol <- a
ex <- ex[!is.na(ex$Gene.symbol),]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))
ex <- ex[,gset$geo_accession[gset$liverOutcome != 'trash']]

GSE62232 <- list("ex"=ex, "meta"=pData(gset))
GSE62232$meta <- GSE62232$meta[GSE62232$meta$liverOutcome != 'trash',]
GSE62232$meta$GSE <- 'GSE62232'

#Human ALD
#GSE94417 (superseries)
# GSE94397
# GSE94399

#Human ALD
#GSE103580
gset <- getGEO("GSE103580", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13667", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'Cir'
gset$liverOutcome[grepl('Mild',gset$source_name_ch1)] <- 'Hep'
gset$liverOutcome[grepl('steato',gset$source_name_ch1)] <- 'NAFLD'

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
a <- sapply(fData(gset)$Gene.Symbol, function(x) trimws(unlist(strsplit(x,"///"))[1]), USE.NAMES = F)
ex$Gene.symbol <- a
ex <- ex[!is.na(ex$Gene.symbol),]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

GSE103580 <- list("ex"=ex, "meta"=pData(gset))
GSE103580$meta$GSE = 'GSE103580'

# Human NASH
# GSE83452
gset <- getGEO("GSE83452", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16686", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'trash'
gset$liverOutcome[!grepl('no',gset$`liver status:ch1`)] <- 'NASH'
gset$liverOutcome[gset$`liver status:ch1`=='undefined'] <- 'trash'
# gset <- gset[gset$liverOutcome != 'trash',]

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
# a <- sapply(fData(gset)$Gene.symbol, function(x) trimws(unlist(strsplit(x,"//"))[1]), USE.NAMES = F)
a <- fData(gset); a <- a$GB_ACC; 
b <- mapIds(org.Hs.eg.db, a, "SYMBOL", "ACCNUM")
ex$Gene.symbol <- unname(b)
ex <- ex[ex$Gene.symbol!="NULL",]
ex <- ex[ex$Gene.symbol!="NULL",] #weird...but have to do this again
ex$Gene.symbol <- unlist(ex$Gene.symbol)
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

keep <- rownames(pData(gset))[pData(gset)$liverOutcome != "trash"]
GSE83452 <- list("ex"=ex[,keep], "meta"=pData(gset)[keep,])
# GSE83452$meta <- GSE83452$meta[GSE83452$meta$liverOutcome != 'trash',]
GSE83452$meta$GSE <- 'GSE83452'

# Human NAFLD
# GSE49541
gset <- getGEO("GSE49541", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'NAFLD'
gset$liverOutcome[grepl('advanced',gset$`Stage:ch1`)] <- 'trash'
# gset$liverOutcome[grepl('steato',gset$source_name_ch1)] <- 'NAFLD'

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
a <- sapply(fData(gset)$Gene.Symbol, function(x) trimws(unlist(strsplit(x,"///"))[1]), USE.NAMES = F)
ex$Gene.symbol <- a
ex <- ex[!is.na(ex$Gene.symbol),]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

keep <- rownames(pData(gset))[pData(gset)$liverOutcome == "NAFLD"]
GSE49541 <- list("ex"=ex[,keep], "meta"=pData(gset)[keep,])
GSE49541$meta$GSE = 'GSE49541'

# Human NAFLD, NASH
# GSE33814
gset <- getGEO("GSE33814", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6884", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'NAFLD'
gset$liverOutcome[gset$`diagnosis:ch1`=="steatohepatitis"] <- 'NASH'
gset$liverOutcome[gset$`diagnosis:ch1`=="normal"] <- 'Healthy'

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
a <- fData(gset); ex$Gene.symbol <- a$Gene.symbol
ex <- ex[ex$Gene.symbol != "",]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

# keep <- rownames(pData(gset))[pData(gset)$liverOutcome == "NAFLD"]
GSE33814 <- list("ex"=ex, "meta"=pData(gset))
GSE33814$meta$GSE = 'GSE33814'

# Human hepatitis, cirrhosis
#GSE10140
# Skipped bc there were only 6k genes

# gset <- getGEO("GSE10140", GSEMatrix =TRUE, AnnotGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL5474", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# fvarLabels(gset) <- make.names(fvarLabels(gset))
# ex <- exprs(gset)
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exprs(gset) <- log2(ex) }
# 
# #Rename
# gset$liverOutcome <- 'Cir'
# gset$liverOutcome[gset$`used for prognostic prediction:ch1`=="No"] <- 'trash'
# # gset$liverOutcome[gset$`diagnosis:ch1`=="normal"] <- 'Healthy'
# 
# #Annotate and format
# ex <- exprs(gset)
# ex <- data.frame(ex)
# a <- fData(gset); ex$Gene.symbol <- a$Symbol
# ex <- ex[ex$Gene.symbol != "",]
# ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
# ex <- data.frame(ex)
# rownames(ex) <- ex$Gene.symbol
# ex <- subset(ex, select = -c(Gene.symbol))
# 
# # keep <- rownames(pData(gset))[pData(gset)$liverOutcome == "NAFLD"]
# GSE33814 <- list("ex"=ex, "meta"=pData(gset))
# GSE33814$meta$GSE = 'GSE33814'

# Human HCC
#GSE10141

# Human cirrhosis, HCC
#GSE63898
gset <- getGEO("GSE63898", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13667", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Rename
gset$liverOutcome <- 'Cir'
gset$liverOutcome[gset$source_name_ch1=="hepatocellular carcinoma"] <- 'HCC'
# gset$liverOutcome[gset$`diagnosis:ch1`=="normal"] <- 'Healthy'

#Annotate and format
ex <- exprs(gset)
ex <- data.frame(ex)
a <- sapply(fData(gset)$Gene.Symbol, function(x) trimws(unlist(strsplit(x,"///"))[1]), USE.NAMES = F)
ex$Gene.symbol <- a
ex <- ex[!is.na(ex$Gene.symbol),]
ex <- group_by(ex, Gene.symbol) %>% summarise_if(is.numeric, median)
ex <- data.frame(ex)
rownames(ex) <- ex$Gene.symbol
ex <- subset(ex, select = -c(Gene.symbol))

# keep <- rownames(pData(gset))[pData(gset)$liverOutcome == "NAFLD"]
GSE63898 <- list("ex"=ex, "meta"=pData(gset))
GSE63898$meta$GSE = 'GSE63898'

# Human HCC
#GSE102079

# Human NAFLD
#GSE134438

save.image("D:/Dropbox (SBG)/Analyzing_SCRNAseqData_Seurat/DGS_ALD/all_GSE_data.RData")

###################################
#
# Putting it all together
###################################

# super lazy, slow joining...
g1 <- intersect(rownames(GSE103580$ex),rownames(GSE48452$ex))
g2 <- intersect(rownames(GSE62232$ex), g1)
genes.final <- intersect(g2, rownames(llm.data))

# how about a better version...
genes.final <- Reduce(intersect, list(rownames(GSE103580$ex),
                                       rownames(GSE48452$ex),
                                       rownames(GSE62232$ex),
                                       rownames(GSE33814$ex),
                                       rownames(GSE49541$ex),
                                       rownames(GSE83452$ex),
                                       rownames(GSE63898$ex),
                                      rownames(llm.data)))

length(genes.final) # so 5011 common genes

ex.all <- cbind(GSE103580$ex[genes.final,],
                GSE48452$ex[genes.final,],
                GSE62232$ex[genes.final,],
                GSE33814$ex[genes.final,],
                GSE49541$ex[genes.final,],
                GSE83452$ex[genes.final,],
                GSE63898$ex[genes.final,],
                llm.data[genes.final,])

data.batch <- c(rep(1,ncol(GSE103580$ex)),
                rep(2,ncol(GSE48452$ex)),
                rep(3,ncol(GSE62232$ex)),
                rep(4,ncol(GSE33814$ex)),
                rep(5,ncol(GSE49541$ex)),
                rep(6,ncol(GSE83452$ex)),
                rep(7,ncol(GSE63898$ex)),
                rep(8,sum(grepl("o.",colnames(llm.data)))),
                rep(9,sum(grepl("y.",colnames(llm.data)))))

ex.all.norm <- ComBat(dat=as.matrix(ex.all), batch=data.batch, par.prior=TRUE, prior.plots=FALSE)
ex.all.norm <- data.frame(ex.all.norm)

ex.training <- dplyr::select(ex.all.norm, contains("GSM"))
ex.training <- data.frame(ex.training)

ex.rat <- dplyr::select(ex.all.norm, !contains("GSM"))
ex.rat <- data.frame(ex.rat)

meta.training <- rbind(GSE103580$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE48452$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE62232$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE33814$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE49541$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE83452$meta[,c("geo_accession", "liverOutcome","GSE")],
                       GSE63898$meta[,c("geo_accession", "liverOutcome","GSE")])

# might relabel "Hep" to NASH...
meta.training$liverOutcome[meta.training$liverOutcome=="Hep"] <- "NASH"

write.csv(ex.training, file="humanALD_ex.csv")
write.csv(ex.rat, file="ratLLM_ex.csv")
write.csv(meta.training, file = "humanALD_meta.csv")

plot(density(ex.all.norm[,220]))

heatmap.2(as.matrix(ex.all.norm[sample(1:nrow(ex.all.norm),50),]), 
        scale = 'row', trace='none')

# pathology numbers
table(meta.training$liverOutcome)

#DEG and viz of selected features

humanALD_ex <- read_csv("humanALD_ex.csv")
humanALD_ex <- data.frame(humanALD_ex)
row.names(humanALD_ex) <- humanALD_ex$X1
humanALD_ex <- subset(humanALD_ex, select = -c(X1))

humanALD_meta <- read_csv("humanALD_meta.csv")
humanALD_meta <- data.frame(humanALD_meta)
row.names(humanALD_meta) <- humanALD_meta$X1
humanALD_meta <- subset(humanALD_meta, select = -c(X1))
humanALD_meta$Stage <- 1
humanALD_meta$Stage[grepl('NASH|NAFLD', humanALD_meta$liverOutcome)] <- 2
humanALD_meta$Stage[grepl('Cir|HCC', humanALD_meta$liverOutcome)] <- 3
humanALD_meta$Advanced <- 1
humanALD_meta$Advanced[humanALD_meta$liverOutcome == "NASH"] <- 2
humanALD_meta$Advanced[humanALD_meta$liverOutcome == "HCC"] <- 2
humanALD_meta$Advanced <- factor(humanALD_meta$Advanced)
humanALD_meta$Stage <- factor(humanALD_meta$Stage)
humanALD_meta$Stage <- relevel(humanALD_meta$Stage, ref = "1")
humanALD_meta$Advanced <- relevel(humanALD_meta$Advanced, ref = "1")
humanALD_meta$liverOutcome <- factor(humanALD_meta$liverOutcome)
humanALD_meta$liverOutcome <- relevel(humanALD_meta$liverOutcome, ref="Healthy")

mD.interaction <- model.matrix(~Stage + Stage:Advanced, humanALD_meta)
mD.interaction <- model.matrix(~liverOutcome, humanALD_meta)
fit.path2 <- lmFit(humanALD_ex, mD.interaction)
fit.path2 <- eBayes(fit.path2)

# "(Intercept)"      "Stage2"           "Stage3"           "Stage1:Advanced2" "Stage2:Advanced2" "Stage3:Advanced2"
# "(Intercept)" "liverOutcomeCir"     "liverOutcomeHCC"     "liverOutcomeNAFLD"   "liverOutcomeNASH"
tg<-topTable(fit.path2, coef = "liverOutcomeNASH", n=Inf, p.value = 0.05)
# pathStage <- data.frame(humanALD_meta$Stage); rownames(pathStage) <- rownames(humanALD_meta)
pathStage <- data.frame(humanALD_meta$liverOutcome); rownames(pathStage) <- rownames(humanALD_meta)
pheatmap(as.matrix(humanALD_ex[sigGenes,]), scale = "row",
         show_colnames = F, cluster_cols = T,
         annotation_col = pathStage)

# significant genes
tg<-topTable(fit.path2, coef = "liverOutcomeNASH", n=100, p.value = 0.05)
sigGenes <- rownames(tg) #425
tg<-topTable(fit.path2, coef = "liverOutcomeNAFLD", n=100, p.value = 0.05)
sigGenes <- union(sigGenes, rownames(tg)) #427
tg<-topTable(fit.path2, coef = "liverOutcomeHCC", n=100, p.value = 0.05)
sigGenes <- union(sigGenes, rownames(tg)) #1755
tg<-topTable(fit.path2, coef = "liverOutcomeCir", n=100, p.value = 0.05)
sigGenes <- union(sigGenes, rownames(tg)) #2237

write.csv(sigGenes, file="sigDEgenes.csv")
write.csv(sigGenes, file="sigDEgenes_balanced.csv")
write.csv(sigGenes, file="NAFLDgenes.csv")

