#######################################################
# PART 2
# AGGREGATION AND PROCESSING OF HUMAN ALD MICROARRAY DATA
# 
# The purpose of this script is to download and clean
# various microarray datasets for the purpose of downstream ALD modeling
#
#
#######################################################

getwd()
# "D:/Dropbox (SBG)/Analyzing_SCRNAseqData_Seurat/DGS_ALD"


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

load("all_GSE_data.RData")

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