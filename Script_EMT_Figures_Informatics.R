########################################################################
########################################################################
##     Paper: Statin-induced cancer cell death can be mechanistically
##            uncoupled from protein prenylation
##     
##     Code to generate Figures 5A,B,C,D,E and S5A,B,C
##     Bimodality analysis and EMT association with Statins' sensitivity
##
##     Written by: Wail Ba alawi
##     Date: 17 Feb 2016
########################################################################
########################################################################





library(Biobase)
library(PharmacoGx)
library(mixtools)
library(mclust)



##################
## set working directory to where this script is downloaded

StatinEMT_Directory <- NULL
#change the above line to:
# StatinEMT_Directory <- "./pathToStatinEMT-master_Directory"


if(!is.null(StatinEMT_Directory)){
  setwd(StatinEMT_Directory)

##################

CCLE <- downloadPSet("CCLE",saveDir = "./data/")
CTRPv2 <- downloadPSet("CTRPv2",saveDir = "./data/")




#load("./data/CCLE.RData")
#load("./data/CTRPv2.RData")





#####################################
#####################################
# Merge the two data sets [CTRPv2 and CCLE] (Code by Benjamin)

### fix lapatinib name
#drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pattern="N-{3-Chloro-4-[(3-fluorobenzyl)oxy]phenyl}-6-[5-({[2-(methylsulfonyl)ethyl]amino}methyl)-2-furyl]-4-quinazolinamine", replacement="lapatinib", fixed=TRUE)
rownames(CTRPv2@curation$tissue) <- rownames(CTRPv2@curation$cell)
CTRPv2@curation$tissue[!is.na(CTRPv2@curation$tissue[ , "unique.tissueid"]) & CTRPv2@curation$tissue[ , "unique.tissueid"] == "", "unique.tissueid"] <- NA
cellInfo(CTRPv2)[!is.na(cellInfo(CTRPv2)[ , "tissueid"]) & cellInfo(CTRPv2)[ , "tissueid"] == "", "tissueid"] <- NA


### molecular profiles
CTRPv2.CCLE <- CTRPv2
CTRPv2.CCLE@molecularProfiles <- CCLE@molecularProfiles

ucell <- union(rownames(cellInfo(CTRPv2)), rownames(cellInfo(CCLE)))

### cell line annotations
new.cell.info <- c(union(colnames(cellInfo(CTRPv2)), colnames(cellInfo(CCLE))), "dataset")
cell.info.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.info), dimnames=list(ucell, new.cell.info)))
cell.info.df[rownames(cellInfo(CCLE)), colnames(cellInfo(CCLE))] <- cellInfo(CCLE)
cell.info.df[rownames(cellInfo(CTRPv2)), colnames(cellInfo(CTRPv2))] <- cellInfo(CTRPv2)
cell.info.df[setdiff(rownames(cellInfo(CTRPv2)), rownames(cellInfo(CCLE))), "dataset"] <- "CTRPv2"
cell.info.df[setdiff(rownames(cellInfo(CCLE)), rownames(cellInfo(CTRPv2))), "dataset"] <- "CCLE"
cell.info.df[intersect(rownames(cellInfo(CCLE)), rownames(cellInfo(CTRPv2))), "dataset"] <- "CTRPv2///CCLE"
cellInfo(CTRPv2.CCLE) <- cell.info.df

### curation of cell line names
new.cell.curation <- c(union(colnames(CTRPv2@curation$cell), colnames(CCLE@curation$cell)), "dataset")
cell.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.curation), dimnames=list(ucell, new.cell.curation)))
cell.curation.df[rownames(CCLE@curation$cell), colnames(CCLE@curation$cell)] <- CCLE@curation$cell
cell.curation.df[rownames(CTRPv2@curation$cell), colnames(CTRPv2@curation$cell)] <- CTRPv2@curation$cell
cell.curation.df[setdiff(rownames(CTRPv2@curation$cell), rownames(CCLE@curation$cell)), "dataset"] <- "CTRPv2"
cell.curation.df[setdiff(rownames(CCLE@curation$cell), rownames(CTRPv2@curation$cell)), "dataset"] <- "CCLE"
cell.curation.df[intersect(rownames(CCLE@curation$cell), rownames(CTRPv2@curation$cell)), "dataset"] <- "CTRPv2///CCLE"
CTRPv2.CCLE@curation$cell <- cell.curation.df

### curation of tissue names
new.tissue.curation <- c(union(colnames(CTRPv2@curation$tissue), colnames(CCLE@curation$tissue)), "dataset")
tissue.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.tissue.curation), dimnames=list(ucell, new.tissue.curation)))
tissue.curation.df[rownames(CCLE@curation$tissue), colnames(CCLE@curation$tissue)] <- CCLE@curation$tissue
tissue.curation.df[rownames(CTRPv2@curation$tissue), colnames(CTRPv2@curation$tissue)] <- CTRPv2@curation$tissue
tissue.curation.df[setdiff(rownames(CTRPv2@curation$tissue), rownames(CCLE@curation$tissue)), "dataset"] <- "CTRPv2"
tissue.curation.df[setdiff(rownames(CCLE@curation$tissue), rownames(CTRPv2@curation$tissue)), "dataset"] <- "CCLE"
tissue.curation.df[intersect(rownames(CCLE@curation$tissue), rownames(CTRPv2@curation$tissue)), "dataset"] <- "CTRPv2///CCLE"
CTRPv2.CCLE@curation$tissue <- tissue.curation.df

###############################
##############################
###############################
# intializing EMT related variables

coreSet <- c("ENSG00000039068", "ENSG00000170558", "ENSG00000115414", "ENSG00000026025", "ENSG00000148516")
names(coreSet) <- c("CDH1", "CDH2", "FN1", "VIM","ZEB1")


# FULL expanded EMT set 
expandedSet_Full <- c("ENSG00000039068", "ENSG00000170558", "ENSG00000115414", "ENSG00000026025", "ENSG00000176692", "ENSG00000019549", "ENSG00000124216", "ENSG00000148516", "ENSG00000122691", "ENSG00000169554", "ENSG00000233608", "ENSG00000133937", "ENSG00000149948")
names(expandedSet_Full) <- c("CDH1", "CDH2", "FN1", "VIM", "FOXC2", "SNAI2", "SNAI1", "ZEB1", "TWIST1", "ZEB2", "TWIST2", "GSC", "HMGA2")

expandedSet_Full_Direction <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
names(expandedSet_Full_Direction) <- expandedSet_Full

# control genes for bimodality (POLR2A: Negative, ESR1:Positive)
control <- c("ENSG00000181222", "ENSG00000091831")
names(control) <- c("POLR2A", "ESR1")

j <- grep("statin", drugNames(CTRPv2))
statins <- drugNames(CTRPv2)[j][c(1,2,6)]




############################################
############################################
############################################

## Function to compute BiModality Index (adopted from genefu)
getBiModalScore <- function(data=data){
  
  data <- data[!is.na(data)]
  if(var(data)==0){
    return(0)
  }
  rr2 <- mclust::Mclust(data = data, modelNames="V", G=2)
  res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
  
  if(is.null(rr2[[1]])) { ## EM algorithm did not converge
   # message("EM didn't converage")
    return(0)
  }
  res[1, ] <- rr2$parameters$mean
  res[2, ] <- rr2$parameters$variance$sigmasq
  res[3, ] <- rr2$parameters$pro
  
  ## bimodality index (BI)
  smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2] + res[2, 1]) / 2)
  bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd
  
  return(bi)
  
}



# retrieve all RNAseq expression data from CCLE for all cell lines
RNAseq_CTRPv2_All_CCLE <- t(exprs(summarizeMolecularProfiles(CTRPv2.CCLE,mDataType = "rnaseq",fill.missing = F)))

#################
#### NOTE: The following command will take some time to finish. If you want to skip then comment the coming four commands and uncomment the 5th one ####
################

# calculate Bimodality index for all genes using RNA-seq expression of each gene across all cell lines in CCLE
#BiModalScores_CTRPv2_All_CCLE <- apply(FUN = getBiModalScore,MARGIN = 2,X = RNAseq_CTRPv2_All_CCLE)
#BiModalScores_CTRPv2_All_CCLE <- BiModalScores_CTRPv2_All_CCLE[order(BiModalScores_CTRPv2_All_CCLE,decreasing = T)]
#BiModalScores_CTRPv2_CCLE_expandedSet_Full <- BiModalScores_CTRPv2_All_CCLE[expandedSet_Full]
#BiModalScores_CTRPv2_CCLE_expandedSet_Full <- BiModalScores_CTRPv2_CCLE_expandedSet_Full[order(BiModalScores_CTRPv2_CCLE_expandedSet_Full)]

load("./data/BimodalScores_All_Genes.RData")

#save(list =  c("BiModalScores_CTRPv2_All_CCLE", "BiModalScores_CTRPv2_CCLE_expandedSet_Full"),file = "./BimodalScores_All_Genes.RData")

# Fig 5A data: Bimodal index for all the genes (based on expression level) 

pdf(file = "./BimodalityScoresForAllGenes.pdf",width = 12, height = 7)
par(font=2,font.axis=2,font.lab=2,lwd=4,cex=1.5)
hist(BiModalScores_CTRPv2_All_CCLE,xlab = "Bimodal Index (BI)",freq = T,breaks = 30,col = "gray",main = NULL,xlim = c(0,3.5))
axis(side = 1, lwd = 4)
axis(side = 2, lwd = 4)
Flip <- T
for(i in 1:length(expandedSet_Full)){
  
  if(names(BiModalScores_CTRPv2_CCLE_expandedSet_Full)[i] %in% expandedSet_Full & BiModalScores_CTRPv2_CCLE_expandedSet_Full[i]!=0){ 
    points(x = BiModalScores_CTRPv2_CCLE_expandedSet_Full[i],y=-550,pch=20,col="darkgreen")
    
    if(names(BiModalScores_CTRPv2_CCLE_expandedSet_Full)[i] %in% coreSet){
      points(x = BiModalScores_CTRPv2_CCLE_expandedSet_Full[i],y=-550,pch=20,col="blue")
    }
    
    Flip <- !Flip
  }
}
points(x=BiModalScores_CTRPv2_All_CCLE[control],y=c(rep(-550,length(control))),pch=18)


legend(x = c(2.2,3.1), y=c(36500,41000),legend = c("EMT related genes", "Top five EMT related genes", "Other genes"),col = c("darkgreen","blue","black"),bty="n", pch = c(20,20,18))
dev.off()



############################################
############################################
############################################


# order the top 5 EMT genes' set by Bimodality index
newSet <- coreSet

# reorder top 5 EMT genes by their bimodality score 
newSet <- newSet[c(4,1,5,3,2)]


# Figure 5A: plots for top 5 EMT genes bimodality plus the controls
pdf("./Fig_5A.pdf", height = 17, width = 13)
par(mfrow=c(4,2))

cutoffs <- NULL
for(i in 1:length(newSet)){
  
  data <- RNAseq_CTRPv2_All_CCLE[,newSet[i]][!is.na(RNAseq_CTRPv2_All_CCLE[,newSet[i]])]
  if(var(data)==0){
    cutoffs <- c(cutoffs,max(data,na.rm = T))
    next
  }
  rr2 <- mclust::Mclust(data = data, modelNames="V", G=2)
  if(is.null(rr2[[1]])) { ## EM algorithm did not converge
    message("didn't converage")
    cutoffs <- c(cutoffs,max(data,na.rm = T))
    next
  }
  cells_1 <- names(rr2$classification)[rr2$classification==1] 
  cells_2 <- names(rr2$classification)[rr2$classification==2] 
  
  
  cutoff_tmp <- mean(c(min(data[cells_2]), max(data[cells_1])))
  cutoffs <- c(cutoffs,cutoff_tmp)
  
  symbol <- names(expandedSet_Full)[expandedSet_Full==newSet[i]]
  A1 <- density(data[cells_1],adjust=2)
  max.y_tmp <- max(A1$y)
  A <- density(data)
  max.y <- max(A$y)
  par(lwd = 4)
  foo <- hist(data,breaks = 40,freq = F,axes = F, ylab = NULL, xlab = NULL,col = "lightgray", main = NULL,cex=2, xaxt='n',ann=FALSE)
  max.y <- max(foo$density)
  #  axis(side = 2, lwd = 4,ann=FALSE)
  title(main = symbol,col.main="blue",cex.main=4)
  
  max.y <- max(foo$density)
  A1 <- density(data[cells_1],adjust=2)
  
  lines(A1$x,scales::rescale(A1$y, to=c(0,max.y)),col="indianred4",lwd=4)
  A2 <- density(data[cells_2],adjust=2)
  max.y <- max(A$y)
  lines(A2$x,scales::rescale(A2$y, to=c(0,max.y)), col="mediumorchid4",lwd=4)
  abline(v=cutoff_tmp, col=c("red"), lty=2)
  
}

# positive control

data <- RNAseq_CTRPv2_All_CCLE[,control[2]][!is.na(RNAseq_CTRPv2_All_CCLE[,control[2]])]
if(var(data)==0){
  cutoffs <- c(cutoffs,max(data,na.rm = T))
  next
}
rr2 <- mclust::Mclust(data = data, modelNames="V", G=2)
if(is.null(rr2[[1]])) { ## EM algorithm did not converge
  message("didn't converage")
  cutoffs <- c(cutoffs,max(data,na.rm = T))
  next
}
cells_1 <- names(rr2$classification)[rr2$classification==1] 
cells_2 <- names(rr2$classification)[rr2$classification==2] 


cutoff_tmp <- mean(c(min(data[cells_2]), max(data[cells_1])))
cutoffs <- c(cutoffs,cutoff_tmp)

symbol <- names(control)[2]
A1 <- density(data[cells_1],adjust=2)
max.y_tmp <- max(A1$y)
A <- density(data)
max.y <- max(A$y)
par(lwd = 4)
foo <- hist(data,breaks = 50,freq = F,axes = F, ylab = NULL, xlab = NULL,col = "lightgray", main = NULL,cex=2, xaxt='n',ann=FALSE)
max.y <- max(foo$density)
#  axis(side = 2, lwd = 4,ann=FALSE)
title(main = symbol,col.main="blue",cex.main=4)

max.y <- max(foo$density)
A1 <- density(data[cells_1],adjust=2)

lines(A1$x,scales::rescale(A1$y, to=c(0,max.y)),col="indianred4",lwd=4)
A2 <- density(data[cells_2],adjust=2)
max.y <- max(A$y)
lines(A2$x,scales::rescale(A2$y, to=c(0,max.y/4)), col="mediumorchid4",lwd=4)
abline(v=cutoff_tmp, col=c("red"), lty=2)




# negative control
data <- RNAseq_CTRPv2_All_CCLE[,control[1]][!is.na(RNAseq_CTRPv2_All_CCLE[,control[1]])]
rr2 <- mclust::Mclust(data = data, modelNames="V", G=2)
if(is.null(rr2[[1]])) { ## EM algorithm did not converge
  message("didn't work")
  cutoffs <- c(cutoffs,max(data,na.rm = T))
  next
}
cells_1 <- names(rr2$classification)[rr2$classification==1] 
cells_2 <- names(rr2$classification)[rr2$classification==2] 
#foo <- hist((data[cells_1]),ylab = NULL, xlab = NULL)
#plot(foo$mids, foo$counts, type="l") 
#lines (density(data[cells_2]))

cutoff_tmp <- mean(c(min(data[cells_2]), max(data[cells_1])))
symbol <- "POLR2A"
A <- density(data)
max.y <- max(A$y)
par(lwd = 4)
foo <- hist(data,breaks = 40,freq = F,axes = F, ylab = NULL, xlab = NULL,col = "lightgray", main = NULL,cex=4,ylim = c(0,max.y+0.1))
#axis(side = 2, lwd = 4,ann=T)
title(main = symbol,col.main="blue",cex.main=4)
max.y <- max(foo$density)

A2 <- density(data[cells_2],adjust=2)
max.y <- max(A$y)
lines(A2$x,scales::rescale(A2$y, to=c(0,max.y)), col="mediumorchid4",lwd=4)
#abline(v=cutoff_tmp, col=c("red"), lty=2)

names(cutoffs) <- c(newSet,control[2])

dev.off()



# legends for Figure A
pdf("./Fig_5A_legends.pdf", height = 3, width = 7)
par(font=2)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, legend = c("Non-expressed cell lines","Expressed cell lines","cutoff"), col=c("indianred4","mediumorchid4","red"),lwd=c(2,2,2), cex=2, xjust=0.5, yjust=0.5,bty="n",lty = c(1,1,2))
dev.off()


############################################
############################################
############################################



# Figure 5BCD: top 5 EMT genes association with Statin response


newSet_Dir <- expandedSet_Full_Direction[newSet]
CIs <- list()
pvalues <- list()
pdf("./Fig_5BCD.pdf",height = 18,width = 6)
par(mfrow=c(3,1))
for(drug in statins){
  
  cellsToKeep <- rownames(CTRPv2.CCLE@cell)[CTRPv2.CCLE@cell$tissueid!="haematopoietic_and_lymphoid_tissue"]
  cellsToKeep <- intersect(cellsToKeep, unique(CTRPv2.CCLE@sensitivity$info$cellid[CTRPv2.CCLE@sensitivity$info$drugid==drug]) )
  CTRPv2_test_subset <- subsetTo(CTRPv2.CCLE, drugs = drug,cells = cellsToKeep)
  
  RNAseq_CTRPv2_Drug <- t(exprs(summarizeMolecularProfiles(CTRPv2_test_subset,mDataType = "rnaseq",fill.missing = F)))
  RNAseq_CTRPv2_newSet <- RNAseq_CTRPv2_Drug[,newSet]
  
  for( gene in newSet){ 
    RNAseq_CTRPv2_newSet[,gene] <- ifelse(RNAseq_CTRPv2_newSet[,gene]>cutoffs[gene],1,0)
  }
  label <- NULL
  label <- apply(RNAseq_CTRPv2_newSet,MARGIN = 1,FUN = function(x){
    
    if ( any(as.numeric(!xor(x , expandedSet_Full_Direction[newSet]))==1,na.rm = T)){
      return(1)
    }else{
      return(0)
    }
  })
  
  RNAseq_CTRPv2_newSet <- cbind(RNAseq_CTRPv2_newSet, "label"=label)
  
  AUC_subset <- summarizeSensitivityProfiles(CTRPv2_test_subset,sensitivity.measure = "auc_recomputed",fill.missing = F)*100
  
  cellLines_intersect <- intersect(rownames(RNAseq_CTRPv2_newSet),names(AUC_subset))
  
  
  expressedCellLines <- rownames(RNAseq_CTRPv2_newSet)[RNAseq_CTRPv2_newSet[,"label"]==1]
  expressedCellLines <- intersect(cellLines_intersect,expressedCellLines)
  notExpressedCellLines <- rownames(RNAseq_CTRPv2_newSet)[RNAseq_CTRPv2_newSet[,"label"]==0]
  notExpressedCellLines <- intersect(cellLines_intersect,notExpressedCellLines)
  par(font.axis = 2,font.lab=2,cex=1.5, lwd=4)
  boxplot(cbind("Not Enriched"=AUC_subset[notExpressedCellLines],"Enriched"=AUC_subset[expressedCellLines]),ylab="AUC [%]",col=c("indianred2","lightskyblue"))
  
  integrCindex <- Hmisc::rcorr.cens(S=AUC_subset[c(notExpressedCellLines,expressedCellLines)], x = as.numeric(RNAseq_CTRPv2_newSet[c(notExpressedCellLines,expressedCellLines),"label"]),outx = T )
  C <- integrCindex["C Index"]

  Pvalue <-kruskal.test(list(AUC_subset[notExpressedCellLines], AUC_subset[expressedCellLines]))$p.value

  title(main = drug)
  
  par(font=2,cex=1.3)
  legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("n:",length(notExpressedCellLines)),paste("n:",length(expressedCellLines)) ),cex = 1,bty="n",pch = c(NA,NA,15,15),col=c("indianred2","lightskyblue") )
}
dev.off()


###############################
##############################
###############################
# Figure S5A Tissue Types in CCLE rnaseq data

# retrieve all distribution of all tissue types in CCLE 
tissues <- table(CTRPv2.CCLE@curation$tissue[,"unique.tissueid"])

tissues <- as.array(tissues)

tissues <- tissues[names(tissues)!=""]
names(tissues) <- sub("_", " ", names(tissues))
names(tissues) <- sub("_", " ", names(tissues))
names(tissues) <- sub("_", " ", names(tissues))

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

names(tissues) <- firstup(names(tissues))

# we combined small and large intestine in one category for ease of presentation
tissues["Small and large intestine"] <- sum(tissues["Large intestine"]+tissues["Small intestine"])

tissues <- tissues[names(tissues)!="Large intestine"]
tissues <- tissues[names(tissues)!="Small intestine"]


cols2 <- c(
  "#a6cee3"
  ,"#1f78b4"
  ,"#b2df8a"
  ,"#33a02c"
  ,"#fb9a99"
  ,"#e31a1c"
  ,"#fdbf6f"
  ,"#ff7f00"
  ,"#cab2d6"
  ,"#6a3d9a"
  ,"#ffff99"
  ,"#b15928"
)


# rearrange the order of the tissues for better presentation 
tissues <- tissues[c(1,2,3,4,16,5,15,7,23,14,19,11,17,13,6,9,8,18,10,20,21,22,12)]
pdf("./Fig_S5A.pdf",height = 8, width = 11)
par(font=2,lwd=4)
pie(tissues, col = cols2,main = NULL,cex=1.5)
dev.off()


########################################
# Figure S5B bimodal distribution of ESR1 in breast cancer cell lines. 

pdf("./Fig_S5B.pdf",width = 6,height = 6)
BC_cellines <- rownames(CTRPv2.CCLE@cell)[which(CTRPv2.CCLE@cell$tissueid == "breast")]
BC_cellines <- intersect(BC_cellines,rownames(RNAseq_CTRPv2_All_CCLE))
data <- RNAseq_CTRPv2_All_CCLE[BC_cellines,control[2]][!is.na(RNAseq_CTRPv2_All_CCLE[BC_cellines,control[2]])]

rr2 <- mclust::Mclust(data = data, modelNames="V", G=2)

cells_1 <- names(rr2$classification)[rr2$classification==1] 
cells_2 <- names(rr2$classification)[rr2$classification==2] 


cutoff_tmp <- mean(c(min(data[cells_2]), max(data[cells_1])))

symbol <- names(control)[2]
A1 <- density(data[cells_1],adjust=2)
max.y_tmp <- max(A1$y)
A <- density(data)
max.y <- max(A$y)
par(lwd = 4)
foo <- hist(data,breaks = 50,freq = F,axes = F, ylab = NULL, xlab = NULL,col = "lightgray", main = NULL,cex=2, xaxt='n',ann=FALSE)
max.y <- max(foo$density)
#  axis(side = 2, lwd = 4,ann=FALSE)
title(main = paste(symbol,"[Breast only]",sep = " "),col.main="blue",cex.main=1)

max.y <- max(foo$density)
A1 <- density(data[cells_1],adjust=2)

lines(A1$x,scales::rescale(A1$y, to=c(0,max.y)),col="indianred4",lwd=4)
A2 <- density(data[cells_2],adjust=2)
max.y <- max(A$y)
lines(A2$x,scales::rescale(A2$y, to=c(0,max.y*2)), col="mediumorchid4",lwd=4)
abline(v=cutoff_tmp, col=c("red"), lty=2)

dev.off()



###############################
##############################
###############################
# Figure S5DEF ESR1 association with Statin response


CIs <- list()
pvalues <- list()
pdf("./Fig_S5DEF.pdf",height = 23,width = 6)
par(mfrow=c(3,1))
for(drug in statins){
  
  cellsToKeep <- rownames(CTRPv2.CCLE@cell)[CTRPv2.CCLE@cell$tissueid!="haematopoietic_and_lymphoid_tissue"]
  cellsToKeep <- intersect(cellsToKeep, unique(CTRPv2.CCLE@sensitivity$info$cellid[CTRPv2.CCLE@sensitivity$info$drugid==drug]) )
  CTRPv2_test_subset <- subsetTo(CTRPv2.CCLE, drugs = drug,cells = cellsToKeep)
  
  RNAseq_CTRPv2_Drug <- t(exprs(summarizeMolecularProfiles(CTRPv2_test_subset,mDataType = "rnaseq",fill.missing = F)))
  RNAseq_CTRPv2_newSet <- RNAseq_CTRPv2_Drug[,names(cutoffs)[6],drop=F]
  
  RNAseq_CTRPv2_newSet[,1] <- ifelse(RNAseq_CTRPv2_newSet[,1,drop=F]>cutoffs[6],1,0)
  
  label <- RNAseq_CTRPv2_newSet[,1]
  
  RNAseq_CTRPv2_newSet <- cbind(RNAseq_CTRPv2_newSet, "label"=label)
  
  AUC_subset <- summarizeSensitivityProfiles(CTRPv2_test_subset,sensitivity.measure = "auc_recomputed",fill.missing = F)*100
  
  cellLines_intersect <- intersect(rownames(RNAseq_CTRPv2_newSet),names(AUC_subset))
  
  
  expressedCellLines <- rownames(RNAseq_CTRPv2_newSet)[RNAseq_CTRPv2_newSet[,"label"]==1]
  expressedCellLines <- intersect(cellLines_intersect,expressedCellLines)
  notExpressedCellLines <- rownames(RNAseq_CTRPv2_newSet)[RNAseq_CTRPv2_newSet[,"label"]==0]
  notExpressedCellLines <- intersect(cellLines_intersect,notExpressedCellLines)
  par(font.axis = 2,font.lab=2,cex=1.5, lwd=4)
  boxplot(list("Not Enriched"=AUC_subset[notExpressedCellLines],"Enriched"=AUC_subset[expressedCellLines])
          ,ylab="AUC [%]",col=c("indianred2","lightskyblue"),ylim=c(0,100)
  )
  
  integrCindex <- Hmisc::rcorr.cens(S=AUC_subset[c(notExpressedCellLines,expressedCellLines)], x = as.numeric(RNAseq_CTRPv2_newSet[c(notExpressedCellLines,expressedCellLines),"label"]),outx = T )
  C <- integrCindex["C Index"]
  
  Pvalue <-kruskal.test(list(AUC_subset[notExpressedCellLines], AUC_subset[expressedCellLines]))$p.value
  
  title(main = drug)
  
  par(font=2,cex=1.3)
  legend("topleft",legend = c(paste("CI:", sprintf("%.2g",C)),paste("P-value:", sprintf("%.1E",Pvalue)),paste("n:",length(notExpressedCellLines)),paste("n:",length(expressedCellLines)) ),cex = 1,bty="n",pch = c(NA,NA,15,15),col=c("indianred2","lightskyblue") )
}
dev.off()





}else{
  stop("Please set the working directory to where this script is downloaded\nInstructions are at the beginning of the script file")
}




