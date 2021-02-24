library(SummarizedExperiment)
library(pvca)
library(lme4)
library(multcomp)
###
###prepare psi matrix from exonic_parts.psi files###
sampleroster <- read.table("~/samples.txt", sep="\t", header = FALSE, stringsAsFactors=FALSE)
location <- paste("~/",sampleroster[1,1],".exonic_parts.psi", sep="")
samplecounts <- as.data.frame(read.table(file=location, sep="\t", header=TRUE, row.names=1, col.names=c("exon_bin", "length", "inclusion", "exclusion", "psi")))
samplecounts$psi[samplecounts$inclusion >= 10 & samplecounts$exclusion < 10] <- ceiling(samplecounts$psi[samplecounts$inclusion >= 10 & samplecounts$exclusion < 10] * 10)/10 
samplecounts$psi[samplecounts$inclusion < 10 & samplecounts$exclusion < 10] <- NA
colnames(samplecounts) <- c("length", "inclusion", "exclusion", sampleroster[1,1])
psi <- samplecounts[,4, drop=FALSE]
for (m in 2:nrow(sampleroster)){
	location <- paste("~/",sampleroster[m,1],".exonic_parts.psi", sep="")
	samplecounts <- as.data.frame(read.table(file=location, sep="\t", header=TRUE, row.names=1, col.names=c("exon_bin", "length", "inclusion", "exclusion", "psi")))
	samplecounts$psi[which(samplecounts$inclusion >= 10 & samplecounts$exclusion < 10)] <- ceiling(samplecounts$psi[samplecounts$inclusion >= 10 & samplecounts$exclusion < 10] * 10)/10
	samplecounts$psi[which(samplecounts$inclusion < 10 & samplecounts$exclusion < 10)] <- NA
	colnames(samplecounts) <- c("length", "inclusion", "exclusion", sampleroster[m,1])
	psi <- cbind(psi, samplecounts[,4, drop=FALSE])
}
##
sampleinfo <- read.csv("~/sampleinfo.csv", sep=",", row.names = 1, header=TRUE)
fullpsi <- na.omit(psi)
###PCA/PVCA for all non NA exon bins###
full.psi.pca <- prcomp(t(na.omit(fullpsi)))
full.psi.pca.proportionvariances <- ((full.psi.pca$sdev^2) / (sum(full.psi.pca$sdev^2)))*100
pct_threshold <- 0.35
batch.factors <- c("ancestry", "Location", "Disease")
full.psieset <- ExpressionSet(fullpsi,as(sampleinfo,"AnnotatedDataFrame"))
full.pvcaobj <- pvcaBatchAssess(full.psieset, batch.factors, pct_threshold)
full.pvca.res <- data.frame(label=as.character(full.pvcaobj$label),wmpv=round(as.numeric(full.pvcaobj$dat),3))
###PSI filtering###
filteredpsi <- fullpsi[ rowSums(fullpsi, na.rm=TRUE) > 0, ] 
filteredpsi <- filteredpsi[ rowSums(filteredpsi, na.rm=TRUE) < ncol(filteredpsi), ] 
remzeros <- rowSums(filteredpsi > 0.05) >= 45 
filteredpsi <- filteredpsi[remzeros,]
keep <- rowSums(filteredpsi < 0.95) >= 45 
filteredpsi <- filteredpsi[keep,]
###Filtered PSI PCA/PVCA###
filtered.psi.pca <- prcomp(t(filteredpsi))
filtered.psi.pca.proportionvariances <- ((filtered.psi.pca$sdev^2) / (sum(filtered.psi.pca$sdev^2)))*100
filtered.psieset <- ExpressionSet(filteredpsi,as(sampleinfo,"AnnotatedDataFrame"))
filtered.pvcaobj <- pvcaBatchAssess(filtered.psieset, batch.factors, pct_threshold)
filtered.pvca.res <- data.frame(label=as.character(filtered.pvcaobj$label),wmpv=round(as.numeric(filtered.pvcaobj$dat),3))
###Run linear mixed models to get significant PSI by location, ancestry, disease, & interaction###
psidftrans <- t(filteredpsi)
psiinfo <- cbind(sampleinfo[,c(3,5,6,7,9,14)], psidftrans[,1])
colnames(psiinfo) <- c("ind", "disease", "ancestry", "site", "tissue", "group", "psi")
model.full = lmer(psi ~ disease + tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
model.null = lmer(psi ~ (1|site) + (1|ind), data=psiinfo, REML=FALSE)
model.d = lmer(psi ~ tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
model.t = lmer(psi ~ disease + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
model.a = lmer(psi ~ disease + tissue + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
model.i = lmer(psi ~ disease*tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
pall <- anova(model.null,model.full)[2,8]
pd <- anova(model.d,model.full)[2,8]
pt <- anova(model.t,model.full)[2,8]
pa <- anova(model.a,model.full)[2,8]
pint <- anova(model.full,model.i)[2,8]
pvals <- matrix(c(pall,pd,pt,pa,pint),nrow=1,ncol=5)
for (m in 2:ncol(psidftrans)){
	psiinfo <- cbind(sampleinfo[,c(3,5,6,7,9,14)], psidftrans[,m])
	colnames(psiinfo) <- c("ind", "disease", "ancestry", "site", "tissue", "group", "psi")
	model.full = lmer(psi ~ disease + tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	model.null = lmer(psi ~ (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	model.d = lmer(psi ~ tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	model.t = lmer(psi ~ disease + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	model.a = lmer(psi ~ disease + tissue + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	model.i = lmer(psi ~ disease*tissue + ancestry + (1|site) + (1|ind), data=psiinfo, REML=FALSE)
	pall <- anova(model.null,model.full)[2,8]
	pd <- anova(model.d,model.full)[2,8]
	pt <- anova(model.t,model.full)[2,8]
	pa <- anova(model.a,model.full)[2,8]
	pint <- anova(model.full,model.i)[2,8]
	pvals <- rbind(pvals, c(pall,pd,pt,pa,pint))
}
rownames(pvals) <- rownames(filteredpsi)
colnames(pvals) <- c("overall", "disease", "tissue", "ancestry", "disease.location")
overallpadj <- p.adjust(pvals[,1], method="fdr")
dpadj <- p.adjust(pvals[,2], method="fdr")
tpadj <- p.adjust(pvals[,3], method="fdr")
apadj <- p.adjust(pvals[,4], method="fdr")
intpadj <- p.adjust(pvals[,5], method="fdr")
padjs <- cbind(pvals, overallpadj, dpadj, tpadj, apadj, intpadj )
###relaxed criteria for spliceopathy analysis###
spliceopathypsi <- psi[,which(sampleinfo$category =='SI')]
ileumpsi <- psi[,which(sampleinfo$category =='Ileum')]
rectumpsi <- psi[,which(sampleinfo$category =='Rectum')]
#filtering#
spliceopathypsi <- na.omit(spliceopathypsi)
ileumpsi <- ileumpsi[rownames(spliceopathypsi),]
rectumpsi <- rectumpsi[rownames(spliceopathypsi),]
ileumpsi <- ileumpsi[ rowSums(is.na(ileumpsi)) < 6, ]
rectumpsi <- rectumpsi[rownames(ileumpsi),]
rectumpsi <- rectumpsi[ rowSums(is.na(rectumpsi)) < 6, ]
ileumpsi <- ileumpsi[rownames(rectumpsi),]
for (i in 1:nrow(ileumpsi)){
	ileumpsi[i,][is.na(ileumpsi[i,])] <- mean(ileumpsi[i,], na.rm=TRUE)
}
for (i in 1:nrow(rectumpsi)){
	rectumpsi[i,][is.na(rectumpsi[i,])] <- mean(rectumpsi[i,], na.rm=TRUE)
}
spliceopathypsi <- spliceopathypsi[rownames(ileumpsi),]
sripsi <- cbind(ileumpsi, rectumpsi, spliceopathypsi)
sriinfo <- rbind(sampleinfo[which(sampleinfo$category =='Ileum'),],sampleinfo[which(sampleinfo$category =='Rectum'),],sampleinfo[which(sampleinfo$category =='SI'),])
sripsi <- sripsi[ rowSums(sripsi, na.rm=TRUE) > 0, ]
sripsi <- sripsi[ rowSums(sripsi, na.rm=TRUE) < ncol(sripsi), ]
sriremzeros <- rowSums(sripsi > 0.05) >= 45 
sripsi <- sripsi[sriremzeros,]
srikeep <- rowSums(sripsi < 0.95) >= 45 
sripsi <- sripsi[srikeep,]
###Run linear mixed models to get spliceopathy PSIs significant from Ileum###
sripsitrans <- t(sripsi2)
sripsiinfo <- cbind(sriinfo[,c(3,5,6,7,9,14)], sripsitrans[,1])
colnames(sripsiinfo) <- c("ind", "disease", "ancestry", "site", "tissue", "group", "psi")
sri.model.full = lmer(psi ~ group + ancestry + (1|site) + (1|ind), data=sripsiinfo, REML=FALSE)
sri.model.null = lmer(psi ~ ancestry + (1|site) + (1|ind), data=sripsiinfo, REML=FALSE)
sri.model <- lmer(psi ~ group + ancestry + (1|site) + (1|ind), data=sripsiinfo)
sripcomps <- summary(glht(sri.model,mcp(group="Tukey")))$test$pvalues[1:3]
sripall <- anova(sri.model.null,sri.model.full)[2,8]
sripvals <- matrix(c(sripall,sripcomps),nrow=1,ncol=4)
for (m in 2:ncol(sripsitrans)){
	sripsiinfo <- cbind(sriinfo[,c(3,5,6,7,9,14)], sripsitrans[,m])
	colnames(sripsiinfo) <- c("ind", "disease", "ancestry", "site", "tissue", "group", "psi")
	sri.model.full = lmer(psi ~ group + ancestry + (1|site) + (1|ind), data=sripsiinfo, REML=FALSE)
	sri.model.null = lmer(psi ~ ancestry + (1|site) + (1|ind), data=sripsiinfo, REML=FALSE)
	sripall <- anova(sri.model.null,sri.model.full)[2,8]
	if (sripall <= 0.05){
		sri.model <- lmer(psi ~ group + ancestry + (1|site) + (1|ind), data=sripsiinfo)
		sripcomps <- summary(glht(sri.model,mcp(group="Tukey")))$test$pvalues[1:3]
	} else {
		sripcomps <- c(1,1,1)
	}
	sripvals <- rbind(sripvals, c(sripall,sripcomps))
}
rownames(sripvals) <- rownames(sripsi)
sri.padjall <- p.adjust(sripvals[,1], method="fdr")
R.I <- p.adjust(sripvals[,2], method="fdr")
S.I <- p.adjust(sripvals[,3], method="fdr")
S.R <- p.adjust(sripvals[,4], method="fdr")
sri.padjs <- cbind(sripvals, sri.padjall, R.I, S.I, S.R)
rownames(sri.padjs) <- rownames(sripsi)
colnames(sri.padjs) <- c("pvaloverall", "Re.Il", "Sp.Il", "Sp.Re", "padjoverall","R.I.padj", "S.I.padj", "S.R.padj")
