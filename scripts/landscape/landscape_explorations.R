library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)

annot=read.delim("~/metadata/annotations/gtex_tissue_colours.txt",as.is=T) 
data=readRDS("2024_04_28_gt_rna_var_annotated.rds")
data=as.data.frame(data)
data$tissue2=data$tissue
data$annot=annot$annotation[match(data$tissue2,annot$tissue2)]
data$tissue=annot$tissue[match(data$tissue2,annot$tissue2)]
data$het_id=paste0("mt.",data$Pos)
data$key1=paste(data$SUBJID,data$tissue,data$het_id,sep=".")
data$i=1

## get number of multi-allelic sites 
d1=data.frame(table(data$key1))
data$nhetalleles=d1$Freq[match(data$key1,d1$Var1)]

## get number of het sites and n tissues per person 
d1=ddply(data,.(het_id,tissue,nhetalleles,molecular_process),summarise,count=sum(i))
d1$nind=d1$count
d1$nind[which(d1$nhetalleles==2)]=d1$count[which(d1$nhetalleles==2)]/2
d1$nind[which(d1$nhetalleles==3)]=d1$count[which(d1$nhetalleles==3)]/3
d1$key2=paste(d1$tissue,d1$het_id,sep=".")
d1$keep="NO"
d1$keep[which(d1$nind>=5)]="YES"

## plot all apparent heteroplasmy across all tissues and individuals ------------------------------------
d1$i=1
d1a=ddply(d1,.(tissue,keep),summarise,count=sum(i))
d1a=d1a[which(d1a$keep=="YES"),]
d1a=d1a[order(d1a$count,decreasing=T),]
min(d1a$count) #[1] 28
max(d1a$count) #[1] 331
d1b=ddply(d1,.(tissue,keep,molecular_process),summarise,count=sum(i))
## remove tissue hets that have fewer than 5 inds
d1b=d1b[which(d1b$keep=="YES"),]
d1b$molecular_process[which(d1b$molecular_process=="dna_mutation")]="mtDNA_het"
d1b$molecular_process[which(d1b$molecular_process=="rna_modification")]="mtRNA_mod"
d1b$molecular_process=factor(d1b$molecular_process,levels=c("mtDNA_het","mtRNA_mod"))
d1b$tissue=factor(d1b$tissue,levels=d1a$tissue)
pdf("allhets.pdf",height=5, width=10,useDingbats=FALSE)
p=ggplot(d1b, aes(tissue, count, fill = molecular_process)) + geom_col(colour="black") + theme_bw() +ylab("N apparent heteroplasmies")+xlab("Tissue")
p=p+scale_fill_manual(values=hetcols)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## plot number of alleles for mtDNA hets and mtRNA mods --------------------------------------------------
d2=data.frame(table(d1$nhetalleles,d1$molecular_process))
colnames(d2)=c("nalleles","molecular_process","count")
d2$molecular_process=as.character(d2$molecular_process)
d2$molecular_process[which(d2$molecular_process=="dna_mutation")]="mtDNA_het"
d2$molecular_process[which(d2$molecular_process=="rna_modification")]="mtRNA_mod"
d2$molecular_process=factor(d2$molecular_process)
cols=c("#edf8b1","#7fcdbb","#1d91c0")
pdf("allhets_nalleles.pdf",height=5, width=3,useDingbats=FALSE)
p=ggplot(d2, aes(molecular_process,count,fill=nalleles)) + geom_col(colour="black") + theme_bw() +ylab("Count")+xlab("Molecular process")
p=p+scale_fill_manual(values=cols) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## remove tissue hets that have fewer than 5 inds from future analyses
data$key2=paste(data$tissue,data$het_id,sep=".")
data=data[which(data$key2%in%d1$key2[which(d1$keep=="YES")]),]

## plot boxplots of individual mean varhets per tissue --------------------------------------------------
d2a=ddply(data,.(SUBJID,tissue,molecular_process,nhetalleles),summarise,count=sum(i))
d2a$nhet=d2a$count
d2a$nhet[which(d2a$nhetalleles==2)]=d2a$count[which(d2a$nhetalleles==2)]/2
d2a$nhet[which(d2a$nhetalleles==3)]=d2a$count[which(d2a$nhetalleles==3)]/3
d2b=ddply(d2a,.(SUBJID,tissue,molecular_process),summarise,nhet=sum(nhet))
d2b$annot=annot$annotation[match(d2b$tissue,annot$tissue)]
d2b$tissue=factor(d2b$tissue,levels=d1a$tissue)
annotcols=unique(annot$colours)
d2b$annot=factor(d2b$annot,levels=annot$annot[match(annotcols,annot$colours)])
pdf("allhets_alltissues.pdf",height=8, width=10,useDingbats=FALSE)
p=ggplot(d2b, aes(tissue,nhet,fill=annot)) + geom_boxplot()+ facet_wrap(vars(molecular_process),nrow=2,scales="free_y")
p=p+theme_bw()+ylab("N apparent heteroplasmy per sample")+xlab("Tissue")
p=p+scale_fill_manual(values=annotcols) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## now only work with mtDNA heteroplasmies -------------------------------------------------------------
## get somatic vs inherited ----------------------------------------------------------------------------
data1=data[which(data$molecular_process=="dna_mutation"),]
## get number of tissue types  each pos is found in 
d3=ddply(data1,.(SUBJID,het_id,annot),summarise,count=sum(i))
d3$i =1 ## not taking into account multi-allelic
d3a=ddply(d3,.(SUBJID,het_id),summarise,count_annot=sum(i))
d3b=ddply(data1,.(SUBJID,het_id,nhetalleles),summarise,count=sum(i))
d3b$nhet=d3b$count
d3b$nhet[which(d3b$nhetalleles==2)]=d3b$count[which(d3b$nhetalleles==2)]/2
d3b$nhet[which(d3b$nhetalleles==3)]=d3b$count[which(d3b$nhetalleles==3)]/3
d3b$key=paste(d3b$SUBJID,d3b$het_id,sep=".")
d3c=data.frame(table(d3a$SUBJID,d3a$type))
d3d=ddply(data1,.(SUBJID),summarise,nhet=length(unique(het_id)))
d3d=d3d[order(d3d$nhet,decreasing=T),]
d3c$Var1=factor(d3c$Var1,levels=d3d$SUBJID)
d3c$Var2=factor(d3c$Var2,levels=c("Inherited","Somatic"))

## plot number of somatic/inherited variations per donor 
pdf("mthets_perindtypes.pdf",height=5, width=5,useDingbats=FALSE)
p=ggplot(d3c, aes(Var1,Freq,fill=Var2)) + geom_col(position="stack") + theme_bw() +ylab("mtDNA heteroplasmy per donor across tissue")+xlab("Donor")
p=p+scale_fill_manual(values=cols) +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## get number of somatic/inherited variations per gene -------------------------------------------------
data1$key3=paste(data1$SUBJID,data1$het_id,sep=".")
data1$mtdna_het_type=d3a$type[match(data1$key3,d3a$key)]
data1$gene_name[which(is.na(data1$gene_name))]="D-loop"
d3b=ddply(data1,.(gene_name,tissue,het_id,nhetalleles,mtdna_het_type),summarise,count=sum(i))
d3b$nhet=d3b$count
d3b$nhet[which(d3b$nhetalleles==2)]=d3b$count[which(d3b$nhetalleles==2)]/2
d3b$nhet[which(d3b$nhetalleles==3)]=d3b$count[which(d3b$nhetalleles==3)]/3
d3b$i=1
d3e=ddply(d3b,.(gene_name,mtdna_het_type),summarise,count=sum(nhet),nhet=length(unique(het_id)))
d3f=ddply(d3e,.(gene_name),summarise,count=sum(nhet))
d3f=d3f[order(d3f$count,decreasing=T),]
d3e$gene_name=factor(d3e$gene_name,levels=rev(d3f$gene_name))
cols=c("#deebf7","#084594")
pdf("mthets_pergenetypes.pdf",height=5, width=4,useDingbats=FALSE)
p=ggplot(d3e, aes(gene_name,nhet,fill=mtdna_het_type)) + geom_col(colour="black",position="stack") + theme_bw() +ylab("Count")+xlab("mtDNA gene")
p=p+scale_fill_manual(values=cols) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## transition traversion ------------------------------------------------------------------------------
data1$hetchange = paste0(data1$inherited_allele,data1$heteroplasmic_base)
data1$titv = "Transversion"
data1$titv[which(data1$hetchange%in%c("AG","GA","CT","TC"))] = "Transition"
table(data1$titv)
#   Transition Transversion 
#        69164         9321 
# > 69164/(69164+9321)
# [1] 0.8812385
data1$strand="H-strand"
data1$strand[which(data1$gene_name%in%c("MT-ND6","MT-TQ","MT-TA","MT-TN","MT-TC","MT-TS2"))]="L-strand"
data1$strand[which(data1$gene_name%in%c("D-loop"))]="Intergenic"
# > table(data1$hetchange[which(data1$mtdna_het_type=="Somatic")],data1$strand[which(data1$mtdna_het_type=="Somatic")])
    
#      H-strand Intergenic L-strand
#   AC      377        104        0
#   AG     5996       2382      125
#   AT      139        105        5
#   CA      552        275        0
#   CG      251         74        0
#   CT     1361       3130      155
#   GA    10100       2025      141
#   GC      120          5        0
#   GT     1606         35        0
#   TA     1336        542        7
#   TC     2989       3430      202
#   TG       47        349        0
# > table(data1$strand[which(data1$mtdna_het_type=="Somatic")])

#   H-strand Intergenic   L-strand 
#      24874      12456        635 
# > 10100/24874
# [1] 0.4060465

a=data.frame(table(data1$hetchange[which(data1$mtdna_het_type=="Somatic")],data1$strand[which(data1$mtdna_het_type=="Somatic")]))
a$titv = "Transversion"
a$titv[which(a$Var1%in%c("AG","GA","CT","TC"))] = "Transition"
a=a[order(a$Var2,a$Freq,decreasing=T),]
a$Var1=factor(a$Var1,levels=a$Var1[which(a$Var2=="H-strand")])
a$Var2=factor(a$Var2,levels=c("H-strand","L-strand","Intergenic"))
a$titv=factor(a$titv,levels=c("Transition","Transversion"))
titvcols=c("#fa9fb5","#ae017e")

pdf("mthets_somatictitv.pdf",height=7, width=4,useDingbats=FALSE)
p=ggplot(a, aes(Var1,Freq,fill=titv)) + geom_col(colour="black") + theme_bw() +xlab("Somatic mtDNA heteroplasmy allele change")
p=p+scale_fill_manual(values=titvcols) + facet_wrap(vars(Var2),nrow=3,scales="free_y") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## get correlation between somatic mutations with mt-CN -----------------------------------------------
mtcn=read.csv("~/metadata/annotation/rath_pnas_mtcn.csv",as.is=T)
d3e=ddply(d3b,.(tissue,mtdna_het_type),summarise,count=sum(nhet),nhet=length(unique(het_id)))
d3e$tissue2=data1$tissue2[match(d3e$tissue,data1$tissue)]
d3e$mtcn_median=mtcn$Median.mtCN[match(d3e$tissue2,mtcn$tissue2)]
d3e1=d3e[which(d3e$mtdna_het_type=="Inherited"),]
d3e2=d3e[which(d3e$mtdna_het_type=="Somatic"),]
# > cor.test(d3e1$nhet,d3e1$mtcn_median)
# 	Pearson's product-moment correlation
# data:  d3e1$nhet and d3e1$mtcn_median
# t = -2.977, df = 46, p-value = 0.004631
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6157374 -0.1329699
# sample estimates:
#        cor 
# -0.4019195 
# > cor.test(d3e2$nhet,d3e2$mtcn_median)
# 	Pearson's product-moment correlation
# data:  d3e2$nhet and d3e2$mtcn_median
# t = -1.8762, df = 46, p-value = 0.06698
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.51196849  0.01895222
# sample estimates:
#        cor 
# -0.2666182 
d3e$mtdna_het_type=factor(d3e$mtdna_het_type,levels=c("Inherited","Somatic"))
cols=c("#9ecae1","#084594")
pdf("mthets_mtcn.pdf",height=4, width=5,useDingbats=FALSE)
p=ggplot(d3e, aes(mtcn_median,nhet,colour=mtdna_het_type)) + geom_point() + geom_smooth(method="lm",se=FALSE) 
p=p+scale_colour_manual(values=cols)+ theme_bw() +xlab("Median mtDNA-CN (Rath SP et al)")+ylab("N mtDNA heteroplasmic sites")
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

## examine tissue distribution of one pathogenic het site -----------------------------------------------
forplot=data[which(data$Pos==3243),]
forplot$annot=factor(forplot$annot)
annotcols=annot$colours[match(levels(forplot$annot),annot$annotation)]
pdf(paste0("mt3243.pdf"),height=4,width=5,useDingbats=FALSE)
p=ggplot(forplot,aes(tissue,sum_heteroplasmic_level,fill=annot))+geom_boxplot()+scale_fill_manual(values=annotcols)
p=p+theme_bw()+ylab("VAF")+xlab("Tissue")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()
