
cd ~/new_scratch/FinnGen


module load gcc-5.4.0-gcc-4.8.5-fis24gg
module load r-3.6.1-gcc-5.4.0-zrytncq


library(data.table)
library(qqman)
d=fread("finngen_processed_single.txt",data.table=F)

i="FinnGen"

d = d[,c("chr","pos","chr_pos","pval")]
names(d) = c("CHR","BP","SNP","P")
d=d[which(d$CHR %in%c(1:22)),]
d$CHR=as.numeric(as.character(d$CHR))
d$BP=as.numeric(as.character(d$BP))
d$P=as.numeric(as.character(d$P))

observed <- sort(d$P)
lobs <- -(log10(observed))
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
obs.chisq = qchisq(1-d$P,1)
lambda = round(median(obs.chisq)/qchisq(0.5,1),3)

png(paste(i,"_QQ.png",sep=""),width=6,height=6, units="in",res=300)
par(omi=c(0,0.5,0,0))
plot(lexp,lobs,main=paste("Lambda",lambda,sep="="),xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",cex.axis=0.8)
abline(a=0,b=1)
dev.off()

png(paste(i,"_Manhattan.png",sep=""),units="in",width=12,height=6,res=300)
par(omi=c(0,0.5,0,0))
manhattan(d,suggestiveline = F,cex=0.6,cex.axis=0.8,genomewideline=-log10(5e-8))
dev.off()




cd ~/new_scratch/META


more 100KGP_saige4_send.txt | awk 'NR>1'|  cut -d " " -f 1 | sort | uniq -c | awk '$1==1{print $2}' > 100KGP_single.txt
more SAIGE_out_b38_ready_single_select.txt | cut -d " " -f 1 > BRIDGE_markers.txt
cat 100KGP_single.txt BRIDGE_markers.txt | sort | uniq -c | awk '$1==2 {print $2}' > UK_shared_markers.txt
more finngen_processed_single.txt | cut -d " " -f 1 > FinnGen_markers.txt
cat 100KGP_single.txt BRIDGE_markers.txt FinnGen_markers.txt | sort | uniq -c | awk '$1==3 {print $2}' > shared_markers.txt

~/software/METAL/generic-metal/metal UKmeta.txt
~/software/METAL/generic-metal/metal WithFinnMeta.txt
~/software/METAL/generic-metal/metal WithFinnMeta_verbose.txt > WithFinnMeta_log.txt

echo -e "# MARKER\tAL1\tAL2\tEFF\tSTDERR\tPVAL\tFREQ\tSOURCE" > FinnMeta_FinnData.txt
grep finngen_processed_single.txt WithFinnMeta_log.txt | grep -v Processing  >> FinnMeta_FinnData.txt
more FinnMeta_FinnData.txt | cut  -d " " -f 2- > FinnMeta_FinnData_noPound.txt

module load gcc-5.4.0-gcc-4.8.5-fis24gg
module load r-3.6.1-gcc-5.4.0-zrytncq



##Combining UK only

library(data.table)
library(qqman)


b=fread("SAIGE_out_b38_ready_single_select.txt",data.table=F)
g=fread("100KGP_saige4_send.txt",data.table=F)

x=fread("UK_shared_markers.txt",data.table=F)[,1]

names(g)=paste("G",names(g),sep="-")
names(b)=paste("B",names(b),sep="-")

g=g[which(g[,1] %in% x),]
b=b[which(b[,1] %in% x),]

m1 = merge(b,g,by.x="B-CHR_POS",by.y="G-CHR_POS")

a1=m1$'B-Allele1'
a2=m1$'B-Allele2'
a3=m1$'G-Allele1'
a4=m1$'G-Allele2'

matched = which( ( a1==a3 & a2==a4) | (a1==a4 & a2==a3))
m1_good = m1[matched,]

write.table(m1_good,"UK_Combined_data.txt",quote=F,row.names=F,col.names=T)

f=fread("UKmeta_1.tbl",data.table=F)
x=merge(f,m1_good,by.x="MarkerName",by.y="B-CHR_POS")

bad = which(x$MaxFreq > x$MinFreq + 0.05 | x$HetPVal<1e-5)

png("UKmeta_MinMaxFreq.png")
plot(x$MinFreq,x$MaxFreq)
points(x$MinFreq[bad],x$MaxFreq[bad],col=2)
dev.off()

x2=x[-bad,]
write.table(x2,"UK_Combined_data_good.txt",quote=F,row.names=F,col.names=T)
write.table(x2[,1],"UK_Combined_data_good_markers.txt",quote=F,row.names=F,col.names=F)

grow=c()
for (i in 1:22){
	print(i)
	z=fread(paste("../dbSNP/b38_common_chr",i,".txt",sep=""),data.table=F)
	z$ID = paste(z[,1],z[,2],sep="_")
	zu=z[which(z$ID %in% names(which(table(z$ID)==1))),]
	keep = zu[which(zu$ID %in% x2$MarkerName),]
	grow=rbind(grow,keep)
}
names(grow)=c("CHR","POS","rsID","REF","ALT","ID")

m3 = merge(x2,grow,by.x="MarkerName",by.y="ID",all.x=T)

write.table(m3,"UK_Combined_data_good_rsIDs.txt",quote=F,row.names=F,col.names=T)


##PLOT UK only results!

#Read in first
x=fread("UK_Combined_data_good_rsIDs.txt",data.table=F)
m3=x

i="UK_META"

d = m3[,c("B-CHR","B-POS","MarkerName","P-value")]
names(d) = c("CHR","BP","SNP","P")

observed <- sort(d$P)
lobs <- -(log10(observed))
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
obs.chisq = qchisq(1-d$P,1)
lambda = round(median(obs.chisq)/qchisq(0.5,1),3)

png(paste(i,"_QQ.png",sep=""),width=6,height=6, units="in",res=300)
par(omi=c(0,0.5,0,0))
plot(lexp,lobs,main=paste("Lambda",lambda,sep="="),xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",cex.axis=0.8)
abline(a=0,b=1)
dev.off()

png(paste(i,"_Manhattan.png",sep=""),units="in",width=12,height=6,res=300)
par(omi=c(0,0.5,0,0))
manhattan(d,suggestiveline = F,cex=0.6,cex.axis=0.8,genomewideline=-log10(5e-8))
dev.off()



##Combining with Finnish

#library(data.table)
#x=fread("UK_Combined_data_good_rsIDs.txt",data.table=F)
x=m3

#f_input=fread("finngen_processed_single.txt",data.table=F) #use output from META log instead as the AF/EFFECT seem to be correct with the UK data
f_input=fread("FinnMeta_FinnData_noPound.txt",data.table=F)
f_meta=fread("WithFinnMeta_1.tbl",data.table=F)
f=merge(f_input,f_meta,by.x="MARKER",by.y="MarkerName")
names(f)=paste("F",names(f),sep="_")
z=fread("shared_markers.txt",data.table=F)
f2 = f[which(f$F_MARKER %in% z[,1]),]
x2 = x[which(x$MarkerName %in% z[,1]),]
m1 = merge(x2,f,by.x="MarkerName",by.y="F_MARKER")

a1=m1$'Allele1'
a2=m1$'Allele2'
a3=m1$'F_Allele1'
a4=m1$'F_Allele2'

matched = which( ( a1==a3 & a2==a4) | (a1==a4 & a2==a3))
m1_good = m1[matched,]

bad = which(m1_good$F_MaxFreq > m1_good$F_MinFreq + 0.25 | m1_good$F_HetPVal<1e-5)

png("WithFinMeta_MinMaxFreq.png")
plot(m1_good$F_MinFreq,m1_good$F_MaxFreq)
points(m1_good$F_MinFreq[bad],m1_good$F_MaxFreq[bad],col=2)
dev.off()

m2=m1_good[-bad,]

write.table(m2,"UKFinn_Combined_data_good.txt",quote=F,row.names=F,col.names=T)

b37=fread("b37_b38_ready.txt",data.table=F)
b37$MarkerName=paste(b37$CHR38,b37$POS38,sep="_")
b37$CHR_POS_37=paste(b37$CHR37,b37$POS37,sep="_")

m3=merge(m2,b37)

keep=c("B-CHR",
"B-POS",
"MarkerName",
"POS37",
"CHR_POS_37",
"REF37",
"ALT37",
"rsID",
"Allele1",
"Allele2",
"Freq1",
"Effect",
"StdErr",
"P-value",
"F_AL1",
"F_AL2",
"F_FREQ",
"F_EFF",
"F_STDERR",
"F_PVAL",
"F_Allele1",
"F_Allele2",
"F_Freq1",
"F_Effect",
"F_StdErr",
"F_P-value")

m4=m3[,keep]

names(m4)=c("CHR",
"B38",
"CHR_B38",
"B37",
"CHR_B37",
"REF",
"ALT",
"rsID",
"UK_A1",
"UK_A2",
"UK_Freq",
"UK_Beta",
"UK_StdErr",
"UK_P",
"F_A1",
"F_A2",
"F_Freq",
"F_Beta",
"F_StdErr",
"F_P-value",
"Meta_A1",
"Meta_A2",
"Meta_Freq",
"Meta_Beta",
"Meta_StdErr",
"Meta_P")

write.table(m4,"UKFinn_Combined_data_good_select_b38b37.txt",quote=F,row.names=F,col.names=T)

#Read in
x=fread("UKFinn_Combined_data_good_select_b38b37.txt",data.table=F)
m4=x

i="FINAL"
d = m4[,c("CHR","B38","CHR_B38","Meta_P")]
names(d) = c("CHR","BP","SNP","P")

observed <- sort(d$P)
lobs <- -(log10(observed))
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
obs.chisq = qchisq(1-d$P,1)
lambda = round(median(obs.chisq)/qchisq(0.5,1),3)

png(paste(i,"_QQ.png",sep=""),width=6,height=6, units="in",res=300)
par(omi=c(0,0.5,0,0))
plot(lexp,lobs,main=paste("Lambda",lambda,sep="="),xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",cex.axis=0.8)
abline(a=0,b=1)
dev.off()

png(paste(i,"_Manhattan.png",sep=""),units="in",width=12,height=6,res=300)
par(omi=c(0,0.5,0,0))
manhattan(d,suggestiveline = F,cex=0.6,cex.axis=0.8,genomewideline=-log10(5e-8))
dev.off()


#READ IN to export top loci

library(data.table)
x=fread("UKFinn_Combined_data_good_select_b38b37.txt",data.table=F)
d=as.character(read.table("b38_top_loci.txt",header=F)[,1])


z = x[which(x$CHR_B38 %in% d),]
grow=c()
for (i in 1:length(z[,1])){
	chr=z[i,"CHR"]
	pos=z[i,"B38"]
	min = pos - 1e6
	max = pos + 1e6
	found=x[which(x$CHR==chr & x$B38<max & x$B38>min),]
	grow=rbind(grow,found)
}

grow2 = grow[order(grow$CHR,grow$B38),]
write.table(grow2,"Top_variants_1MB.txt",row.names=F,col.names=T,sep="\t",quote=F)

