library(plyr)
library(dplyr)
library(tidyr)


library(msa)
library(ape)
library(ips)
library(seqinr)
library("Biostrings")

auto<-read.table("auto/file.txt", stringsAsFactors = F)
z<-read.table("z/file.txt", stringsAsFactors = F )


i<-auto$V1[1]
for (i in auto$V1){
  fastaFile <- readDNAStringSet(paste0("auto/",i), format="FASTA")
  sequence = paste(fastaFile)
  seqname=names(fastaFile)
  seq<-data.frame(seqname, sequence, stringsAsFactors = F)
  fa <- c(paste0(">", seq$seqname[1]), as.character(paste0(seq$sequence[1])),sprintf(">%s", seq$seqname[2]),as.character(paste0(seq$sequence[2])))
  writeLines(fa, paste0("pamlqueryseq.fas")) 
  system("codeml")
  template=readLines("pamloutput.txt")
  writeLines(template, paste0("aGGA/", i))
  
  fa <- c(paste0(">", seq$seqname[3]), as.character(paste0(seq$sequence[3])),sprintf(">%s", seq$seqname[4]),as.character(paste0(seq$sequence[4])))
  writeLines(fa, paste0("pamlqueryseq.fas")) 
  system("codeml")
  template=readLines("pamloutput.txt")
  writeLines(template, paste0("aASI/", i))
}


for (i in z$V1){
  fastaFile <- readDNAStringSet(paste0("z/",i), format="FASTA")
  sequence = paste(fastaFile)
  seqname=names(fastaFile)
  seq<-data.frame(seqname, sequence, stringsAsFactors = F)
  fa <- c(paste0(">", seq$seqname[1]), as.character(paste0(seq$sequence[1])),sprintf(">%s", seq$seqname[2]),as.character(paste0(seq$sequence[2])))
  writeLines(fa, paste0("pamlqueryseq.fas")) 
  system("codeml")
  template=readLines("pamloutput.txt")
  writeLines(template, paste0("zGGA/", i))
  
  fa <- c(paste0(">", seq$seqname[3]), as.character(paste0(seq$sequence[3])),sprintf(">%s", seq$seqname[4]),as.character(paste0(seq$sequence[4])))
  writeLines(fa, paste0("pamlqueryseq.fas")) 
  system("codeml")
  template=readLines("pamloutput.txt")
  writeLines(template, paste0("zASI/", i))
}

#change directory




out<-read.table("zASI/file.txt", stringsAsFactors = F)

file= data.frame(V1=character(), stringsAsFactors = F)
for (i in out$V1){
    fas1<-paste0( as.character(i))
  output<-readLines(paste0("zASI/",fas1))
  n=length(output)
  x=output[n]
  if (x!="Codon usage in sequences") {
    file[nrow(file)+1,]<-c(i)
  }
  
}


Autodata<-data.frame(gene.ID=character(),t=numeric(),S=numeric(),N=numeric(),dnds=numeric(),dN=numeric(),ds=numeric(),lnl=numeric(), stringsAsFactors = F)
for(i in file$V1){
  
  fas1<-paste0("zASI/", i)
  output<-readLines(paste0(fas1))
  n<-length(output)
  
  txt<-output[n-2]
  wnum<-as.numeric(unlist(regmatches(txt,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",txt))))
  for (j in 1:n) { #looking for the line that has the log likelihood value
    subset<-substr(output[j],1,3)
    if (subset=="lnL") {
      l<-j
    }
  }
  ltxt<-output[l]
  like<-as.numeric(unlist(regmatches(ltxt,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",ltxt))))
  Autodata[nrow(Autodata)+1,]<-c(substr(i,7,nchar(i)-6),wnum,(-1*like))
}

write.csv(Autodata, "zmankASI.csv")


#permutation tests

autoset<-read.csv("automankGGA.csv", stringsAsFactors = F)
xset<-read.csv("zmankGGA.csv", stringsAsFactors = F)

xset<-xset[xset$ds<2,]

xset$Ndn=xset$dN*xset$N
xset$Sds=xset$ds*xset$S
xset$gene.ID<-NULL
zsum=colSums(xset)
zdn=zsum[9]/zsum[4]
zds=zsum[10]/zsum[3]
zdnds=zdn/zds

autoset<-autoset[autoset$ds<2,]

autoset$Ndn=autoset$dN*autoset$N
autoset$Sds=autoset$ds*autoset$S
autoset$gene.ID<-NULL
asum=colSums(autoset)

adn=asum[9]/asum[4]
ads=asum[10]/asum[3]
adnds=adn/ads
data<-rbind(autoset,xset)


p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(autoset)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  padnds<-as.numeric((aasum[9]/aasum[4])/(aasum[10]/aasum[3]))
  pzdnds<-as.numeric((zzsum[9]/zzsum[4])/(zzsum[10]/zzsum[3]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzdnds-padnds)>abs(zdnds-adnds)){
    p1=p1+1
  }
  if ((pzdnds-padnds)<(-abs(zdnds-adnds))){
    p2=p2+1
  }
}
p1
p2
pvaluednds<-(p1+p2)/10000
#
#
#
#
#
#
#dn
p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(autoset)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  padn<-as.numeric((aasum[9]/aasum[4]))
  pzdn<-as.numeric((zzsum[9]/zzsum[4]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzdn-padn)>abs(zdn-adn)){
    p1=p1+1
  }
  if ((pzdn-padn)<(-abs(zdn-adn))){
    p2=p2+1
  }
}
p1
p2
pvaluedn<-(p1+p2)/10000


#ds
p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(autoset)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  pads<-as.numeric((aasum[10]/aasum[5]))
  pzds<-as.numeric((zzsum[10]/zzsum[5]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzds-pads)>abs(zds-ads)){
    p1=p1+1
  }
  if ((pzds-pads)<(-abs(zds-ads))){
    p2=p2+1
  }
}
p1
p2
pvalueds<-(p1+p2)/10000



pvaluednds
pvaluedn
pvalueds


#bootstrap
aboot<-data.frame(dnds=as.numeric(),dn=as.numeric(), ds=as.numeric(), stringsAsFactors = F)
zboot<-data.frame(dnds=as.numeric(),dn=as.numeric(), ds=as.numeric(),  stringsAsFactors = F)


aa<-autoset
zz<-xset

aa$taxa<-NULL
aa$gene.ID<-NULL

zz$taxa<-NULL
zz$gene.ID<-NULL

for (i in 1:10000){
  n<-sample(1:nrow(aa),replace = T)
  ap<-as.matrix.data.frame(aa[c(n),])
  n<-sample(1:nrow(zz),replace = T)
  zp<-as.matrix.data.frame(zz[c(n),])
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  badnds<-c(dnds=as.numeric((aasum[9]/aasum[4])/(aasum[10]/aasum[3])),dn=as.numeric((aasum[9]/aasum[4])),ds=as.numeric((aasum[10]/aasum[3])))
  bzdnds<-c(dnds=as.numeric((zzsum[9]/zzsum[4])/(zzsum[10]/zzsum[3])),dn=as.numeric((zzsum[9]/zzsum[4])),ds=as.numeric((zzsum[10]/zzsum[3])))
  aboot[nrow(aboot)+1,]<-c(badnds)
  zboot[nrow(zboot)+1,]<-c(bzdnds)
  
}
a<-aboot$dnds[order(aboot$dnds)]
a[250]
a[10000-250]

a<-aboot$dn[order(aboot$dn)]
a[250]
a[10000-250]

a<-aboot$ds[order(aboot$ds)]
a[250]
a[10000-250]


z<-zboot$dnds[order(zboot$dnds)]
z[250]
z[10000-250]

z<-zboot$dn[order(zboot$dn)]
z[250]
z[10000-250]

z<-zboot$ds[order(zboot$ds)]
z[250]
z[10000-250]

