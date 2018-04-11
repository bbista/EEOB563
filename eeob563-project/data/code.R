Z<- read.table("../Zlinked/file.txt")
a=""

for (i in Z$V1){
  fas<-as.character(paste0(i))
  fas1<-as.character(paste0("../Zlinked/",fas))
  system (paste0("prank -d=",fas1, " -o=output -codon -once"))
  template<-readLines(paste0("output.best.fas"))
  writeLines(template,paste0("Alignment/prank_",fas))
  
  writeLines(template,"pamlqueryseq.fas")
  system("codeml")
  template1<-readLines(paste0("pamloutput.txt"))
  writeLines(template1,paste0("Output/",substr(fas,1,(nchar(fas)-6)),".txt"))
  writeLines(a,paste0("output.best.fas"))
}