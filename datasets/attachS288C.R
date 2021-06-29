# attach s288C protein seq at the very beginning of each protein seq file

library(seqinr)

setwd("E:/SeqMoid/datasets")

file.list<-dir('E:/SeqMoid/datasets/allrefgene_translated')

s288c<-read.csv('S288C_ProteinSequence.csv',header=T)

# for each gene in the translated folder, append s288c aas seq

for (f in file.list) {
  
  trans.file<-read.fasta(paste0('allrefgene_translated/',f),seqtype = 'AA')
  
  gene.name<-gsub('.fasta','',f)
  
  s288c.seq<-s288c$ProteinSeq[which(s288c$SystematicName==gene.name)]
  
  annot.raw<-as.character(getAnnot(trans.file[1]))
  
  annot.raw<-strsplit(annot.raw,'\t')
  
  strain.name<-annot.raw[[1]][1]
  
  strain.name<-strsplit(strain.name,'_')
  
  s.len<-length(strain.name[[1]])
  
  s288c.name<-paste('S288C', strain.name[[1]][s.len-1],strain.name[[1]][s.len],sep='_')
  
  annot.s288c<-paste(s288c.name,annot.raw[[1]][2],annot.raw[[1]][3],sep='\t')
  
  write.fasta(sequences = s288c.seq, names=annot.s288c,
              file.out=paste0('lab_1011_translated/',f))
  
  trans.seq<-getSequence(trans.file,as.string=T)
  
  trans.annot<-getAnnot(trans.file)
  
  trans.annot<-gsub('>','',trans.annot)
  
  write.fasta(sequences = trans.seq,names=trans.annot,
              file.out=paste0('lab_1011_translated/',f),open='a')
  
}
