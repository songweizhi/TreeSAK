library(Biostrings)
options(digits=8)
protein <- readAAStringSet("/lomi_home/wenxiu/RIS_virus/workdir/27_ortho_RIStree/ortho_tree/v1/tree_workdir/test.aln")
outgroup <- c()
protein <- protein[setdiff(names(protein),outgroup)]
library(stringr)
chi2_2 <- function(set,taxa_number,all_amino,amino_number){
  O <- str_count(set[[taxa_number]],all_amino[amino_number])##某个物种中(已经除去一个位置的氨基酸之后)的全部氨基酸
  E <- sum(str_count(set,all_amino[amino_number]))/sum(nchar(str_replace_all(set,"-","")))*nchar(str_replace_all(set[[taxa_number]],"-",""))###
  return((O-E)^2/E)
}

chi2_1 <- function(set,taxa_number){
  all_amino <- unique(strsplit(as.character(str_replace_all(set,"-","")),"")[[1]])
  return(sum(sapply(1:length(all_amino),chi2_2,taxa_number=taxa_number,set=set,all_amino=all_amino)))
}

untrim <- sum(sapply(1:length(protein),chi2_1,set=protein))

#13746.9

trimmed_chi2 <- function(x){
  library(Biostrings)
  protein <- readAAStringSet("/lomi_home/wenxiu/RIS_virus/workdir/27_ortho_RIStree/ortho_tree/v1/tree_workdir/test.aln")
  outgroup <- c()
  protein <- protein[setdiff(names(protein),outgroup)]
  library(stringr)
  chi2_2 <- function(set,taxa_number,all_amino,amino_number){
    O <- str_count(set[[taxa_number]],all_amino[amino_number])#某一个氨基酸(传进来的)在这个物种出现的频次
    E <- sum(str_count(set,all_amino[amino_number]))/sum(nchar(str_replace_all(set,"-","")))*nchar(str_replace_all(set[[taxa_number]],"-",""))##氨基酸出现在配对中出现的的总次数/总氨基酸数*这个序列的氨基酸数
    return((O-E)^2/E)
  }
  
  chi2_1 <- function(set,taxa_number){###
    all_amino <- unique(strsplit(as.character(str_replace_all(set,"-","")),"")[[1]])##除去补空位的,在该位置的所有氨基酸
    return(sum(sapply(1:length(all_amino),chi2_2,taxa_number=taxa_number,set=set,all_amino=all_amino)))###对于其中的每一种氨基酸
  }
  
  end_pos_1 <- x-1
  start_pos_2 <- x+1
  if (end_pos_1>=1 & start_pos_2<=width(protein)[1]){
    trimmed_protein <- str_c(substr(protein,1,end_pos_1),substr(protein,start_pos_2,width(protein)[1]))###str_c元素对元素的合并列表
    print(trimmed_protein)
  }else if(start_pos_2<=width(protein)[1]){
    trimmed_protein <- substr(protein,start_pos_2,width(protein)[1])#从start到全长
  }else {
    trimmed_protein <- substr(protein,1,end_pos_1)#从起始到end
  }#到此截取了想要的氨基酸（所有序列的）
  return(sum(sapply(1:length(trimmed_protein),chi2_1,set=trimmed_protein)))###对于每一个物种
}


library(parallel)

cl<-makeCluster(48)

trimmed_chi2_set <- parLapply(cl,1:width(protein)[1],trimmed_chi2)

stopCluster(cl)

save.image("/lomi_home/wenxiu/RIS_virus/workdir/27_ortho_RIStree/ortho_tree/v1/tree_workdir/test_removed_bias.aln")


###删除对应的位点
trim_site<-function(percentage){
  protein1<-protein
  realnum=c()
  for (i in 1:length(trimmed_chi2_set)) {
    realnum= c(realnum,trimmed_chi2_set[[i]])
  }
  #realnum=abs(realnum-untrim)
  print(percentage*length(realnum))
  kafang<-realnum[order(realnum,decreasing = F)[round(percentage*length(realnum))]]
  site_to_be_removed=which(realnum<=kafang)
  site<-IRanges(start=site_to_be_removed, end=site_to_be_removed, width=1)
  protein
  replaceAt(protein1,site,'')
  protein1<-replaceAt(protein1,site,'')
  substr(protein1[1],start=1,stop=width(protein1)[1])
  #10：5984
  #print(paste0('~/removed_bias_dataset3_concatenate_',percentage,'.txt'))
  writeXStringSet(protein1,paste0('/lomi_home/wenxiu/RIS_virus/workdir/27_ortho_RIStree/ortho_tree/v1/tree_workdir/test_removed_bias.trimmed.concat_',percentage*100,'%.fas'))
}

  
for (cutoff in c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,0.9)) {
  trim_site(cutoff)
}








###先分成两份，再取5%
# protein1<-protein
# realnum=c()
# for (i in 1:length(trimmed_chi2_set)) {
#   realnum= c(realnum,trimmed_chi2_set[[i]])
# }
# realnumplus<-realnum[which(realnum-untrim>0)]
# realnummiuis<-realnum[which(realnum-untrim<0)]
# plus<-realnumplus[order(realnumplus,decreasing = T)[1:(0.1*length(realnumplus)+1)]]
# miuis<-realnummiuis[order(realnummiuis,decreasing = F)[1:(0.1*length(realnummiuis)+1)]]
# site_to_be_removed=which(realnum %in% c(plus,miuis))
# site<-IRanges(start=site_to_be_removed, end=site_to_be_removed, width=1)
# protein
# replaceAt(protein1,site,'')
# protein1<-replaceAt(protein1,site,'')
# substr(protein1[1],start=1,stop=5984)
