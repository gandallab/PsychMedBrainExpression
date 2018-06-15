# 01_makeDatExpr.R

options(stringsAsFactors = FALSE)
#source("http://bioconductor.org/biocLite.R")
#biocLite("oligo")

setwd("/u/project/gandalm/sparhami/PsychotropicTx")

library(affy)
library(GEOquery)
library(oligo)


# function to extract metadata and expression for raw Affymetrix sets
ExtractAffy <- function(accession,experiment,oligo=NULL){
  gse <- getGEO(accession, GSEMatrix=TRUE,getGPL=FALSE)
  datMeta <- pData(gse[[1]])
  rownames(datMeta) <- datMeta[,2]
  colnames(datMeta) <- gsub(":ch1","",colnames(datMeta))
  colnames(datMeta) <- gsub("_ch1.*","",colnames(datMeta))
  
  celpath <- paste0("raw_data/",experiment)
  listdir <- list.files(celpath)
  celfiles <- listdir[grep("CEL",listdir,ignore.case=TRUE)]
  
  if (oligo){
    data.affy <- read.celfiles(paste0(celpath,"/",celfiles))
  }else{
    data.affy <- ReadAffy(celfile.path = celpath, filenames = celfiles)
  }
  
  datExpr <- exprs(data.affy)
  filemeta <- paste0(celpath,"/datMeta_",accession,".RData")
  save(datMeta,file = filemeta)
  fileexpr <- paste0(celpath,"/datExpr_",accession,".RData")
  save(datExpr,file = fileexpr)
}


ExtractAffy("GSE28644","Benton_Ms_SSRI_2012")
ExtractAffy("GSE2880","Hassel_RatFC_Anticonvulsant_2010")
ExtractAffy("GSE66276","Lanz_Rat_MoodStabilizer2015")
ExtractAffy("GSE43748","Noh_RatPFC_Stimulant_2013")
ExtractAffy("GSE33619","Sadasivan_MsSN_stimulant_2012")
ExtractAffy("GSE43261","Samuels_MsHC_SSRI_2013")
ExtractAffy("GSE28435","Spradling_Rat_Anticholinergic_2011")
ExtractAffy("GSE56780","Zhao_MsHC_Anticonvulsant_2017")

ExtractAffy("GSE76700","Furukawa_MsCortex_Benzo_2016",oligo=TRUE)
ExtractAffy("GSE93041","Orozco_Ms_Ketamine_2017",oligo=TRUE)
ExtractAffy("GSE56028","Patr├нcio_Rat_SSRI_2015",oligo=TRUE)
