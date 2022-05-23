!#/usr/bin/env R


suppressPackageStartupMessages(library(tidyverse))

work_dir <- "~/projects/secretome/"
outcyte_dir <- paste0(work_dir,"Results_Outcyte/")

# load(paste0(work_dir,".RData")) #--> set1, set2, set3 for D.E-genes tables

myset = list.dirs(outcyte_dir)


###
### Read files with Concatenated results :
### (outcyte is run on multiple sequences for one gene-set --> results must be concatenated to one )


##> read all results of chunks of peptide-sequences
res <- read_delim( my_set_results, full.names = TRUE),
                  delim=" ", col_names = c("EntrezID","Prediction","Score"),
                  col_types = cols("EntrezID"="c","Prediction"="c","Score"="n"))


#####-----------Add GeneName to Outcyte Results --------------####


source("~/PROJECTS/useful_scripts/convertID2Symbol.R")

res$GeneName <- convert_ID2Symbol(as.character(res$EntrezID))


####------------------------  Add Prot-ID  ---------------------------####

.get_prot_id_and_descr <- function(ensembl, attributes_prot, input_df ){

    library(biomaRt)
    ensembl <- useMart("ensembl")
    ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl )
    attributes_prot <- c("entrezgene_id", "uniprotswissprot", "uniprotsptrembl")


    ##
    ## Group by EntrezID to have unique lines
    ##
    input_df <- input_df %>%
                group_by(EntrezID)  %>%
                summarise(Prediction = paste(unique(Prediction),collapse="//"),
                          Score=mean(Score),
                          GeneName = paste(unique(GeneName),collapse="//") )

    # res1 %>% head()
    # EntrezID Prediction       Score  GeneName
    # 10024    Intracellular    0.909  TROAP
    # 1006     Transmembrane    0.703  CDH8


    ## Extract swiss/embl proteinID for.each Entrez ID,
    ## add it in results-data-frame
    ##
    prot_id_df <- getBM(attributes = attributes_prot,
                     filters = "entrezgene_id",
                     values = as.character(input_df$EntrezID),
                     mart = ensembl )
    names(prot_id_df)[1] <- "EntrezID"

    ## Keep only one prot. per gene,group_by again :

    prot_id_df <- prot_id_df %>%
                  group_by(EntrezID) %>%
                  summarize(uniprot_embl = paste(uniprotsptrembl,collapse=","),
                            swiss_prot   = paste(uniprotswissprot,collapse=","))

    #> prot_id.res1 %>% head
    # entrezgene_id uniprot_embl
    #  41   //F8VSK4//H7C1W9//H0YHD6
    #  89   A0A087WSZ2////D6RH00
    #  215  //A6NEP8

    ##> merge both by matching the lines-EntrezIDs :
    mrg_df <- merge(input_df, prot_id_df ,by="EntrezID")


    genedesc <- getBM(attributes=c('external_gene_name','description'),
                      filters = 'external_gene_name',
                      values = mrg_df$GeneName,
                      mart = ensembl )
    names(genedesc)[1] <- "GeneName"

    if (any(duplicated(genedesc$GeneName) )) {
      genedesc <- genedesc[-duplicated(genedesc$GeneName),]
    }

    ##> merge again with Gene-Name :
    mrg2_df <- merge(mrg_df, genedesc, by="GeneName")

    return(mrg2_df)


}

res.mdf <- .get_prot_id_and_descr(ensembl,attributes_prot, res1 )


####-------------- Add Expression Values to Outcyte Results ------------####

res.mdf <- inner_join(res1.mdf, set1) %>% select(-c("AveExpr","t"))

###
### Save in EXCEL file
### (in an Excel Workbook,in diff.sheets)
###


library(readxl)
library(openxlsx)

write_in_excel <- function( wb, res_df, sheet_name ) {
    addWorksheet( wb=wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, res_df, colNames=TRUE)
}

wb <- createWorkbook(title="Outcyte_results_and_geneInfo")

# map2( .x=list(res1.mdf, res2.mdf, res3.mdf),
#       .y=list("DMD_vs_SETDB1", "WT-TGFb_vs_SETDB1", "WT-SETDB1_vs_DMD-SETDB1"),
#       ~write_in_excel(wb, .x,.y ))


write_in_excel(wb, res.mdf, "my_gene_set")
#Save Workbook
saveWorkbook(wb, paste0(outcyte_dir,"All_Outcyte_Results_and_Gene_Info.xlsx", overwrite = TRUE)
