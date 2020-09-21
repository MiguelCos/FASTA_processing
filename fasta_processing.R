### FASTA processing -----

## Load packages ----

library(tidyverse)
library(seqinr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UniprotR)
library(UniProt.ws)

## Load the fasta ----

fasta <- read.fasta(file = "EBI_Human_reference_1protein_1gene_RELEASE_2020_02_w_iRTs_plus_scrambled.fasta", seqtype = "AA")

## Checking ----

head(fasta)

head(names(fasta))

getAnnot(head(fasta))

getAnnot(fasta) %>% unlist() %>% str_detect(., "SCRAMBLED") %>% table(.) # there are scrambled sequences here

getName(fasta) %>% head()

tnames <- getName(fasta) 

## Get only the scrambled sequences ----

scrafast <- Filter(function(x) str_detect(getAnnot(x), "SCRAMBLED"),
                   fasta)

tannot <- getAnnot(scrafast) %>% unlist()

descram <- str_remove(tannot, ".*-- ")

reannot <- str_split(descram, "\\|", n = 3, simplify = TRUE)

reannot <- as.data.frame(reannot)

reannot <- mutate(reannot,
                  proteinID = paste0(V2,"_SCRAMBLED")) %>%
      mutate(annotation = paste0(V1,"|",proteinID,"|",V3,"_SCRAMBLED"))


length(scrafast)

names(scrafast) <- reannot$annotation

head(scrafast)

## Get non-scrambled sequences ----

goodfast <- Filter(function(x) str_detect(getAnnot(x), "SCRAMBLED", negate = TRUE),
                   fasta)

length(goodfast)
head(goodfast)

## Join both ----

fastanew <- c(goodfast, scrafast)

## Every sequence to capital letters ----

fastanew2 <- lapply(fastanew, str_to_upper)

tail(fastanew2)

## Save new fasta ----  

"EBI_Human_reference_1protein_1gene_RELEASE_2020_02_w_iRTs_plus_wSCRAMBLED_named2.fasta"

write.fasta(fastanew2, names = names(fastanew2), 
            file.out = "EBI_Human_reference_1protein_1gene_RELEASE_2020_02_w_iRTs_plus_wSCRAMBLED_named2.fasta")

## Filter-out "empty" elements of the fasta ----
# Some IDs don't contain any sequences, those should be eliminated

fastall <- Filter(function(x) length(x) > 0,
                  fasta)

head(fastall, 3)

names(fastall) <- getAnnot(fastall)

## Filter-out redundant sequences ----

sequences <- sapply(fastall,
                    function(x) paste0(x, collapse = ""))

logidup <- duplicated(sequences)

noduplifast <- fastall[!logidup]

table(str_detect(names(noduplifast), ">ENST00000623834 Ref"))

## Extract Ensembl ID ----

uniproting <- uniprot_get("9606")


ensemblIDs <- getName(noduplifast)



ensembl2uniprot2 <- UniprotR::ConvertID(ProteinAccList = ensemblIDs,
                                        ID_from = , ID_to = "ACC+ID")

ensembl2uniprot <- bitr(ensemblIDs, fromType = "ENSEMBLTRANS",
                        toType = c("UNIPROT","SYMBOL","UNIGENE"),
                        OrgDb = "org.Hs.eg.db", drop = FALSE)

## Subset only variant sequences -----

varfast <- Filter(function(x) str_detect(getAnnot(x), "Variant"),
                  noduplifast)
length(varfast)

write.fasta(varfast,names = str_remove(names(varfast), "^>"), 
            file.out = "GBM_specific_variant_proteins.fasta")

## Subset only Reference sequences ---- 

refast <- Filter(function(x) str_detect(getAnnot(x), "Ref$"),
                 noduplifast)

write.fasta(refast,names = str_remove(names(refast), "^>"), 
            file.out = "GBM_specific_reference_proteins.fasta")


head(refast)

length(refast)