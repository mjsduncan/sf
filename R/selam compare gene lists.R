# Last time we didn't send you the moses gene list for rat(it took long time), here we attached yours and ours moses gene list for rat with combo programs.. and we found 243 genes in common.. and from your list 26 of them are found in Genage database and 50 of ours are in that database.. and we have 11 genes in common found in Genage. (ratMik_genes, ratUS_genes... refers yours and ours list of rat genes respectively)

length(intersect(ratMik_genes$feature ,ratUS_genes$feature))

#[1] 243

sort(intersect(ratMik_mapped$symbol, genage_human$aliases))

# [1] "ARNTL" "BDNF" "CDKN2B" "FEN1" "GHRHR" "IGF2" "IGFBP3" "INSR" "LRP2" "NFKBIA" "NGFR" "PDPK1" "PIK3R1" "PLCG2" "RAD52" "TCF3" "TGFB1"

# [18] "TNF" "TOP2A" "TP53" "TPP2" "UCHL1" "UCP2" "UCP3" "VEGFA" "XRCC5"

sort(intersect(ratUS_mapped$symbol, genage_human$aliases))

# [1] "ADCY5" "AGTR1" "AKT1" "ARNTL" "ATP5O" "BUB1B" "CACNA1A" "CDK1" "EEF1A1" "EGFR" "FLT1" "GPX1" "GRB2" "GSK3A" "GSK3B" "H2AFX"

# [17] "HIF1A" "HOXB7" "HOXC4" "HTT" "IGFBP3" "JAK2" "JUN" "KCNA3" "MAPK8" "MAPK9" "MDM2" "MLH1" "NFE2L2" "PCNA" "PLAU" "PLCG2"

# [33] "PTPN11" "RAD52" "RGN" "S100B" "SOD1" "SOD2" "STAT3" "TBP" "TCF3" "TNF" "TOP1" "TPP2" "UCP2" "UCP3" "VCP" "VEGFA"

# [49] "XRCC5" "YWHAZ"

intersect(intersect(ratUS_mapped$symbol, genage_human$aliases), intersect(ratMik_mapped$symbol, genage_human$aliases))

# [1] "UCP3" "PLCG2" "VEGFA" "TNF" "RAD52" "UCP2" "XRCC5" "TCF3" "IGFBP3" "TPP2" "ARNTL"

# So as you can see, we compare list of genes found in this case with our and yours list of moses genes for all organism. (humanEach_genes, mouseEach_genes and atEach_genes refers the binded gene moses list where each list found from moses fun on each dataset).
# 
# so found that we have a few genes in common(between binded moses list of individual run by organism, and yours & ours list of genes for all organism) . One gene from moses gene list for human found by binding individual moses run, is found in Genage and 2 from mouses and 1 gene from rat also found to be age related.



length(intersect(humanEach_genes$feature ,humanUS_genes$feature)) #[1] 1 "GPR20"

length(intersect(humanEach_genes$feature ,humanMik_genes$feature)) #[1] 1 "GCNT1"

intersect(humanEach_genes$feature, genage_human$aliases) # "GCLM"

##### compare mouses

length(intersect(mouseEach_genes$feature ,mouseUS_genes$feature)) #[1] 5 "Gcnt1" "Gcnt2" "Gcsam" "Gcsh" "Gdf10"

length(intersect(mouseEach_genes$feature ,mouseMik_genes$feature)) #[1] 3 "Gdf15" "Gdf5" "Gja3"

intersect(mouseEach_genes$feature, subset(genage_models, organism == "Mus musculus")$symbol)

#[1] "Ghr" "Ghrh"

intersect(mouseEach_mapped$symbol, genage_human$aliases)

#[1] "GHRH" "GHR"

#### compare rat

length(intersect(ratEach_genes$feature ,ratUS_genes$feature)) #[1] 3 "Fst" "Fstl1" "Papss2"

length(intersect(ratEach_genes$feature ,ratMik_genes$feature)) #[1] 4 "Frmd5" "Pafah1b1" "Parp2" "Pax4"

intersect(ratEach_mapped$symbol, genage_human$aliases)

# [1] "PARP1"