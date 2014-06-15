### http://genomics.senescence.info data sets

# genge_human needed final column name added to downloaded csv file "homologue"
genage_human <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_human.csv", stringsAsFactors = FALSE)
genage_models <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_models.csv", stringsAsFactors = FALSE)
gendr_manipulations <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/gendr_manipulations.csv", stringsAsFactors = FALSE)
longevity <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/longevity.csv", stringsAsFactors = FALSE)
anage_data <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/anage_data.txt", stringsAsFactors = FALSE)

# rm(genage_human, genage_models, gendr_manipulations, longevity, anage_data)

### tcm data base http://www.megabionet.org/tcmid/download/
compound_protein <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/tcmid/compound_protein.csv", header = FALSE, stringsAsFactors = FALSE)
names(compound_protein) <- c("compound", "uniprot", "source")
compounds <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/tcmid/compounds.csv", stringsAsFactors = FALSE)
herb <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/tcmid/herb.csv", stringsAsFactors = FALSE)
prescription <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/tcmid/prescription.csv", stringsAsFactors = FALSE)

# rm(compounds, herb, prescription, compound_protein)

### set up postgresql db for comparative tox db http://ctdbase.org/
library("RPostgreSQL")
m <- dbDriver("PostgreSQL")
con <- dbConnect(m, user="coguser", password="cog", dbname="ctdb", port = 5433)

chemicals <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/chemicals.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "chemicals", chemicals)

diseases <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/diseases.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "diseases", diseases)

genes <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/genes.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "genes", genes)

pathways <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/pathways.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "pathways", pathways)

genes_pathways <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/genes_pathways.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "genes_pathways", genes_pathways)

diseases_pathways <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/diseases_pathways.csv", stringsAsFactors = FALSE)
dbWriteTable(con, "diseases_pathways", diseases_pathways)

chemicals_diseases <- read.clump("C:/Users/user/Desktop/biomind/artificial biologist/stevia/ctdb/CTD_chemicals_diseases.csv", lines = 35, clump = 1)
names(chemicals_diseases) <- c("ChemicalName", "ChemicalID", "CasRN", "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubMedIDs")
dbWriteTable(con, "chemicals_diseases", chemicals_diseases)

# rm(chemicals, diseases, genes, pathways, genes_pathways, diseases_pathways, chemicals_diseases)

##### functions

read.clump <- function(file, lines, clump, readFunc=read.csv,
    skip=(lines*(clump-1))+ifelse((header) & (clump>1) & (!inherits(file, "connection")),1,0),
    nrows=lines,header=TRUE,...){
    if(clump > 1){
            colnms<-NULL
            if(header)
            {
                colnms<-unlist(readFunc(file, nrows=1, header=FALSE))
                print(colnms)
            }
      p = readFunc(file, skip = skip,
          nrows = nrows, header=FALSE,...)
            if(! is.null(colnms))
            {
        colnames(p) = colnms
            }
    } else {
        p = readFunc(file, skip = skip, nrows = nrows, header=header)
    }
    return(p)
}