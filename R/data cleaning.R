### import datasets

# import geo data sets
library("GEOquery")
library("GEOmetadb")
# from de Magalhães 2009 supplement
human <- c(
  h_brain = "GDS707", # Lu, et al., 2004
  muscle1 = "GDS472", muscle2 = "GDS473", # Welle, et al., 2004
  muscle3 = "GDS287", muscle4 = "GDS288", # Welle, et al.,2003
  muscle5 = "GDS156" # Welle, et al., 2002
)
mouse <- c(
  muscle = "GDS2612", # Edwards, et al., 2007
  kidney1 = "GDS355", kidney2 = "GDS356",
  m_brain = "GDS1311", # Godbout, et al., 2005
  m_hippo = "GDS2082", # Verbitsky, et al., 2004
  liver = "GDS2019", # Papaconstantinou, et al., 2005
  m_heart = "GDS40",
  lung = "GDS2929", # Misra, et al., 2007
# not in GEO ftp site #  eye1 = "GDS396", eye2 = "GDS397", #Ida, et al., 2003
  cochlea = "GDS2681", # Someya, et al., 2007
  hemato_stem = "GDS1803", # Rossi, et al., 2005
  myo_progen = "GDS1079" # Beggs, et al., 2004
)
rat <- c(
  r_hippo = "GDS2639", # Rowe, et al., 2007
  stromal = "GDS2231", # Akavia, et al., 2006
  spinal_cord = "GDS1280",
  oculomotor = "GDS1280",
  skeletal_ms = "GDS1279",
  extraoc_ms = "GDS1279",
  laryngeal_ms = "GDS1278", # McMullen and Andrade, 2006
  r_heart = "GDS399", # Dobson, et al., 2003
# not in GEO ftp site # CA1_hipp1 = "GDS2315", # Burger, et al., 2007
  CA1_hipp2 = "GDS520" # Blalock, et al., 2003
)
geodsets <- c(human, mouse, rat)

# get annotation packages for geo datasets
getSQLiteFile(destdir = 'C:/Users/user/Documents')
gpls <- geoConvert(geodsets, out_type = "gpl", sqlite = 'C:/Users/user/Documents/GEOmetadb.sqlite')[[1]]
row.names(gpls) <- names(geodsets[match(gpls[,1], geodsets)])
con <- dbConnect(SQLite(),'C:/Users/user/Documents/GEOmetadb.sqlite')
arrays <- merge(gpls, getBiocPlatformMap(con, bioc='all'), by.x = "to_acc", by.y = "gpl", all.x = TRUE)
dbDisconnect(con)
# file.remove('C:/Users/user/Documents/GEOmetadb.sqlite')
BiocInstaller::biocLite(paste(unique(arrays[["bioc_package"]]), "db", sep = "."))

# dowlnoad geo datasets into list of GEOData objects
gdsData <- vector("list",length(geodsets))
names(gdsData) <- names(geodsets)
for(i in 1:length(geodsets)){gdsData[[i]] <- getGEO(GEO = geodsets[i], destdir = '~/github/stevia/data')}

# make list of expression set data frames
gdsExp <- vector("list",length(geodsets))
names(gdsExp) <- names(geodsets)
for(i in 1:length(geodsets)){gdsExp[[i]] <- Table(gdsData[[i]])}
gdsExp <- lapply(gdsExp, GEOfix)
save(gdsExp, file = "~/github/stevia/Exp/gdsExp.rdata")

### functions

# convert  class "GEOData" Table method output to matrix
GEOfix <- function(df) {
 out <- lapply(df, function (vec) if(class(vec) == "character") vec <- as.numeric(vec) else vec <- as.character(vec))
 out <- as.data.frame(out, row.names = as.character(out[[1]]))[, -1]
 as.matrix(out[, -1])
}

# apply median norm to df columns
med.normalize <- function(mat) {
  out <- mat
  for (i in seq(dim(mat)[2])) { 
    vect <- mat[,i]
    med <- median(vect, na.rm = TRUE)
    out[,i] <- as.numeric(vect >= med)
  }
  return(out)
}
