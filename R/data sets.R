### http://genomics.senescence.info data sets

# genge_human needed final column name added to downloaded csv file "homologue"
genage_human <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_human.csv", stringsAsFactors = FALSE)
genage_models <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/genage_models.csv", stringsAsFactors = FALSE)
gendr_manipulations <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/gendr_manipulations.csv", stringsAsFactors = FALSE)
longevity <- read.csv("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/longevity.csv", stringsAsFactors = FALSE)
anage_data <- read.delim("C:/Users/user/Desktop/biomind/artificial biologist/stevia/sendat/anage_data.txt", stringsAsFactors = FALSE)

### tcm data base http://www.megabionet.org/tcmid/download/
temp <- tempfile()
download.file("http://www.megabionet.org/tcmid/static/downloads/prescriptions.zip",temp, mode="wb")
unzip(temp, "prescriptions.csv")
prescriptions <- read.csv("prescriptions.csv", stringsAsFactors = FALSE)