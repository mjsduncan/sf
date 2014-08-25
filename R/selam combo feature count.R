###hey mike, we count the number of appearance of each gene for mouse , pls check it if its correct we can do it for the rest organisms and off course for the whole runs.


moses_outRun1 <- lapply(mouses_1_hx5.train, Mout2str)
moses_outRun2 <- lapply(mouses_2_hx5.train, Mout2str)

Mout_Mus_1_hx5.train  <- vector("list")
Mout_Mus_2_hx5.train  <- vector("list")
Mout_all_Mus_hx5.train  <- vector("list")

Mout_Mus_1_hx5.train <- moses_outRun1[[1]]
for(i in 2:length(moses_outRun1)){
  Mout_Mus_1_hx5.train <- lbind(Mout_Mus_1_hx5.train, moses_outRun1[[i]]) }

Mout_Mus_2_hx5.train <- moses_outRun2[[1]]
for(i in 2:length(moses_outRun2)){
  Mout_Mus_2_hx5.train <- lbind(Mout_Mus_2_hx5.train, moses_outRun2[[i]]) }

# bind all freatures and results from 2 runs
Mout_all_Mus_hx5.train <- lbind(Mout_Mus_1_hx5.train, Mout_Mus_2_hx5.train)
length(unique(Mout_all_Mus_hx5.train$features$feature))
# [1] 251
# [1] 242  where does discrepency come from?  is it combo filtered via test-set-scoring column name rejection?

length(Mout_all_Mus_hx5.train$features$feature)
# [1] 333
# [1] 319  where does discrepency come from?

#### scoring the gene names
Gene_score <- matrix(nrow = length(unique(Mout_all_Mus_hx5.train$features$feature)), ncol = 2)
dimnames(Gene_score)[[2]]  <- c("GeneSymbol", "Score")
for (i in 1: length(unique(Mout_all_Mus_hx5.train$features$feature))){
 
    count <- sum(str_detect (Mout_all_Mus_hx5.train$combo, unique(Mout_all_Mus_hx5.train$features$feature)[i]))
    Gene_score[i, "GeneSymbol"] <- unique(Mout_all_Mus_hx5.train$features$feature)[i]
    Gene_score[i, "Score"]      <- count
   
}

Gene_score[1:10 ,]
#       GeneSymbol Score
#  [1,] "Fancg"    "1" 
#  [2,] "Gli3"     "100"
#  [3,] "Gltp"     "30"
#  [4,] "H13"      "9" 
#  [5,] "Glt25d1"  "100"
#  [6,] "Bphl"     "10"
#  [7,] "Glrx3"    "100"
#  [8,] "Gpr132"   "10"
#  [9,] "Fam83f"   "1" 
# [10,] "Hcrt"     "1"

# where does discripenecy come from?
# GeneSymbol      Score
#  [1,] "1810058I24Rik" "1"  
#  [2,] "2010010I01Rik" "1"  
#  [3,] "Zfp14"         "20" 
#  [4,] "Gm10229"       "20" 
#  [5,] "Zzz3"          "9"  
#  [6,] "Glt25d1"       "120"
#  [7,] "Zfp444"        "10" 
#  [8,] "Slc6a9"        "10" 
#  [9,] "Bola2"         "10" 
# [10,] "Glrx3"         "130"

### try with combo count versions
# make df of features from combo strings: calculate count of feature in all combo strings "combo count"
combo2flist2 <- function(cstr) {
require(stringr)
probe <- str_replace_all(cstr, "and+", "")
probe <- str_replace_all(probe, "or+", "")
probe <- str_replace_all(probe, "[()$]+", "")
flist <- str_split(probe, pattern = " ")
fdf <- data.frame(feature = unlist(flist), stringsAsFactors = FALSE)
fdf$low <- str_detect(fdf$feature, "!")
fdf$combo <- paste(rep(seq(length(probe)), vapply(flist, length, integer(1))), sep = "")
fdf$feature <- str_replace(fdf$feature, "!", "")
fdf <- aggregate(rep(1, dim(fdf)[1]) ~ feature + !low, data = fdf, length)
names(fdf) <- c("feature", "upregulated", "combo count")
return(fdf[order(fdf[, 1]),])
}

# make combo strings and feature dfs using combo2flist2
Mout2str2 <- function(ostr) {
require(stringr)
out <- vector("list", 3)
names(out) <- c("combo", "features", "ranks")
out[[1]] <- str_trim(str_split_fixed(ostr, " ", 2)[,2])
rank <- as.numeric(str_split_fixed(ostr, " ", 2)[,1])
out[[2]] <- combo2flist2(out[[1]])
out[[3]] <- rank
return(out)
}

# bind together results by chip using combo counting arguement
lbind2 <- function(list1, list2) {
  if(sum(sapply(list1, class) != sapply(list2, class)) > 0) return("lists element classes do not match")
  out <- vector("list", length(list1))
  names(out) <- names(list1)
  for(i in seq_along(out)) {
    if(is.atomic(list1[[i]])) {out[[i]] <- c(list1[[i]], list2[[i]])}
    if(!is.atomic(list1[[i]]))  {
      fdf <- rbind(list1[[i]], list2[[i]])
      out[[i]] <- aggregate(. ~ feature + upregulated, data = fdf, sum)
      out[[i]] <- out[[i]][order(out[[i]][, 1]),]
    }
  }
  return(out)
}

# redo combination
moses_outRun12 <- lapply(mouses_1_hx5.train, Mout2str2)
moses_outRun22 <- lapply(mouses_2_hx5.train, Mout2str2)

Mout_Mus_1_hx5.train2  <- vector("list")
Mout_Mus_2_hx5.train2  <- vector("list")
Mout_all_Mus_hx5.train2  <- vector("list")

Mout_Mus_1_hx5.train2 <- moses_outRun12[[1]]
for(i in 2:length(moses_outRun12)){
  Mout_Mus_1_hx5.train2 <- lbind2(Mout_Mus_1_hx5.train2, moses_outRun12[[i]]) }

Mout_Mus_2_hx5.train2 <- moses_outRun22[[1]]
for(i in 2:length(moses_outRun22)){
  Mout_Mus_2_hx5.train2 <- lbind2(Mout_Mus_2_hx5.train2, moses_outRun22[[i]]) }

# bind all freatures and results from 2 runs
Mout_all_Mus_hx5.train2 <- lbind(Mout_Mus_1_hx5.train2, Mout_Mus_2_hx5.train2)

merge(Mout_all_Mus_hx5.train2[[2]], as.data.frame(Gene_score), by.x = "feature", by.y = "GeneSymbol")
          feature upregulated combo count Score
1   1110057K04Rik       FALSE           9     9
2   1810058I24Rik        TRUE           1     1
3   2010010I01Rik        TRUE           1     1
4   2010107E04Rik       FALSE           1     1
5   2010111I01Rik       FALSE           1     1
6   2410002F23Rik       FALSE           1     1
7   2610002J02Rik       FALSE           1     1
8   2610201A13Rik       FALSE           1     1
9   2610507B11Rik       FALSE           1     1
10  2610524H06Rik       FALSE           1     1
11  2700094K13Rik        TRUE           1     1
12  4933415E08Rik        TRUE           1    10
13  4933415E08Rik       FALSE          10    10
14          Abhd8       FALSE           1     1
15          Acsm3        TRUE          20    10
16         Actr1a        TRUE           1     1
17          Acvr1        TRUE           1     1
18         Acvr2b       FALSE           1     1
19           Acy1       FALSE           1     1
20          Acyp1       FALSE           2    11
21          Acyp1        TRUE          10    11
22         Adam11       FALSE           1     1
23         Adam22        TRUE           1    11
24         Adam22       FALSE           1    11
25         Adam22        TRUE          10    11
26          Adcy8       FALSE          11    10
27          Adrb2       FALSE           1     1
28         Adssl1        TRUE           1     1
29         Afg3l1        TRUE           5     5
30            Afp       FALSE           2    11
31            Afp        TRUE           1    11
32            Afp        TRUE          10    11
33          Agap2        TRUE          10    10
34           Aif1        TRUE          10    10
35          Aifm1       FALSE          10    20
36          Aifm1        TRUE          10    20
37          Aifm1        TRUE          10    20
38            Aip       FALSE          20    20
39            Aip        TRUE           1    20
40           Aire        TRUE          10    10
41            Ak2        TRUE          21    20
42            Ak3        TRUE          19    19
43            Ak4       FALSE           1    10
44            Ak4        TRUE          10    10
45          Akap2       FALSE          10    10
46         Akr1b8       FALSE          12    10
47        Akr1c13        TRUE          10    10
48        Akr1c13       FALSE           2    10
49        Aldh9a1        TRUE          10    10
50         Aloxe3        TRUE          10    10
51           Ambn        TRUE          20    10
52           Amz2       FALSE          10    10
53           Amz2        TRUE           1    10
54         Anapc2       FALSE          10    10
55         Anapc2        TRUE           1    10
56           Ank3       FALSE          11    10
57        Ankrd17        TRUE           1    10
58        Ankrd17       FALSE          10    10
59          Apoa2        TRUE          10    10
60          Art2b       FALSE          11    10
61       AU020206        TRUE           1     1
62            Blm        TRUE          11    10
63          Bola2        TRUE          10    10
64         C76336        TRUE          10    10
65         C76554       FALSE           9     9
66         C77815        TRUE           7     7
67         C79242       FALSE          10    10
68         C80719        TRUE          11    10
69        Cacna1e       FALSE           8     8
70        Cacna1e        TRUE           1     8
71         Cacnb3       FALSE           8     8
72         Cacng1        TRUE           2     2
73         Cacng2        TRUE          10    10
74         Cacng2       FALSE           1    10
75          Calm4       FALSE           2     2
76         Camk2a        TRUE          11    10
77         Cenpc1       FALSE          11    10
78         Col2a1        TRUE           1    10
79         Col2a1       FALSE          11    10
80          Cpne3        TRUE           2     2
81            Cpq        TRUE           1     1
82          Cpt1a        TRUE           1     1
83          Ctrb1        TRUE          11    10
84          Ctrb1       FALSE           1    10
85        Cyp19a1        TRUE          10     9
86         Cyp2f2       FALSE          11    10
87        Cyp3a13        TRUE          10    10
88           Dio1        TRUE          10    10
89         Dnajc7        TRUE          11    10
90      DXErtd11e        TRUE          12    10
91        Dync1h1        TRUE          10    10
92         Dyrk1a        TRUE          10    10
93  E130012A19Rik       FALSE          10    10
94           Ebf1       FALSE          10    10
95           Ece1        TRUE          11    10
96            Ell       FALSE          10    10
97            Ell        TRUE           1    10
98          Ephx1        TRUE          11    10
99          Ercc2       FALSE           2     2
100         Ercc4       FALSE           6     6
101          Ern2        TRUE           6     6
102         Erp44        TRUE           1     1
103        Errfi1        TRUE          10     9
104         Esrp1       FALSE           4     4
105         Fbln1        TRUE           1     1
106         Fbln2       FALSE           1     1
107         Fgf10        TRUE          11    10
108         Fhod3        TRUE          10    10
109         Foxb1       FALSE          11    10
110         Gkap1       FALSE          11    10
111          Gli1        TRUE          29    30
112          Gli1        TRUE          11    30
113          Gli2       FALSE          10    10
114          Gli3       FALSE          51    80
115          Gli3       FALSE          40    80
116        Glipr2       FALSE          10    10
117         Glrp1        TRUE          10    30
118         Glrp1       FALSE          10    30
119         Glrp1        TRUE          10    30
120         Glrp1       FALSE          10    30
121          Glrx        TRUE          13   130
122         Glrx3       FALSE          61   130
123         Glrx3       FALSE          72   130
124       Glt25d1        TRUE          30   120
125       Glt25d1        TRUE          10   120
126       Glt25d1       FALSE          41   120
127       Glt25d1       FALSE          50   120
128          Gltp       FALSE          10    30
129          Gltp        TRUE          20    30
130          Gltp        TRUE          10    30
131          Gltp       FALSE           1    30
132         Glud1       FALSE          51    70
133         Glud1        TRUE          10    70
134         Glud1       FALSE          20    70
135          Glul        TRUE          52    40
136          Glul       FALSE           1    40
137       Glycam1       FALSE          10    30
138       Glycam1       FALSE          10    30
139       Glycam1        TRUE          20    30
140         Glyr1       FALSE          20    20
141       Gm10229        TRUE          10    20
142       Gm10229       FALSE           1    20
143       Gm10229        TRUE          10    20
144       Gm11944        TRUE          10    10
145       Gm16793        TRUE          16    30
146       Gm16793       FALSE          10    30
147       Gm16793       FALSE          10    30
148       Gm17066       FALSE          11    10
149        Hiatl1       FALSE          11    10
150         Hoxa3       FALSE          11    10
151         Htr2c       FALSE           1     1
152         Irak1        TRUE          10    10
153          Isl1        TRUE          10    10
154         Kcne1        TRUE          10    10
155         Kcne1       FALSE           1    10
156        Klk1b3        TRUE          10    10
157      Krtap6-2        TRUE          11    10
158          Lcp2        TRUE           1     1
159         Letm1       FALSE           4     4
160          Lgr5        TRUE           4     4
161          Lhpp       FALSE           7     9
162          Lhpp        TRUE           5     9
163          Lhx5        TRUE           5     5
164          Lhx6        TRUE           7     7
165          Lsm4        TRUE           1     1
166        Luc7l3        TRUE           3     3
167          Maea        TRUE          10    10
168          Mapt        TRUE          12    10
169         Matn3        TRUE          10    10
170          Melk        TRUE           1     1
171         Memo1        TRUE          10    10
172         Mesp1       FALSE           2     2
173           Met        TRUE           1     1
174         Mfap4        TRUE          10    10
175           Mia       FALSE           2     2
176        Mogat2       FALSE          10    10
177        Mrpl16        TRUE          10    10
178        Mrpl33       FALSE           1    10
179        Mrpl33        TRUE          10    10
180        Mtfr1l       FALSE          10    10
181         Ninj1       FALSE           1     1
182          Nit2       FALSE           1     1
183       Nkiras2       FALSE           1     1
184        Nkx2-3       FALSE           1     1
185        Notch2        TRUE          10    10
186         Nptx2        TRUE          11    10
187         Olfr2        TRUE           4    10
188         Olfr2       FALSE          10    10
189       Onecut1        TRUE          10    10
190         Pde1b       FALSE           1     1
191         Podxl        TRUE           4     4
192         Polg2       FALSE          10    10
193        Ppp3cc        TRUE          10    10
194         Prss2       FALSE           4     4
195         Psma5       FALSE           2     2
196         Psmc6        TRUE          16    10
197        Ptpn21       FALSE          10    10
198         Ptpn6       FALSE          10    10
199         Ptprd        TRUE           9     9
200         Ptrh1        TRUE           1     1
201         Rasd1        TRUE          10    10
202           Rb1        TRUE           2     2
203         Rbpjl        TRUE           1     1
204          Rcn2       FALSE           4     4
205           Rdx        TRUE           1     1
206           Rfk       FALSE           1     1
207         Rfwd2       FALSE           9     9
208           Rgn        TRUE          12    10
209         Rgs14        TRUE          10    10
210         Rgs19       FALSE           1     8
211         Rgs19        TRUE           7     8
212          Rgs7       FALSE          10    10
213          Rhou        TRUE          13    10
214          Rhou       FALSE          10    10
215       Sertad1       FALSE           6     6
216         Sesn1       FALSE           1     1
217        Setd1a       FALSE          10    10
218       Slc29a1       FALSE           1     1
219        Slc2a1        TRUE           1     1
220        Slc5a1       FALSE          10    10
221        Slc6a2        TRUE          11    10
222        Slc6a9        TRUE          10    10
223          Sncg        TRUE          20    20
224       Snrnp25       FALSE          10    10
225         Snx10        TRUE          10    10
226        Sprr1b        TRUE          10    10
227    St6galnac3       FALSE           1     1
228         Tesk1       FALSE           8     8
229        Tubb2a       FALSE          11    10
230          Uba7        TRUE           2     2
231         Ubac1        TRUE           1     1
232       Ube2d2a        TRUE           1     1
233         Ube2f       FALSE           3     3
234         Ube2t       FALSE          10    10
235         Ube4b        TRUE           1     1
236         Ubfd1       FALSE           4     4
237          Ubl3        TRUE           1     1
238          Vapa        TRUE           1     1
239       Vipas39        TRUE           1     1
240         Vprbp       FALSE           1     1
241         Wbp1l        TRUE           1     1
242         Wdr77       FALSE           1     1
243          Wnt3       FALSE          10    20
244         Wnt3a        TRUE          11    10
245          Wnt4        TRUE          10    10
246         Wnt8a       FALSE           1     1
247          Xpot       FALSE           4     4
248         Ywhag        TRUE          10    10
249         Zap70       FALSE          10    10
250        Zbtb14       FALSE           9     9
251        Zbtb16        TRUE           7     7
252        Zbtb17       FALSE           2     2
253        Zbtb22        TRUE          19    19
254        Zbtb25       FALSE          10    10
255        Zbtb46       FALSE           1     1
256        Zbtb48       FALSE          15    13
257        Zbtb7a       FALSE           7     7
258        Zbtb7b        TRUE           1     1
259       Zc3h12c       FALSE           8     8
260          Zfp1       FALSE           1    21
261         Zfp14       FALSE          20    20
262        Zfp213        TRUE          10    10
263         Zfp26       FALSE           1     1
264        Zfp277       FALSE           3     3
265         Zfp28        TRUE          10    10
266        Zfp358       FALSE           1     1
267        Zfp398       FALSE          10    10
268        Zfp444       FALSE          10    10
269          Zic4       FALSE          10    10
270        Zmynd8       FALSE           1     3
271        Zmynd8        TRUE           2     3
272         Znfx1       FALSE           1     1
273        Znhit1        TRUE           2    13
274        Znhit1        TRUE          10    13
275        Znhit1       FALSE           1    13
276        Znhit2        TRUE          10    15
277        Znhit2        TRUE           5    15
278         Znrf1       FALSE           1    12
279         Znrf1        TRUE          11    12
280         Znrf4       FALSE           5     6
281         Znrf4        TRUE           1     6
282           Zp1        TRUE           2     2
283           Zp2        TRUE           2     3
284           Zp2       FALSE           1     3
285           Zp3        TRUE           2    11
286           Zp3        TRUE           3    11
287          Zp3r        TRUE           6     6
288        Zranb1        TRUE           1     1
289         Zrsr1        TRUE           1     1
290         Zrsr2        TRUE           1     1
291        Zscan2       FALSE           1     1
292        Zswim1       FALSE           1     2
293        Zswim1        TRUE           1     2
294        Zswim7        TRUE          11    11
295          Zw10        TRUE           2     1
296           Zyx       FALSE           1     1
297          Zzz3       FALSE           9     9

#### scoring the gene names
Gene_score1 <- matrix(nrow = length(unique(moses_outRun1[[2]]$features$feature)), ncol = 2)
dimnames(Gene_score1)[[2]]  <- c("GeneSymbol", "Score")
for (i in 1: length(unique(moses_outRun1[[2]]$features$feature))){
 
    count <- sum(str_detect (moses_outRun1[[2]]$combo, unique(moses_outRun1[[2]]$features$feature)[i]))
    Gene_score1[i, "GeneSymbol"] <- unique(moses_outRun1[[2]]$features$feature)[i]
    Gene_score1[i, "Score"]      <- count
   
}
merge(moses_outRun12[[2]][[2]], as.data.frame(Gene_score1), by.x = "feature", by.y = "GeneSymbol")
   feature upregulated combo count Score
1      Blm        TRUE          11    10
2   Dnajc7        TRUE          11    10
3     Ece1        TRUE          11    10
4      Ell       FALSE          10    10
5      Ell        TRUE           1    10
6    Glyr1       FALSE          10    10
7  Gm16793       FALSE          10    10
8   Ptpn21       FALSE          10    10
9    Ptpn6       FALSE          10    10
10   Ptprd        TRUE           9     9
11   Ptrh1        TRUE           1     1
12   Ube4b        TRUE           1     1
13   Ubfd1       FALSE           4     4
14    Ubl3        TRUE           1     1

moses_outRun12[[2]]
$combo
 [1] "or(and(or(and(!$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21 !$Ubfd1) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"               
 [2] "or(and(or(and(or(!$Ptpn6 !$Ubfd1) $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"           
 [3] "or(and(or(and(!$Ptpn6 $Ptprd) and($Ube4b !$Ubfd1) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"   
 [4] "or(and(or(and(or(and(!$Ubfd1 $Ubl3) !$Ptpn6) $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"
 [5] "or(and(or(and(!$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"                       
 [6] "or(and(or(and(!$Ptpn6 $Ptrh1) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"                       
 [7] "or(and(or(and($Blm !$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"                  
 [8] "or(and(or(and($Dnajc7 !$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"               
 [9] "or(and(or(and($Ece1 !$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"                 
[10] "or(and(or(and($Ell !$Ptpn6 $Ptprd) !$Gm16793 !$Ptpn21) or($Blm $Dnajc7 $Ece1)) and(!$Ell !$Glyr1))"                  

$features
   feature upregulated combo count
7      Blm        TRUE          11
8   Dnajc7        TRUE          11
9     Ece1        TRUE          11
1      Ell       FALSE          10
10     Ell        TRUE           1
2    Glyr1       FALSE          10
3  Gm16793       FALSE          10
4   Ptpn21       FALSE          10
5    Ptpn6       FALSE          10
11   Ptprd        TRUE           9
12   Ptrh1        TRUE           1
13   Ube4b        TRUE           1
6    Ubfd1       FALSE           4
14    Ubl3        TRUE           1
