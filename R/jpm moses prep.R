### start with jpmExp from jpmData.rdata.  6 arrays cover 20 data sets
# make moses data sets
# fix ages$lung
ages$lung <- ages$lung[1:15]

mosesExp <- vector("list", 26)
names(mosesExp) <- names(jpmExp)
for(n in names(mosesExp)) {
  mosesExp[[n]] <- rbind(age = t(as.matrix(ages[[n]])), jpmExp[[n]])
  mosesExp[[n]] <- med.normalize(t(mosesExp[[n]]))
  dimnames(mosesExp[[n]])[[2]][1] <- "age"
}

# check age binary
mapply(function(x, y) rbind(x, y[,1]), ages, mosesExp)
# $h_brain
#   GSM27015 GSM27016 GSM27018 GSM27021 GSM27023 GSM27024 GSM27025 GSM27027 GSM27028 GSM27031 GSM27032 GSM27034 GSM27035 GSM27036 GSM27038 GSM27040 GSM27042
# x       26       26       29       37       40       42       45       52       53       66       70       73       77       80       85       90       91
#          0        0        0        0        0        0        0        0        0        1        1        1        1        1        1        1        1
#   GSM27043 GSM27017 GSM27019 GSM27020 GSM27022 GSM27026 GSM27029 GSM27030 GSM27033 GSM27037 GSM27039 GSM27041 GSM27044
# x       95       27       30       36       38       48       56       61       71       81       87       90      106
#          1        0        0        0        0        0        0        1        1        1        1        1        1
mosesExp$h_brain[c(5:9, 23:24),1] <- 1
# 
# $muscle1
#   GSM10374 GSM10375 GSM10376 GSM10377 GSM10378 GSM10379 GSM10380 GSM10381 GSM10382 GSM10383 GSM10384 GSM10385 GSM10386 GSM10387 GSM10388
# x       23       24       21       26       20       29       22       70       71       66       65       68       68       66       67
#          0        0        0        0        0        0        0        1        1        1        1        1        1        1        1
# 
# $muscle2
#   GSM10354 GSM10355 GSM10356 GSM10359 GSM10360 GSM10361 GSM10362 GSM10363 GSM10364 GSM10365 GSM10366 GSM10367 GSM10368 GSM10369 GSM10370
# x       23       24       21       26       20       29       22       70       71       66       65       68       68       66       67
#          0        0        0        0        0        0        0        1        1        1        1        1        1        1        1
# 
# $muscle3
#   GSM5254 GSM5255 GSM5256 GSM5257 GSM5259 GSM5260 GSM5261 GSM5262 GSM5263 GSM5264 GSM5265 GSM5266 GSM5268 GSM5270 GSM5271
# x      25      27      24      26      23      25      21      69      69      68      73      67      74      71      75
#         0       0       0       0       0       0       0       1       1       1       1       1       1       1       1
# 
# $muscle4
#   GSM5300 GSM5301 GSM5302 GSM5303 GSM5305 GSM5306 GSM5307 GSM5308 GSM5309 GSM5310 GSM5311 GSM5312 GSM5313 GSM5314 GSM5315
# x      25      27      24      26      23      25      21      69      69      68      73      67      74      71      75
#         0       0       0       0       0       0       0       1       1       1       1       1       1       1       1
# 
# $muscle5
#   GSM2390 GSM2391 GSM2392 GSM2393 GSM2394 GSM2395 GSM2396 GSM2397 GSM2398 GSM2399 GSM2400 GSM2401
# x      26      26      26      26      26      26    69.5    69.5    69.5    69.5    69.5    69.5
#         0       0       0       0       0       0     1.0     1.0     1.0     1.0     1.0     1.0
# 
# $muscle
#   GSM146340 GSM146341 GSM146342 GSM146343 GSM146344 GSM146345 GSM146346 GSM146347 GSM146348 GSM146349
# x         5         5         5         5         5        25        25        25        25        25
#           0         0         0         0         0         1         1         1         1         1
# 
# $kidney1
#   GSM7467 GSM7468 GSM7469 GSM7470 GSM7471 GSM7457 GSM7459 GSM7461 GSM7463 GSM7465
# x       5       5       5       5       5      30      30      30      30      30
#         0       0       0       0       0       1       1       1       1       1
# 
# $kidney2
#   GSM7472 GSM7473 GSM7474 GSM7475 GSM7476 GSM7458 GSM7460 GSM7462 GSM7464 GSM7466
# x       5       5       5       5       5      30      30      30      30      30
#         0       0       0       0       0       1       1       1       1       1
# 
# $m_brain
#   GSM73001 GSM73014 GSM73015 GSM73002 GSM73016 GSM73017
# x      4.5      4.5      4.5       22       22       22
#        0.0      0.0      0.0        1        1        1
mosesExp$m_brain[1:3, 1] <- 0
# 
# $m_hippo
#   GSM114426 GSM114427 GSM114428 GSM114430 GSM114431 GSM114433 GSM114434 GSM114436 GSM114437 GSM114391 GSM114393 GSM114394 GSM114396 GSM114397 GSM114400
# x         2         2         2         2         2         2         2         2         2        15        15        15        15        15        15
#           0         0         0         0         0         0         0         0         0         1         1         1         1         1         1
#   GSM114402 GSM114409 GSM114410 GSM114412 GSM114414 GSM114417 GSM114424 GSM114425
# x        15        15        15        15        15        15        15        15
#           1         1         1         1         1         1         1         1
# 
# $liver
#   GSM69713 GSM69714 GSM69715 GSM69716 GSM69707 GSM69708 GSM69709
# x        6        6        6        6       22       22       22
#          1        1        1        1        1        1        1
mosesExp$liver[1:4, 1] <- 0
# 
# $m_heart
#   GSM2334 GSM2335 GSM2336 GSM2186 GSM2187 GSM2188 GSM2180 GSM2181 GSM2182 GSM2337 GSM2338 GSM2339
# x       3       3       3       5       5       5      12      12      12      12      12      12
#         0       0       0       0       0       0       1       1       1       1       1       1
# 
# $lung
#   GSM152256 GSM152257 GSM152258 GSM152259 GSM152260 GSM152261 GSM152262 GSM152263 GSM152264 GSM152265 GSM152266 GSM152267 GSM152268 GSM152269 GSM152270
# x         2         2         2        18        18        18        26        26        26        26        26        26         2         2         2
#           0         0         0         1         1         1         1         1         1         1         1         1         0         0         0
# 
# $cochlea
#   GSM108106 GSM108107 GSM108108 GSM108103 GSM108104 GSM108105
# x         4         4         4        15        15        15
#           0         0         0         1         1         1
# 
# $hemato_stem
#   GSM98881 GSM98882 GSM98883 GSM98876 GSM98877 GSM98878 GSM98879 GSM98880
# x      2.5      2.5      2.5     22.5     22.5     22.5     22.5     22.5
#        0.0      0.0      0.0      1.0      1.0      1.0      1.0      1.0
mosesExp$hemato_stem[, 1] <- c(0, 0, 0, 1, 1, 1, 1, 1)
# 
# $myo_progen
#   GSM3810 GSM3811 GSM3812 GSM3813
# x       8       8      23      23
#         0       0       1       1
# 
# $r_hippo
#   GSM132501 GSM132509 GSM132510 GSM132511 GSM132525 GSM132526 GSM132527 GSM132528 GSM132529 GSM132530 GSM132486 GSM132505 GSM132506 GSM132507 GSM132544
# x         5         5         5         5         5         5         5         5         5         5         5         5         5         5         5
#           0         0         0         0         0         0         0         0         0         0         0         0         0         0         0
#   GSM132545 GSM132546 GSM132547 GSM132548 GSM132549 GSM132489 GSM132490 GSM132491 GSM132492 GSM132493 GSM132502 GSM132503 GSM132504 GSM132543 GSM132500
# x         5         5         5         5         5         5         5         5         5         5         5         5         5         5        25
#           0         0         0         0         0         0         0         0         0         0         0         0         0         0         1
#   GSM132518 GSM132519 GSM132523 GSM132524 GSM132557 GSM132558 GSM132559 GSM132560 GSM132561 GSM132488 GSM132495 GSM132496 GSM132497 GSM132498 GSM132499
# x        25        25        25        25        25        25        25        25        25        25        25        25        25        25        25
#           1         1         1         1         1         1         1         1         1         1         1         1         1         1         1
#   GSM132521 GSM132537 GSM132539 GSM132540 GSM132484 GSM132485 GSM132494 GSM132512 GSM132513 GSM132520 GSM132522 GSM132533 GSM132536 GSM132541 GSM132487
# x        25        25        25        25        25        25        25        25        25        25        25        25        25        25        25
#           1         1         1         1         1         1         1         1         1         1         1         1         1         1         1
#   GSM132508 GSM132515 GSM132538 GSM132542 GSM132550 GSM132551 GSM132552 GSM132554 GSM132556 GSM132514 GSM132516 GSM132517 GSM132531 GSM132532 GSM132534
# x        25        25        25        25        25        25        25        25        25        25        25        25        25        25        25
#           1         1         1         1         1         1         1         1         1         1         1         1         1         1         1
#   GSM132535 GSM132553 GSM132555
# x        25        25        25
#           1         1         1
# 
# $stromal
#   GSM75444 GSM75445 GSM75446
# x        3        3       15
#          1        1        1
mosesExp$stromal[1:2, 1] <- 0
# 
# $spinal_cord
#   GSM74342 GSM74343 GSM74344 GSM74345 GSM74346 GSM74347 GSM74348 GSM74349 GSM74350
# x        6        6        6       18       18       18       30       30       30
#          0        0        0        1        1        1        1        1        1
# 
# $oculomotor
#   GSM74333 GSM74334 GSM74335 GSM74336 GSM74337 GSM74338 GSM74339 GSM74340 GSM74341
# x        6        6        6       18       18       18       30       30       30
#          0        0        0        1        1        1        1        1        1
# 
# $skeletal_ms
#   GSM74432 GSM74433 GSM74434 GSM74435 GSM74436 GSM74437 GSM74438 GSM74439 GSM74440 GSM74441 GSM74442 GSM74443
# x        6        6        6        6       18       18       18       18       30       30       30       30
#          0        0        0        0        1        1        1        1        1        1        1        1
# 
# $extraoc_ms
#   GSM74444 GSM74445 GSM74446 GSM74447 GSM74448 GSM74449 GSM74450 GSM74451 GSM74452 GSM74453 GSM74454 GSM74455
# x        6        6        6        6       18       18       18       18       30       30       30       30
#          0        0        0        0        1        1        1        1        1        1        1        1
# 
# $laryngeal_ms
#   GSM74456 GSM74457 GSM74458 GSM74459 GSM74460 GSM74461 GSM74462 GSM74463 GSM74464
# x        6        6        6       18       18       18       30       30       30
#          0        0        0        1        1        1        1        1        1
# 
# $r_heart
#   GSM6174 GSM6175 GSM6176 GSM6177 GSM6178 GSM6168 GSM6169 GSM6170 GSM6171 GSM6172 GSM6173
# x     3.5     3.5     3.5     3.5     3.5      21      21      21      21      21      21
#       0.0     0.0     0.0     0.0     0.0       1       1       1       1       1       1
mosesExp$r_heart[1:5, 1] <- 0
# 
# $CA1_hipp2
#   GSM13323 GSM13324 GSM13325 GSM13326 GSM13327 GSM13328 GSM13329 GSM13330 GSM13331 GSM13313 GSM13314 GSM13315 GSM13316 GSM13317 GSM13318 GSM13319 GSM13320
# x        4        4        4        4        4        4        4        4        4       14       14       14       14       14       14       14       14
#          0        0        0        0        0        0        0        0        0        1        1        1        1        1        1        1        1
#   GSM13321 GSM13322 GSM13303 GSM13304 GSM13305 GSM13306 GSM13307 GSM13308 GSM13309 GSM13310 GSM13311 GSM13312
# x       14       14       24       24       24       24       24       24       24       24       24       24
#          1        1        1        1        1        1        1        1        1        1        1        1

# export file
setwd("~/projects/sf_results")
save(mosesExp, file = "mosesExp.rdata")

# merge data sets using same microarray chip
gpl1261 <- merge(t(mosesExp[["m_brain"]]), t(mosesExp[["cochlea"]]), by = "row.names", all = TRUE)
row.names(gpl1261) <- gpl1261$Row.names
gpl1261$Row.names <- NULL
gpl1261 <- merge(gpl1261, t(mosesExp[["hemato_stem"]]), by = "row.names", all = TRUE)
row.names(gpl1261) <- gpl1261$Row.names
gpl1261$Row.names <- NULL
gpl1261 <- merge(gpl1261, t(mosesExp[["lung"]]), by = "row.names", all = TRUE)
row.names(gpl1261) <- gpl1261$Row.names
gpl1261$Row.names <- NULL
gpl1261 <- t(as.matrix(rbind(gpl1261["age",], gpl1261[row.names(gpl1261) != "age",])))
write.csv(gpl1261, file = "~/projects/sf_results/gpl1261.csv", col.names = TRUE, row.names = FALSE)

gpl341 <- merge(t(mosesExp[["stromal"]]), t(mosesExp[["spinal_cord"]]), by = "row.names", all = TRUE)
row.names(gpl341) <- gpl341$Row.names
gpl341$Row.names <- NULL
gpl341 <- merge(gpl341, t(mosesExp[["r_hippo"]]), by = "row.names", all = TRUE)
row.names(gpl341) <- gpl341$Row.names
gpl341$Row.names <- NULL
gpl341 <- merge(gpl341, t(mosesExp[["oculomotor"]]), by = "row.names", all = TRUE)
row.names(gpl341) <- gpl341$Row.names
gpl341$Row.names <- NULL
gpl341 <- t(as.matrix(rbind(gpl341["age",], gpl341[row.names(gpl341) != "age",])))

# gpl341[["X1370846_at"]] <-NULL
write.csv(gpl341, file = "~/projects/sf_results/gpl341.csv", col.names = TRUE, row.names = FALSE)

gpl85 <- merge(t(mosesExp[["laryngeal_ms"]]), t(mosesExp[["skeletal_ms"]]), by = "row.names", all = TRUE)
row.names(gpl85) <- gpl85$Row.names
gpl85$Row.names <- NULL
gpl85 <- merge(gpl85, t(mosesExp[["r_heart"]]), by = "row.names", all = TRUE)
row.names(gpl85) <- gpl85$Row.names
gpl85$Row.names <- NULL
gpl85 <- merge(gpl85, t(mosesExp[["CA1_hipp2"]]), by = "row.names", all = TRUE)
row.names(gpl85) <- gpl85$Row.names
gpl85$Row.names <- NULL
gpl85 <- merge(gpl85, t(mosesExp[["extraoc_ms"]]), by = "row.names", all = TRUE)
row.names(gpl85) <- gpl85$Row.names
gpl85$Row.names <- NULL
gpl85 <- t(as.matrix(rbind(gpl85["age",], gpl85[row.names(gpl85) != "age",])))
write.csv(gpl85, file = "~/projects/sf_results/gpl85.csv", col.names = TRUE, row.names = FALSE)

gpl81 <- merge(t(mosesExp[["myo_progen"]]), t(mosesExp[["liver"]]), by = "row.names", all = TRUE)
row.names(gpl81) <- gpl81$Row.names
gpl81$Row.names <- NULL
gpl81 <- merge(gpl81, t(mosesExp[["m_hippo"]]), by = "row.names", all = TRUE)
row.names(gpl81) <- gpl81$Row.names
gpl81$Row.names <- NULL
gpl81 <- t(as.matrix(rbind(gpl81["age",], gpl81[row.names(gpl81) != "age",])))
write.csv(gpl81, file = "~/projects/sf_results/gpl81.csv", col.names = TRUE, row.names = FALSE)

gpl96 <- merge(t(mosesExp[["muscle3"]]), t(mosesExp[["muscle1"]]), by = "row.names", all = TRUE)
row.names(gpl96) <- gpl96$Row.names
gpl96$Row.names <- NULL
gpl96 <- t(as.matrix(rbind(gpl96["age",], gpl96[row.names(gpl96) != "age",])))
write.csv(gpl96, file = "~/projects/sf_results/gpl96.csv", col.names = TRUE, row.names = FALSE)

gpl97 <- merge(t(mosesExp[["muscle2"]]), t(mosesExp[["muscle4"]]), by = "row.names", all = TRUE)
row.names(gpl97) <- gpl97$Row.names
gpl97$Row.names <- NULL
gpl97 <- t(as.matrix(rbind(gpl97["age",], gpl97[row.names(gpl97) != "age",])))
write.csv(gpl97, file = "~/projects/sf_results/gpl97.csv", col.names = TRUE, row.names = FALSE)

# # example of combining data sets using Reduce
# gpl.data <- vector("list", 6)
# names(gpl.data) <- unique(arrays$to_acc)[c(1, 4, 7, 9, 11:12)]
# for(n in names(gpl.data)) gpl.data[[n]] <- Reduce(cbind, gdsExp[get(n)])
# sapply(gpl.data, dim)
