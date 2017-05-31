library(dplyr)
library(tibble)
library(reshape2)
library(magrittr)
library(tidyr)
library(ggplot2)
library(XLConnect)
library(DataCombine)
library(vegan)
library("corrgram")
library(tibble)
library(DataCombine)
library(scales)
###Algandmete analyys
setwd("C:/Users/kr66t/Dropbox/kroot/Ylevaatamiseks")
algne_tabel<-read.csv2("algne.csv")%>% drop_na()
#v6tan v4lja esimese sekveneerimise tabeli
esimene_sekv<-algne_tabel%>%
  select(-KARS1, -KARS2, -KARS3, -KARS4, -KARS5, -KARK1, -KARK2,
         -KARK3, -KARK4, -KARK5, -KLK1, -KLK2, -KLK3, -KLK4, -KLK5,
         -KLS1, -KLS2, -KLS3, -KLS4, -KLS5, -PEK1, -PEK2, -PEK3, -PEK4,
         -PEK5, -PES1, -PES2, -PES3, -PES4, -PES5, -PORK1, -PORK2, -PORK3,
         -PORK4, -PORK5, -PORS1, -PORS2, -PORS3, -PORS4, -PORS5, -TK1,
         -TK2, -TK3, -TK4, -TK5, -TS1, -TS2, -TS3, -TS4, -TS5) 
#v6tan v4lja teise sekveneerimise tabeli
teine_sekv<-algne_tabel%>%
  select(-kars1, -kars2, -kars3, -kars4, -kars5, -kark1, -kark2, 
         -kark3, -kark4, -kark5, -kalk1, -kalk2, -kalk3, -kalk4, 
         -kalk5, -kals1, -kals2, -kals3, -kals4, -kals5, -pek1, 
         -pek2, -pek3, -pek4, -pek5, -pes1, -pes2, -pes3, -pes4, 
         -pes5, -pork1, -pork2, -pork3, -pork4, -pork5, -pors1, 
         -pors2, -pors3, -pors4, -pors5, -topk1, -topk2, -topk3, 
         -topk4, -topk5, -tops1, -tops2, -tops3, -tops4, -tops5)
#?hendan samade proovide andmed
y2<-teine_sekv%>%
  rename(kars1=KARS1, kars2=KARS2, kars3=KARS3, kars4=KARS4, kars5=KARS5, 
         kark1=KARK1, kark2=KARK2, kark3=KARK3, kark4=KARK4, kark5=KARK5, 
         kalk1=KLK1, kalk2=KLK2, kalk3=KLK3, kalk4=KLK4, kalk5=KLK5, 
         kals1=KLS1, kals2=KLS2, kals3=KLS3, kals4=KLS4, kals5=KLS5, pek1=PEK1, 
         pek2=PEK2, pek3=PEK3, pek4=PEK4, pek5=PEK5, pes1=PES1, pes2=PES2, 
         pes3=PES3, pes4=PES4, pes5=PES5, pork1=PORK1, pork2=PORK2, 
         pork3=PORK3, pork4=PORK4, pork5=PORK5, pors1=PORS1, pors2=PORS2, 
         pors3=PORS3, pors4=PORS4, pors5=PORS5, topk1=TK1, topk2=TK2, topk3=TK3, 
         topk4=TK4, topk5=TK5, tops1=TS1, tops2=TS2, tops3=TS3, tops4=TS4, 
         tops5=TS5 )
y1<-esimene_sekv
tabel1 <- y1[order(y1$X), ]
tabel2 <- y2[order(y2$X), ]
nveerg <- ncol(tabel1)
#liidan kahes tabelis vastavad countid
k6ikkokku <- tabel1[ , 3:nveerg] + tabel2[ , 3:nveerg]
#lisan tabeli lÃµppu 2 veergu: OTUD (X) ja TAXON.ID
k6ikkokku$X <- tabel1$X
k6ikkokku$TAXON.ID <- tabel1$TAXON.ID

##Esimese sekveneerimise tulemuste anal??s, algtabel esimene_sekv
#eemaldan ainult nullidega OTUd
esimene_ridade_rsum<-rowSums(esimene_sekv[-(1:2)])
esimene_ridade_rsum<-as.data.frame(esimene_ridade_rsum)
#liidan ridade summade kolumni esimese sekveneerimise tabelile 
esimene_sekv_sum<-cbind(esimene_sekv, esimene_ridade_rsum)
esimene_sekv_sum<-esimene_sekv_sum%>%rename(summa=esimene_ridade_rsum)
#eemaldan nullidega read
esimene_sekv_0ta<-esimene_sekv_sum%>%filter(esimene_sekv_sum$summa>0)
#eemaldan kolumni summa veeru, saan tabeli, kus ei ole nulli ridasid
esimene_sekveneerimine<-select(esimene_sekv_0ta,-summa)

##leian proovide kohta keskmise lugemite arvu ja mediaani
##selleks teen tabelist maatriksi
esimene<-esimene_sekveneerimine[, 3:52]
b2<-colSums(esimene)%>%as.data.frame()
summary(b2)
#arvutan, mitu erinevat OTU-t on proovides, kus on lugemeid rohkem >0, >10 ja >100,K6IK ANDMED
esimene_OTUde_jaotus<-esimene_sekveneerimine%>%
  melt(id.vars="X", measure.vars=names(esimene_sekveneerimine[,-(1:2)]))%>%
  group_by(variable)%>%
  summarise(length0=length(value[value>0]), 
            length10=length(value[value>10]), 
            length100=length(value[value>100]))

##liidan lugemite summa kolumni OTUde jaotuse tabelile
#kasutan eelnevat tehtud kolumnisummade tabelit b2
#teen kolumni "katse"?
esimene_kolumn_sum<-b2%>%rownames_to_column(var = "katse")
esimene_kolumn_sum$summa<-esimene_kolumn_sum$.
esimene_kolumn_sum<-esimene_kolumn_sum%>%select(katse, summa)
###liidan kaks tabelit- lugemite arvud; lugemite summa
esimene_OTUd_summa<-full_join(esimene_OTUde_jaotus, esimene_kolumn_sum, by=c("variable" = "katse"))
esimene_OTUd_summa<-as.data.frame(esimene_OTUd_summa)

# >100 ja >10 lugemite arvutused
#>100 lugemit
idx=esimene>100
esimenerohkemkui100lugem<-colSums(esimene*idx)
esimenerohkemkui100lugem<-esimenerohkemkui100lugem%>%as.data.frame()
esimenerohkemkui100lugem$sumrohkem100<-esimenerohkemkui100lugem$.
esimenerohkemkui100lugem<-select(esimenerohkemkui100lugem,sumrohkem100)
#teen variable kolumni
esimenerohkemkui100lugem$variable<- rownames(esimenerohkemkui100lugem)
#panen numbrid rownamedeks
rownames(esimenerohkemkui100lugem)<- 1:nrow(esimenerohkemkui100lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 100
esimene_OTUd_summa<-full_join(esimene_OTUd_summa,esimenerohkemkui100lugem)

#>10 lugemit
idx2=esimene>10
esimenerohkemkui10lugem<-colSums(esimene*idx2)
esimenerohkemkui10lugem<-esimenerohkemkui10lugem%>%as.data.frame()
esimenerohkemkui10lugem$sumrohkem10<-esimenerohkemkui10lugem$.
esimenerohkemkui10lugem<-select(esimenerohkemkui10lugem,sumrohkem10)
#teen variable kolumni
esimenerohkemkui10lugem$variable<- rownames(esimenerohkemkui10lugem)
#panen numbrid rownamedeks
rownames(esimenerohkemkui10lugem)<- 1:nrow(esimenerohkemkui10lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 10
esimene_OTUd_summa<-full_join(esimene_OTUd_summa,esimenerohkemkui10lugem)
#panen kolumnid j?rjekorda
esimene_OTUd_summa<-esimene_OTUd_summa[, c(1,2,3,4,5,7,6)]

### arvutan kui palju sekventse j??b proovides alles peale mitokondrite ja kloroplastide eemaldamist
library(DataCombine)
esimene_ilma_klor_mit<-grepl.sub(data=esimene_sekveneerimine, 
                                 pattern = c("Chloroplast","Mitochondria"),
                                 Var= "TAXON.ID", 
                                 keep.found=F)
##leian palju sekventse j?i alles peale mitokondrite ja kloroplastide eemaldamist
esimene_ilma_kolsum<-colSums(esimene_ilma_klor_mit[-(1:2)])
esimene_ilma_kolsum<-as.data.frame(esimene_ilma_kolsum)
esimene_ilma_kolsum$variable<-rownames(esimene_ilma_kolsum)
esimene_OTUd_summa<-full_join(esimene_OTUd_summa,esimene_ilma_kolsum)
library(XLConnect)
#kirjutan esimese kokkuv?tva tabeli excelisse
writeWorksheetToFile("Esimene_sekv_alganalyys.xlsx", data=esimene_OTUd_summa, sheet="tabel", startRow=1, startCol=1)
##k6ik esindatud h6imkonnad esimeses sekveneerimises
library(tidyr)
kogu_ID<-esimene_ilma_klor_mit%>%
  separate(TAXON.ID,
           c("riik", "hmk", "klass", "selts", "sgk", "prk", "liik"), 
           sep=";", 
           extra="merge")
kogu_ID$hmk<-as.factor(kogu_ID$hmk)
levels(kogu_ID$hmk)
###1] "Acidobacteria"       "Actinobacteria"      "Armatimonadetes"     "Bacteroidetes"       "Chlamydiae"         
##[6] "Chloroflexi"         "Deinococcus-Thermus" "Firmicutes"          "Gemmatimonadetes"    "Gracilibacteria"    
###[11] "Planctomycetes"      "Proteobacteria"      "Saccharibacteria"    "Tenericutes"         "Thaumarchaeota"     
###[16] "TM6"                 "Verrucomicrobia"  
#17 h6imkonda


#esimese sekveneerimise keskmine OTUde arv proovides
summary(esimene_OTUd_summa$length0)
#k6ikide sekventside summa 
esimese_sekv_read_SUM<-sum(esimene_OTUd_summa$summa)
#kloroplastide lugemite summa
esimeneklor<-grepl.sub(data = esimene_sekveneerimine, pattern = c("Chloroplast"), Var = "TAXON.ID", keep.found=T)
esimeneklor2<-esimeneklor[,3:52]
esimeneklor_sum<-sum(esimeneklor2)
#Arvutan mitokondrite lugemite summa
esimene_mitokondrid<-grepl.sub(data = esimene_sekveneerimine, pattern = c("Mitochondria"), Var = "TAXON.ID", keep.found=T)
esimene_mitokondrid_2<-esimene_mitokondrid[,3:52]
esimene_mit_sum<-sum(esimene_mitokondrid_2)
#arvutan, mitu OTU-t on kloroplastide seas
esimeneklor_OTUd<-length(esimeneklor$X)
#arvutan mitokondrite OTUd
esimenemitokondrid_OTU<-length(esimene_mitokondrid$X)
#arvutan kui palju on arhesid
esimene_arhed<-grepl.sub(data = esimene_sekveneerimine, pattern = c("Archaea"), Var = "TAXON.ID", keep.found=T)
esimene_arhed<-esimene_arhed[,3:52]
esimene_arhed_SUM<-sum(esimene_arhed)
#k6ik OTU-d
esimene_k6ik_OTUd<-length(esimene_sekveneerimine$X)
#ilma mitokondrite ja kloroplastideta OTU-d
esimene_ilma_mitokon_kloroplastiteta_OTUd<-esimene_k6ik_OTUd-esimenemitokondrid_OTU-esimeneklor_OTUd
#ilma mitokondrite ja kloroplastiteta sekventsid
esimese_sekv_lugemid_ilma_mit_klor_ta<-esimese_sekv_read_SUM-esimene_mit_sum-esimeneklor_sum



##Teise sekveneerimise tulemuste anal??s, algtabel teine_sekv
#eemaldan ainult nullidega OTUd
teine_ridade_rsum<-rowSums(teine_sekv[-(1:2)])
teine_ridade_rsum<-as.data.frame(teine_ridade_rsum)
#liidan ridade summade kolumni esimesele 
teine_sekv_sum<-cbind(teine_sekv, teine_ridade_rsum)
teine_sekv_sum<-teine_sekv_sum%>%rename(summa=teine_ridade_rsum)
#eemaldan nullidega read
teine_sekv_0ta<-teine_sekv_sum%>%filter(teine_sekv_sum$summa>0)
#eemaldan kolumni summa
teine_sekveneerimine<-select(teine_sekv_0ta,-summa)
##leian proovide kohta keskmise lugemite arvu
teine<-teine_sekveneerimine[, 3:52]
b2<-colSums(teine)%>%as.data.frame()
summary(b2)
#arvutan, mitu erinevat OTU-t on proovides, kus on lugemeid rohkem >0, >10 ja >100,K6IK ANDMED
teine_OTUde_jaotus<-teine_sekveneerimine%>%melt(id.vars="X", measure.vars=names(teine_sekveneerimine[,-(1:2)]))%>%
  group_by(variable)%>%
  summarise(length0=length(value[value>0]), length10=length(value[value>10]), length100=length(value[value>100]))

#liidan lugemite summa kolumni OTUde jaotuse tabelile
#kasutan eelnevat tehtud kolumnisummade tabelit b2
#teen kolumni "katse"?
teine_kolumn_sum<-b2%>%rownames_to_column(var = "katse")
teine_kolumn_sum$summa<-teine_kolumn_sum$.
teine_kolumn_sum<-teine_kolumn_sum%>%select(katse, summa)
###liidan kaks tabelit- OTUde arvud ja lugemite summa
teine_OTUd_summa<-full_join(teine_OTUde_jaotus, teine_kolumn_sum, by=c("variable" = "katse"))
teine_OTUd_summa<-as.data.frame(teine_OTUd_summa)

#lugemite arvutamised
#>100 lugemit
idx=teine>100
teinerohkemkui100lugem<-colSums(teine*idx)
teinerohkemkui100lugem<-teinerohkemkui100lugem%>%as.data.frame()
teinerohkemkui100lugem$sumrohkem100<-teinerohkemkui100lugem$.
teinerohkemkui100lugem<-select(teinerohkemkui100lugem,sumrohkem100)
#teen variable kolumni
teinerohkemkui100lugem$variable<- rownames(teinerohkemkui100lugem)
#panen numbrid rownamedeks
rownames(teinerohkemkui100lugem)<- 1:nrow(teinerohkemkui100lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 100
teine_OTUd_summa<-full_join(teine_OTUd_summa,teinerohkemkui100lugem)

#>10 lugemit
idx2=teine>10
teinerohkemkui10lugem<-colSums(teine*idx2)
teinerohkemkui10lugem<-teinerohkemkui10lugem%>%as.data.frame()
teinerohkemkui10lugem$sumrohkem10<-teinerohkemkui10lugem$.
teinerohkemkui10lugem<-select(teinerohkemkui10lugem,sumrohkem10)
#teen variable kolumni
teinerohkemkui10lugem$variable<- rownames(teinerohkemkui10lugem)
#panen numbrid rownamedeks
rownames(teinerohkemkui10lugem)<- 1:nrow(teinerohkemkui10lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 10
teine_OTUd_summa<-full_join(teine_OTUd_summa,teinerohkemkui10lugem)
#panen kolumnid j?rjekorda
teine_OTUd_summa<-teine_OTUd_summa[, c(1,2,3,4,5,7,6)]
### arvutan kui palju sekventse j??b proovides alles peale mitokondrite ja kloroplastide eemaldamist
teine_ilma_klor_mit<-grepl.sub(data=teine_sekveneerimine, pattern = c("Chloroplast","Mitochondria"), Var= "TAXON.ID", keep.found=F)
##leian palju sekventse j?i alles peale mitokondrite ja kloroplastide eemaldamist
teine_ilma_kolsum<-colSums(teine_ilma_klor_mit[-(1:2)])
teine_ilma_kolsum<-as.data.frame(teine_ilma_kolsum)
teine_ilma_kolsum$variable<-rownames(teine_ilma_kolsum)
teine_OTUd_summa<-full_join(teine_OTUd_summa,teine_ilma_kolsum)
#kirjutan teise kokkuv?tva tabeli excelisse
writeWorksheetToFile("Teine_sekv_alganalyys.xlsx", data=teine_OTUd_summa, sheet="tabel", startRow=1, startCol=1)
##k6ik esindatud h6imkonnad teises sekveneerimises
kogu_ID<-teine_ilma_klor_mit%>%separate(TAXON.ID,c("riik", "hmk", "klass", "selts", "sgk", "prk", "liik"), sep=";", extra="merge")
kogu_ID$hmk<-as.factor(kogu_ID$hmk)
levels(kogu_ID$hmk)
##[1] "Acidobacteria"       "Actinobacteria"      "Armatimonadetes"     "Bacteroidetes"      
##[5] "Chlamydiae"          "Chloroflexi"         "Deinococcus-Thermus" "Euryarchaeota"      
##[9] "Firmicutes"          "Gemmatimonadetes"    "Nitrospirae"         "Planctomycetes"     
##[13] "Proteobacteria"      "Saccharibacteria"    "Tenericutes"         "Thaumarchaeota"     
##[17] "Verrucomicrobia" 
##17 h6imkonda

#teise sekveneerimise keskmine OTUde arv proovides
summary(teine_OTUd_summa$length0)
#k6ikide sekventside summa 
teine_sekv_read_SUM<-sum(teine_OTUd_summa$summa)
#kloroplastide lugemite summa
teineklor<-grepl.sub(data = teine_sekveneerimine, pattern = c("Chloroplast"), Var = "TAXON.ID", keep.found=T)
teineklor2<-teineklor[,3:52]
teineklor_sum<-sum(teineklor2)
#Arvutan mitokondrite lugemite summa
teine_mitokondrid<-grepl.sub(data = teine_sekveneerimine, pattern = c("Mitochondria"), Var = "TAXON.ID", keep.found=T)
teine_mitokondrid_2<-teine_mitokondrid[,3:52]
teine_mit_sum<-sum(teine_mitokondrid_2)
#arvutan, mitu OTU-t on kloroplastide seas
teineklor_OTUd<-length(teineklor$X)
#arvutan mitokondrite OTUd
teinemitokondrid_OTU<-length(teine_mitokondrid$X)
#arvutan kui palju on arhesid
teine_arhed<-grepl.sub(data = teine_sekveneerimine, pattern = c("Archaea"), Var = "TAXON.ID", keep.found=T)
teine_arhed<-teine_arhed[,3:52]
teine_arhed_SUM<-sum(teine_arhed)
#k6ik OTU-d
teine_k6ik_OTUd<-length(teine_sekveneerimine$X)
#ilma mitokondrite ja kloroplastideta OTU-d
teine_ilma_mitokon_kloroplastiteta_OTUd<-teine_k6ik_OTUd-teinemitokondrid_OTU-teineklor_OTUd
#ilma mitokondrite ja kloroplastiteta sekventsid
teine_sekv_lugemid_ilma_mit_klor_ta<-teine_sekv_read_SUM-teine_mit_sum-teineklor_sum

##Kogu sekveneerimise tulemuste anal??s, algtabel k6ikkokku
#eemaldan ainult nullidega OTUd
kogu_ridade_rsum<-rowSums(k6ikkokku[-(51:52)])
kogu_ridade_rsum<-as.data.frame(kogu_ridade_rsum)
#liidan ridade summade kolumni esimesele 
kogu_sekv_sum<-cbind(k6ikkokku, kogu_ridade_rsum)
kogu_sekv_sum<-kogu_sekv_sum%>%rename(summa=kogu_ridade_rsum)
#eemaldan nullidega read
kogu_sekv_0ta<-kogu_sekv_sum%>%filter(kogu_sekv_sum$summa>0)
#eemaldan kolumni summa
kogu_sekveneerimine<-select(kogu_sekv_0ta,-summa)
##leian proovide kohta keskmise lugemite arvu
kogu<-kogu_sekveneerimine[, 1:50]
b2<-colSums(kogu)%>%as.data.frame()
summary(b2)
#arvutan, mitu erinevat OTU-t on proovides, kus on lugemeid rohkem >0, >10 ja >100,K6IK ANDMED
kogu_OTUde_jaotus<-kogu_sekveneerimine%>%melt(id.vars="X", measure.vars=names(kogu_sekveneerimine[,-(51:52)]))%>%
  group_by(variable)%>%
  summarise(length0=length(value[value>0]), length10=length(value[value>10]), length100=length(value[value>100]))

#liidan lugemite summa kolumni OTUde jaotuse tabelile
#kasutan eelnevat tehtud kolumnisummade tabelit b2
#teen kolumni "katse"?
kogu_kolumn_sum<-b2%>%rownames_to_column(var = "katse")
kogu_kolumn_sum$summa<-kogu_kolumn_sum$.
kogu_kolumn_sum<-kogu_kolumn_sum%>%select(katse, summa)
###liidan kaks tabelit- OTUde arvud ja lugemite summa 
kogu_OTUd_summa<-full_join(kogu_OTUde_jaotus, kogu_kolumn_sum, by=c("variable" = "katse"))
kogu_OTUd_summa<-as.data.frame(kogu_OTUd_summa)

#lugemite arvutamised
#>100 lugemit
idx=kogu>100
kogurohkemkui100lugem<-colSums(kogu*idx)
kogurohkemkui100lugem<-kogurohkemkui100lugem%>%as.data.frame()
kogurohkemkui100lugem$sumrohkem100<-kogurohkemkui100lugem$.
kogurohkemkui100lugem<-select(kogurohkemkui100lugem,sumrohkem100)
#teen variable kolumni
kogurohkemkui100lugem$variable<- rownames(kogurohkemkui100lugem)
#panen numbrid rownamedeks
rownames(kogurohkemkui100lugem)<- 1:nrow(kogurohkemkui100lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 100
kogu_OTUd_summa<-full_join(kogu_OTUd_summa,kogurohkemkui100lugem)

#>10 lugemit
idx2=kogu>10
kogurohkemkui10lugem<-colSums(kogu*idx2)
kogurohkemkui10lugem<-kogurohkemkui10lugem%>%as.data.frame()
kogurohkemkui10lugem$sumrohkem10<-kogurohkemkui10lugem$.
kogurohkemkui10lugem<-select(kogurohkemkui10lugem,sumrohkem10)
#teen variable kolumni
kogurohkemkui10lugem$variable<- rownames(kogurohkemkui10lugem)
#panen numbrid rownamedeks
rownames(kogurohkemkui10lugem)<- 1:nrow(kogurohkemkui10lugem)
#lisan nende OTU-de sekventside summad, mida on rohkem kui 10
kogu_OTUd_summa<-full_join(kogu_OTUd_summa,kogurohkemkui10lugem)
#panen kolumnid j?rjekorda
kogu_OTUd_summa<-kogu_OTUd_summa[, c(1,2,3,4,5,7,6)]
### arvutan kui palju sekventse j??b proovides alles peale mitokondrite ja kloroplastide eemaldamist
kogu_ilma_klor_mit<-grepl.sub(data=kogu_sekveneerimine, pattern = c("Chloroplast","Mitochondria"), Var= "TAXON.ID", keep.found=F)
##leian palju sekventse j?i alles peale mitokondrite ja kloroplastide eemaldamist
kogu_ilma_kolsum<-colSums(kogu_ilma_klor_mit[-(51:52)])
kogu_ilma_kolsum<-as.data.frame(kogu_ilma_kolsum)
kogu_ilma_kolsum$variable<-rownames(kogu_ilma_kolsum)
kogu_OTUd_summa<-full_join(kogu_OTUd_summa,kogu_ilma_kolsum)
#kirjutan teise kokkuv?tva tabeli excelisse
writeWorksheetToFile("Kogu_sekv_alganalyys.xlsx", data=kogu_OTUd_summa, sheet="tabel", startRow=1, startCol=1)
##k6ik esindatud h6imkonnad teises sekveneerimises
kogu_ID<-kogu_ilma_klor_mit%>%separate(TAXON.ID,c("riik", "hmk", "klass", "selts", "sgk", "prk", "liik"), sep=";", extra="merge")
kogu_ID$hmk<-as.factor(kogu_ID$hmk)
levels(kogu_ID$hmk)
###[1] "Acidobacteria"       "Actinobacteria"      "Armatimonadetes"     "Bacteroidetes"      
###[5] "Chlamydiae"          "Chloroflexi"         "Deinococcus-Thermus" "Euryarchaeota"      
###[9] "Firmicutes"          "Gemmatimonadetes"    "Gracilibacteria"     "Nitrospirae"        
###[13] "Planctomycetes"      "Proteobacteria"      "Saccharibacteria"    "Tenericutes"        
###[17] "Thaumarchaeota"      "TM6"                 "Verrucomicrobia"  
### 19 h6imkonda

#kogu sekveneerimise keskmine OTUde arv proovides
summary(kogu_OTUd_summa$length0)
#k6ikide sekventside summa 
kogu_sekv_read_SUM<-sum(kogu_OTUd_summa$summa)
#kloroplastide lugemite summa
koguklor<-grepl.sub(data = kogu_sekveneerimine, pattern = c("Chloroplast"), Var = "TAXON.ID", keep.found=T)
koguklor2<-koguklor[,1:50]
koguklor_sum<-sum(koguklor2)
#Arvutan mitokondrite lugemite summa
kogu_mitokondrid<-grepl.sub(data = kogu_sekveneerimine, pattern = c("Mitochondria"), Var = "TAXON.ID", keep.found=T)
kogu_mitokondrid_2<-kogu_mitokondrid[,1:50]
kogu_mit_sum<-sum(kogu_mitokondrid_2)
#arvutan, mitu OTU-t on kloroplastide seas
koguklor_OTUd<-length(koguklor$X)
#arvutan mitokondrite OTUd
kogumitokondrid_OTU<-length(kogu_mitokondrid$X)
#arvutan kui palju on arhesid
kogu_arhed<-grepl.sub(data = kogu_sekveneerimine, pattern = c("Archaea"), Var = "TAXON.ID", keep.found=T)
kogu_arhed<-kogu_arhed[,1:50]
kogu_arhed_SUM<-sum(kogu_arhed)
#k6ik OTU-d
kogu_k6ik_OTUd<-length(kogu_sekveneerimine$X)
#ilma mitokondrite ja kloroplastideta OTU-d
kogu_ilma_mitokon_kloroplastiteta_OTUd<-kogu_k6ik_OTUd-kogumitokondrid_OTU-koguklor_OTUd
#ilma mitokondrite ja kloroplastiteta sekventsid
kogu_sekv_lugemid_ilma_mit_klor_ta<-kogu_sekv_read_SUM-kogu_mit_sum-koguklor_sum


###Algandmete kvaliteet.Lisa 8
library(tidyverse)
library(reshape2)
library(ggplot2)
df<-read.csv2("algne.csv")

df_mito<-df %>% filter(grepl("Mitochon", TAXON.ID))
df_chloro<- df %>% filter(grepl("Chloroplast", TAXON.ID))

df_mito1<-df_mito %>% arrange(desc(kark1)) 
df_mito1<- df_mito1[1:3,]
df_mito1<-melt(df_mito1)
ggplot(df_mito1)+ geom_col(aes(variable, log10(value), fill=X))+
  theme(axis.text.x = element_text(size=7, angle=90))+
  xlab("Sekveneerimisproovid")+
  labs(fill = " ")

### I ja II sekveneerimisandmete Shannon, Lisa 4
y1<-read.csv2("algne.csv") %>% na.omit()
##eemaldan t?hjad read 
y11<-y1[!apply(y1 == "", 1, all),]
y12<-y11[complete.cases(y11),]
##eemaldan mitokondrid ja kloroplastid
y2<-grepl.sub(data=y12, pattern = c("Chloroplast","Mitochondria"), Var= "TAXON.ID", keep.found=F)
##eemaldan t?hjad read
y2_sum<-rowSums(y2[-(1:2)])
y2_sum<-as.data.frame(y2_sum)
#liidan ridade summade kolumni esimese sekveneerimise tabelile 
y2_summaga<-cbind(y2, y2_sum)
names(y2_summaga)
y2_summaga<-y2_summaga%>%rename(summa=y2_sum)
#eemaldan nullidega read
y3<-y2_summaga%>%filter(y2_summaga$summa>0)
#eemaldan kolumni summa veeru, saan tabeli, kus ei ole nulli ridasid
y3<-select(y3,-summa)
##muudan nimesid
names(y3)
y3<-y3%>%rename(OTU=X)
b1<-y3%>%select(1:93, TOPK=num_range("TK", 1:5), TOPS=num_range("TS", 1:5))
b1<-b1%>%select(OTU:KARS5, KALK=num_range("KLK", 1:5), KALS=num_range("KLS", 1:5), PEK1:TOPS5)
names(b1)
library(tidyr)
b2<-b1%>%gather(exp, count, 3:ncol(b1))# dim 88300 X 4 long format
b2$exp<-as.factor(b2$exp)
levels(b2$exp)
b2$count[b2$count<3]<- 0
b2[ is.na(b2) ] <- 0 

b3<-b2%>%select(exp, OTU, count)
b3<-b3%>%spread(OTU, count)
b3<-as_data_frame(b3)
b3$exp<-as.character(b3$exp)
b3<- b3%>%mutate(exp1=exp)
names(b3)

b4<-b3%>%separate(exp, into=c("a", "b"), sep= -2)
b4<- b4%>%separate(a, into=c("a", "location"), sep= -2)
names(b4)
b4<-b4%>%select(-a, -b, location, 4:820)
names(b4)
b4<-b4%>%select(exp1, 1:818)
names(b4)
b4$location<- tolower(b4$location)
##teen maatriksi
y100.matrix<-b4[,-(1:2)]
library(vegan)
shannon_y100_OTU<-diversity(y100.matrix, index="shannon", MARGIN=1)%>%as.data.frame()
shannon_y100_OTU<-cbind(ex=b4$exp1, loc=b4$location, shannon_y100_OTU)
names(shannon_y100_OTU)[3]<-paste("indeks") #annab veerule uue nime
#summeerimata andmed
library(ggplot2)
###Lisa I 
sh100<-shannon_y100_OTU
ggplot(data=sh100, aes(x=factor(ex), y=indeks))+
  geom_boxplot(data=sh100[c(1,3,5,7,8,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99),], colour="black")+
  geom_boxplot(data=sh100[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38.40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100),], colour="red")+
  theme_bw()+
  coord_flip()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(size = 12, face = "bold"),axis.text.y = element_text(size = 5))+
  xlab("PROOVID") + ylab("Shannon-Weaver indeks")
##töös olevate joonised
ggplot(data=shannon_y100_OTU[1:50,], aes(x=factor(ex), y=indeks))+
  geom_boxplot(data=shannon_y100_OTU[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49),],colour="black")+
  geom_boxplot(data=shannon_y100_OTU[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38.40,42,44,46,48,50),],colour="red")+
  theme_bw()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(size = 12, face = "bold"),axis.text.y = element_text(size = 7, face = "bold"))+
  coord_flip()+
  xlab("PROOVID") + ylab("Shannon-Weaver indeks")
ggplot(data=shannon_y100_OTU[51:100,], aes(x=factor(ex), y=indeks))+
  geom_boxplot(data=shannon_y100_OTU[c(51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99),],colour="black")+
  geom_boxplot(data=shannon_y100_OTU[c(50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100),], colour="red")+
  theme_bw()+
  theme(axis.title = element_text(size = 12),axis.text.x = element_text(size = 12, face = "bold"),axis.text.y = element_text(size = 7, face = "bold"))+
  coord_flip()+
  xlab("PROOVID") + ylab("Shannon-Weaver indeks")



###Venni diagramm Lisa 7
##teen esimese sekveneerimise tabeli
esimene<-b1[,1:52]
esimene_sum<-rowSums(esimene[(3:52)])
esimene_sum<-as.data.frame(esimene_sum)
#liidan ridade summade kolumni esimese sekveneerimise tabelile 
esimene_sum<-cbind(esimene, esimene_sum)
esimene_sum<-esimene_sum%>%rename(summa=esimene_sum)
#eemaldan nullidega read
esimene_OTU<-esimene_sum%>%filter(esimene_sum$summa>0)
#eemaldan kolumni summa veeru, saan tabeli, kus ei ole nulli ridasid
esimene_OTU<-select(esimene_OTU,-summa)
##709 OTUt
##teine tabel
teine<-subset(b1, select =c (1,2,53:102))
teine_sum<-rowSums(teine[(3:52)])
teine_sum<-as.data.frame(teine_sum)
#liidan ridade summade kolumni esimese sekveneerimise tabelile 
teine_sum<-cbind(teine, teine_sum)
teine_sum<-teine_sum%>%rename(summa=teine_sum)
#eemaldan nullidega read
teine_OTU<-teine_sum%>%filter(teine_sum$summa>0)
#eemaldan kolumni summa veeru, saan tabeli, kus ei ole nulli ridasid
teine_OTU<-select(teine_OTU,-summa)
##636 OTUt
##leian OTUde arvu, mis on m6lemas tabelis
length(intersect(esimene_OTU[[1]], teine_OTU[[1]]))
##yhiseid OTUsid on 529
##teen Venn diagrammi
library(VennDiagram)
require(VennDiagram)
venn.diagram(list(Sekvn1 = 1:1238, Sekvn2 = 710:1345),fill = c("red", "green"),
             alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, 
             filename = "venn.emf");

#################################################################
###Diagrammid###tabeliks b1
#################################################################
names(b1) #[1] "OTU" "TAXON.ID" "kark1"  "kark2"  "kark3" "kark4" jne
#teen teise eksperimendi kolumni, kus on ainult juurviljade sisu ja koor
y2<-b1%>%gather(exp, counts, 3:ncol(y1)) # dim 88300 X 4 long format
names(y2) #[1] "OTU" "TAXON.ID" "exp" "counts"
y2<-y2%>%mutate(exp1=exp)%>%separate(exp1, into = c("exp1", "n1"), sep = -2)%>%select(-n1)
y2$exp1<- as.factor(y2$exp1)
levels(y2$exp1)
#[1] "KARK" "KARS" "KLK"  "KLS"  "PEK"  "PES"  "PORK" "PORS" "TK"   "TS"   "kalk" "kals" "kark"
#[14] "kars" "pek"  "pes"  "pork" "pors" "topk" "tops"
#Normaliseerimine, viskasin v?lja OTUd, kus oli v?hem kui 20 counti
y3<-y2%>%filter(counts>20)%>% 
  group_by(exp)%>%mutate(n.counts=(counts/sum(counts))* 20)
#kontrollin, kas k?ikide runide normeritud summad on v?rdsed
ykontroll<-y3%>%group_by(exp)%>%summarise(SUM=sum(n.counts))

#muudan nimesid exp kolumnis
y3$exp<- as.character(y3$exp)
y3$exp<- tolower(y3$exp)
y3$exp<- as.factor(y3$exp)
levels(y3$exp) #50 erinevat

#muudan nimesid exp1 kolumnis
y3$exp1<- as.character(y3$exp1)
y3$exp1<- tolower(y3$exp1)
y3$exp1<- as.factor(y3$exp1)
levels(y3$exp1)#10 erinevat
##leian kui palju erinevaid OTUsid j?i peale normaliseerimist alles
length(unique(y3$OTU))
#247
##25+25 andmekogu
#50+50 proovi summeerimine, same 25+25 andmekogu
y50<-y3%>%group_by(exp, OTU, TAXON.ID, exp1)%>% summarise(counts.summed= sum(n.counts), Nr.exps= n())
#eraldan TAXON.ID mitmeks kolumniks
y50_1<-y50%>%separate(TAXON.ID,c("riik", "hmk", "klass", "selts", "sgk", "prk", "liik"), sep=';', extra="merge")

###2% OTU diagrammid
#v?tan v?lja k?ik OTUd, mida on rohkem kui kolm protsenti lugemitest
y50_2<-y50%>%filter(counts.summed>2)
#v?tan v?lja k?ik OTUd, mida on v?hem kui kolm protsenti lugemitest
y50_2_2<-y50%>%filter(counts.summed<2)
##arvutan mitu protsenti kuulub OTUde alla, mida on alla 3 protsendi
y50_2_3<-y50_2_2%>%group_by(exp)%>%summarise(counts.summed=sum(counts.summed))
#teen uue kolumni OTU, kus k?ikides ridades on "Ja teised"
y50_2_4<-y50_2_3%>%mutate(OTU="Ja teised")
#kustutan kolumnid, mida mul vaja ei l?he
y50_3<-y50_2%>% select( exp, counts.summed, OTU)
#liidan kaks tabelit
y51<-bind_rows(y50_3, y50_2_4)
#grupeerin katsete kaupa
y51_1<-y51%>%group_by(exp)

####2%  diagramm
##KOOR 25 2%
palette2_koor=c("lemonchiffon3", "deepskyblue", "red", "red3", "#FFC300", "cadetblue1", "yellow", "darkorange1", "indianred1","lightgoldenrod1","black","mediumpurple", "dodgerblue","cornsilk1", "darkorange3", "purple", "turquoise2", "maroon")

y25_koor<-y51%>%filter(exp %in% c("kalk1", "kalk2", "kalk3", "kalk4", "kalk5", "kark1", "kark2", "kark3", "kark4", "kark5", "pek1", "pek2", "pek3", "pek4", "pek5", "pork1", "pork2", "pork3", "pork4", "pork5", "topk1", "topk2", "topk3", "topk4", "topk5"))
y25_koor_2perc<-y25_koor%>%filter(counts.summed>2)%>%ggplot(aes(x=exp, y=counts.summed)) +
  geom_bar(position= "fill", stat= "identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_fill_manual(values = palette2_koor, breaks=c("OTU1","OTU3","OTU14", "OTU52","OTU11", "OTU12", "OTU181", "OTU79","OTU21", "OTU44", "OTU201", "OTU13", "OTU18", "OTU39", "OTU15", "OTU20", "OTU35", "Ja teised" ),
                    labels=c("OTU1:Enterobacteriaceae", "OTU3:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU52:Luiteimonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei", "OTU181:Sphingomonas sp.", "OTU79:Sphingopyxis sp.", "OTU21:Janthinobacterium sp.", "OTU44:Comomonadaceae", "OTU201:Sandaracinaceae", "OTU13:Arthrobacter sp.", "OTU18:Microbacterium sp.", "OTU39:Streptomyces turgidiscabies", "OTU15:Pedobacter sp.", "OTU20:Flavobacterium sp.", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
y25_koor_2perc

##SISU 25 2%
palette2_sisu=c("lemonchiffon3", "deepskyblue", "olivedrab", "red", "red3", "#FFC300", "cadetblue1", "yellow", "magenta", "mediumpurple", "dodgerblue", "paleturquoise3", "cornsilk1", "orange", "skyblue1", "plum", "olivedrab3","olivedrab1", "mediumblue")

y25_sisu<-y51%>%filter(exp %in% c("kals1", "kals2", "kals3", "kals4", "kals5", "kars1", "kars2", "kars3", "kars4", "kars5", "pes1", "pes2", "pes3", "pes4", "pes5", "pors1", "pors2", "pors3", "pors4", "pors5", "tops1", "tops2", "tops3", "tops4", "tops5"))
y25_sisu_2perc<-y25_sisu%>%filter(counts.summed>2)%>%ggplot(aes(x=exp, y=counts.summed)) +
  geom_bar(position= "fill", stat= "identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_fill_manual(values = palette2_sisu, breaks=c("OTU1","OTU45","OTU34","OTU3","OTU8","OTU14", "OTU11", "OTU12", "OTU21", "OTU49", "OTU13", "OTU38", "OTU6", "OTU5", "OTU10", "OTU15", "OTU16", "OTU35", "Ja teised"),
                    labels=c("OTU1:Enterobacteriaceae", "OTU45:Enterobacteriaceae","OTU34:Pantoea sp.", "OTU3:Pseudomonas sp.", "OTU8:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei", "OTU21:Janthinobacterium sp.", "OTU49:Variovorax paradoxus", "OTU13:Arthrobacter sp.", "OTU38:Arthrobacter sp.", "OTU6:Bacillus subtilis.", "OTU5:Paenibacillus polymyxa.", "OTU10:Leuconostoc gelidum", "OTU15:Pedobacter sp.", "OTU16:Saccharibacteria", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
y25_sisu_2perc

##SISU ja KOOR koos
palette2=c("lemonchiffon3", "deepskyblue", "olivedrab", "red", "red3", "#FFC300", "cadetblue1", "yellow", "magenta", "darkorange1", "indianred1", "lightgoldenrod1", "black", "mediumpurple", "dodgerblue", "paleturquoise3", "cornsilk1", "orange", "darkorange3", "purple", "skyblue1", "plum", "olivedrab3", "turquoise2", "olivedrab1", "maroon", "mediumblue")

y50_2perc<-y51%>%filter(counts.summed>2)%>%ggplot(aes(x=exp, y=counts.summed)) +
  geom_bar(position= "fill", stat= "identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_fill_manual(values = palette2, breaks=c("OTU1","OTU45", "OTU34", "OTU3","OTU8", "OTU14", "OTU52","OTU11", "OTU12", "OTU181", "OTU79","OTU21", "OTU44", "OTU49", "OTU201", "OTU13", "OTU38", "OTU18", "OTU39", "OTU6", "OTU5", "OTU10", "OTU15", "OTU20", "OTU16", "OTU35", "Ja teised" ),
                    labels=c("OTU1:Enterobacteriaceae", "OTU45:Enterobacteriaceae","OTU34:Pantoea sp.", "OTU3:Pseudomonas sp.", "OTU8:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU52:Luiteimonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei", "OTU181:Sphingomonas sp.", "OTU79:Sphingopyxis sp.", "OTU21:Janthinobacterium sp.", "OTU44:Comomonadaceae", "OTU49:Variovorax paradoxus", "OTU201:Sandaracinaceae", "OTU13:Arthrobacter sp.", "OTU38:Arthrobacter sp.", "OTU18:Microbacterium sp.", "OTU39:Streptomyces turgidiscabies", "OTU6:Bacillus subtilis.", "OTU5:Paenibacillus polymyxa.", "OTU10:Leuconostoc gelidum", "OTU15:Pedobacter sp.", "OTU20:Flavobacterium sp.", "OTU16:Saccharibacteria", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.position="bottom", legend.key.size=unit(.7, "cm"), legend.text = element_text(size = 10))
y50_2perc

##andmekogu 10 OTU diagramm
#yhendan 25+25 proovi, saan 5+5 andmekogu
y10<-y50%>% group_by(exp1, OTU, TAXON.ID)%>% summarise(counts = sum(counts.summed), Nr.exps= n())
#eraldan TAXON.ID mitmeks kolumniks
y10<-y10%>%separate(TAXON.ID,c("riik", "hmk", "klass", "selts", "sgk", "prk", "liik"), sep=';', extra="merge")
#panen eraldi tabelisse k?ik OTUd, mida on rohkem kui kolm protsenti lugemitsest
y10_2<-y10%>%filter(counts>2)
#panen eraldi tabelisse k?ik OTUd, mida on v?hem kui kolm protsenti lugemitest
y10_2_2<-y10%>%filter(counts<2)
##arvutan mitu protsenti kuulub OTUde alla, mida on alla 2%
y10_2_3<-y10_2_2%>%group_by(exp1)%>%summarise(counts=sum(counts))
#teen uue kolumni OTU, kus k?ikides ridades on "Ja teised"
y10_2_4<-y10_2_3%>%mutate(OTU="Ja teised")
#kustutan kolumnid, mida mul vaja ei l?he
y10_3<-y10_2%>% select(exp1, counts, OTU)
#liidan kaks tabelit
y11<-bind_rows(y10_3, y10_2_4)
#grupeerin katsete kaupa
y11_1<-y11%>%group_by(exp1)

##### 2% OTU diagrammid
palette2_10=c("lemonchiffon3", "deepskyblue", "olivedrab", "red", "red3","coral", "#FFC300", "cadetblue1", "yellow","peru","magenta", "darkorange1", "indianred1","lightcyan", "mediumslateblue","lightgoldenrod1", "black", "mediumpurple","lightsalmon","red4","burlywood4", "palevioletred", "darkorange4", "dodgerblue","brown1" ,"paleturquoise3", "cornsilk1","firebrick" ,"orange", "darkorange3","sandybrown", "pink2","sienna2", "purple", "skyblue1", "orchid1" ,"olivedrab3", "turquoise2","chocolate1", "olivedrab1", "maroon", "mediumblue", "steelblue4","pink3","lightgreen","seagreen","gold")

y10_2perc<-y11_1%>%filter(counts>2)%>%ggplot(aes(x=exp1, y=counts))+
  geom_bar(position="fill", stat="identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette2_10,breaks=c("OTU1","OTU45", "OTU19", "OTU34", "OTU3","OTU8","OTU83", "OTU14", "OTU52","OTU11", "OTU12","OTU41","OTU27", "OTU181", "OTU79","OTU24","OTU32","OTU37","OTU88","OTU21", "OTU44","OTU199","OTU49", "OTU201", "OTU13", "OTU38", "OTU18", "OTU39","OTU25", "OTU56", "OTU159", "OTU29","OTU124","OTU23","OTU25","OTU40","OTU42", "OTU6", "OTU5", "OTU10","OTU9","OTU90", "OTU15", "OTU20","OTU94","OTU16", "OTU35", "Ja teised" ),
                    labels=c("OTU1:Enterobacteriaceae", "OTU45:Enterobacteriaceae", "OTU19:Citrobacter sp.","OTU34:Pantoea sp.", "OTU3:Pseudomonas sp.", "OTU8:Pseudomonas sp.","OTU83:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU52:Luiteimonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei","OTU41:Bradyrhizobium liaoningense","OTU27:Methylobacterium sp.", "OTU181:Sphingomonas sp.", "OTU79:Sphingopyxis sp.", "OTU24:Devosia sp.", "OTU32:Phyllobacterium sp.","OTU37: Rhodopseudomonas sp.", "OTU88:Aureimonas sp.", "OTU21:Janthinobacterium sp.", "OTU44:Comomonadaceae","OTU199:Comamonadaceae","OTU49:Variovorax paradoxus", "OTU201:Sandaracinaceae", "OTU13:Arthrobacter sp.", "OTU38:Arthrobacter sp.", "OTU18:Microbacterium sp.", "OTU39:Streptomyces turgidiscabies", "OTU25:Nocardioides sp.", "OTU56:Nocardioides sp.", "OTU159:Nocardioides sp.", "OTU29:Aeromicrobium ginsengisoli", "OTU124:Nocardioides sp.","OTU23:Sanguibacter inulinus","OTU25:Nocardioides sp.","OTU40:Promicromonospora sp.","OTU42:Rhodococcus sp.", "OTU6:Bacillus sp.", "OTU5:Paenibacillus sp.", "OTU10:Leuconostoc gelidum", "OTU9:Enterococcus faecium","OTU90: Bacillus pocheonensis","OTU15:Pedobacter sp.", "OTU20:Flavobacterium sp.","OTU94:Chryseobacterium sp.","OTU16:Saccharibacteria", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill= " ")+
  theme(legend.position="bottom",legend.key.size=unit(.7, "cm"),legend.text = element_text(size = 8,face = "italic"))
y10_2perc

##KOOR 2% 5 juurvilja
palette2_5_KOOR=c("lemonchiffon3", "deepskyblue","red", "red3","coral", "#FFC300", "cadetblue1", "yellow","peru","darkorange1", "indianred1","lightcyan", "mediumslateblue","lightgoldenrod1","black", "mediumpurple","lightsalmon","red4","burlywood4", "palevioletred", "darkorange4", "dodgerblue","brown1" ,"paleturquoise3", "cornsilk1","firebrick","darkorange3","sandybrown", "pink2","sienna2", "purple","turquoise2","chocolate1", "maroon", "mediumblue","pink3","seagreen","gold")

y11_koor<-y11%>%filter(exp1 %in% c("kalk", "kark", "pek", "pork", "topk"))
y5_koor_2perc<-y11_koor%>%filter(counts>2)%>%ggplot(aes(x=exp1, y=counts))+
  geom_bar(position="fill", stat="identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette2_5_KOOR,breaks=c("OTU1","OTU19", "OTU34", "OTU3","OTU8", "OTU14", "OTU52","OTU11", "OTU12","OTU41","OTU27", "OTU181", "OTU79","OTU24","OTU32","OTU37","OTU88","OTU21", "OTU44","OTU199", "OTU201", "OTU13", "OTU18", "OTU39","OTU25", "OTU56", "OTU159", "OTU29","OTU124","OTU23","OTU25","OTU40","OTU42","OTU90", "OTU15", "OTU20","OTU94", "OTU35", "Ja teised" ),
                    labels=c("OTU1:Enterobacteriaceae", "OTU19:Citrobacter sp.","OTU34:Pantoea sp.", "OTU3:Pseudomonas sp.", "OTU8:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU52:Luiteimonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei","OTU41:Bradyrhizobium liaoningense","OTU27:Methylobacterium sp.", "OTU181:Sphingomonas sp.", "OTU79:Sphingopyxis sp.", "OTU24:Devosia sp.", "OTU32:Phyllobacterium sp.","OTU37: Rhodopseudomonas sp.", "OTU88:Aureimonas sp.", "OTU21:Janthinobacterium sp.", "OTU44:Comomonadaceae","OTU199:Comamonadaceae", "OTU201:Sandaracinaceae", "OTU13:Arthrobacter sp.", "OTU18:Microbacterium sp.", "OTU39:Streptomyces turgidiscabies", "OTU25:Nocardioides sp.", "OTU56:Nocardioides sp.", "OTU159:Nocardioides", "OTU29:Aeromicrobium ginsengisoli", "OTU124:Nocardioides sp.","OTU23:Sanguibacter inulinus","OTU25:Nocardioides sp.","OTU40:Promicromonospora sp.","OTU42:Rhodococcus sp.", "OTU90: Bacillus pocheonensis","OTU15:Pedobacter sp.", "OTU20:Flavobacterium sp.","OTU94:Chryseobacterium", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill= " ")+
  theme(legend.key.size=unit(.7, "cm"),axis.text.x = element_text(size = 14, face = "bold" ))
y5_koor_2perc


###SISU 2% 5 juurvilja
palette2_5_SISU=c("lemonchiffon3", "deepskyblue", "olivedrab", "red", "red3", "#FFC300", "cadetblue1", "yellow","magenta", "mediumpurple","dodgerblue","paleturquoise3","cornsilk1", "orange", "skyblue1", "orchid1" ,"olivedrab3", "olivedrab1", "mediumblue", "steelblue4","lightgreen")
y11_sisu_2<-y11%>%filter(exp1 %in% c("kals", "kars", "pes", "pors", "tops"))

y5_sisu_2perc<-y11_sisu_2%>%filter(counts>2)%>%ggplot(aes(x=exp1, y=counts))+
  geom_bar(position="fill", stat="identity", width = 0.7, aes(fill=OTU), colour="black")+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette2_5_SISU, breaks=c("OTU1","OTU45", "OTU34", "OTU3","OTU8","OTU83", "OTU14","OTU11", "OTU12","OTU21", "OTU49","OTU13","OTU38","OTU6","OTU5","OTU10","OTU9","OTU15","OTU16","OTU35", "Ja teised" ),
                    labels=c("OTU1:Enterobacteriaceae", "OTU45:Enterobacteriaceae","OTU34:Pantoea sp.", "OTU3:Pseudomonas sp.", "OTU8:Pseudomonas sp.","OTU83:Pseudomonas sp.", "OTU14:Stenotrophomonas sp.", "OTU11:Rhizobium sp.", "OTU12:Rhizobium larrymoorei", "OTU21:Janthinobacterium sp.","OTU49:Variovorax paradoxus","OTU13:Arthrobacter sp.", "OTU38:Arthrobacter sp.","OTU6:Bacillus sp.", "OTU5:Paenibacillus sp.", "OTU10:Leuconostoc gelidum", "OTU9:Enterococcus faecium","OTU15:Pedobacter sp.","OTU16:Saccharibacteria", "OTU35:Phytoplasma sp.", "teised OTUd"))+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill= " ")+
  theme(legend.key.size=unit(.7, "cm"),axis.text.x = element_text(size = 14, face = "bold" ))
y5_sisu_2perc


#K?rgema taksoni p?hine tabel
#tabel 25+25 andmekogu kohta
#Proteobacteria tabel, mis tuleb klassi j?rgi
proteobakter<-y50_1%>%filter(hmk=="Proteobacteria")
levels(factor(proteobakter$klass))#vaatan, mis faktorid on kolumnis klass
proteobakter$klass<-as.factor(proteobakter$klass)
proteobakter.klass<-proteobakter%>%group_by(exp, klass)%>%summarise(summa=sum(counts.summed))

#Firmicutese tabel, mis tuleb klassi j?rgi
firmicutes<-y50_1%>%filter(hmk=="Firmicutes")
firmicutes$klass<-as.factor(firmicutes$klass)
firmicutes.klass<-firmicutes%>%group_by(exp,klass)%>%summarise(summa=sum(counts.summed))

#teised h6imkonnad - Viskan Proteobacteria ja Firmicutese tabelist v?lja, h?imkonna tasemel
rest<- y50_1%>%filter(hmk != "Proteobacteria", hmk!= "Firmicutes")
rest$hmk<-as.factor(rest$hmk)
rest.hmk<-rest%>%group_by(exp,hmk)%>%summarise(summa=sum(counts.summed))%>%rename(klass=hmk)

#panen tabelid ?ksteise otsa kokku
y50.class<-rbind(proteobakter.klass, firmicutes.klass, rest.hmk)%>%rename(counts=summa)
#v?tan v?lja tabeli, kus on rohkem kui 1%
y50.class_1<-y50.class%>%filter(counts>1)
#v?tan v?lja tabeli, kus on v?hem kui 1%
y50.class_2<-y50.class%>%filter(counts<1)
y50.class_3<-y50.class_2%>%group_by(exp)%>%summarise(counts=sum(counts))
#teen uue kolumni OTU, kus k?ikides ridades on "Ja teised"
y50.class_3<-y50.class_3%>%mutate(klass="Ja teised")
#liidan kaks tabelit
y50.class<-bind_rows(y50.class_3, y50.class_1)
#grupeerin katsete kaupa
y51.class<-y50.class%>%group_by(exp)


##diagramm 25+25 andmekogu
palette=c("#FFC300", "#DC143C", "olivedrab3", "yellow", "purple", "black", "deepskyblue1","lemonchiffon3", "magenta", "cornsilk")

##KOORE diagramm 25
y25.class_KOOR1<-y50.class%>%filter(exp %in% c("kalk1", "kalk2", "kalk3", "kalk4", "kalk5", "kark1", "kark2", "kark3", "kark4", "kark5", "pek1", "pek2", "pek3", "pek4", "pek5", "pork1", "pork2", "pork3", "pork4", "pork5", "topk1", "topk2", "topk3", "topk4", "topk5"))
d25hmk_1perc_KOOR<-y25.class_KOOR1%>%filter(counts>1)%>%ggplot(aes(x=exp, y=counts))+
  geom_bar(stat="identity",width = 0.7, aes(fill=klass), position = "fill", colour="black")+
  scale_fill_manual(values = palette, breaks=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "Ja teised"),
                    labels=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "teised kõrgemad taksonid"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"),axis.text.x = element_text(size = 14, face = "bold" ))
d25hmk_1perc_KOOR


###SISU diagramm 25
palette_SISU=c("#FFC300", "#DC143C", "olivedrab3", "yellow", "purple", "deepskyblue1","lemonchiffon3", "magenta", "cornsilk")

y25.class_SISU<-y50.class%>%filter(exp %in% c("kals1", "kals2", "kals3", "kals4", "kals5", "kars1", "kars2", "kars3", "kars4", "kars5", "pes1", "pes2", "pes3", "pes4", "pes5", "pors1", "pors2", "pors3", "pors4", "pors5", "tops1", "tops2", "tops3", "tops4", "tops5"))
d25hmk_1perc_SISU<-y25.class_SISU%>%filter(counts>1)%>%ggplot(aes(x=exp, y=counts))+
  geom_bar(stat="identity",width = 0.7, aes(fill=klass), position = "fill", colour="black")+
  scale_fill_manual(values = palette_SISU,breaks=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "Ja teised"),
                    labels=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "teised k?rgemad taksonid"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
d25hmk_1perc_SISU

#andmekogu 10
proteobakter<-y10%>%filter(hmk=="Proteobacteria")
levels(factor(proteobakter$klass))#vaatan, mis faktorid on kolumnis klass
proteobakter$klass<-as.factor(proteobakter$klass)
proteobakter.klass<-proteobakter%>%group_by(exp1, klass)%>%summarise(summa=sum(counts))

#Firmicutese tabel, mis tuleb klassi j?rgi
firmicutes<-y10%>%filter(hmk=="Firmicutes")
firmicutes$klass<-as.factor(firmicutes$klass)
firmicutes.klass<-firmicutes%>%group_by(exp1,klass)%>%summarise(summa=sum(counts))

#teised h6imkonnad - Viskan Proteobacteria ja Firmicutese tabelist v?lja, h?imkonna tasemel
rest<- y10%>%filter(hmk != "Proteobacteria", hmk!= "Firmicutes")
rest$hmk<-as.factor(rest$hmk)
rest.hmk<-rest%>%group_by(exp1,hmk)%>%summarise(summa=sum(counts))%>%rename(klass=hmk)
#panen tabelid ?ksteise otsa kokku
y10.class<-rbind(proteobakter.klass, firmicutes.klass, rest.hmk)%>%rename(counts=summa)
#v?tan v?lja tabeli, kus on rohkem kui 1%
y10.class_1<-y10.class%>%filter(counts>1)
#v?tan v?lja tabeli, kus on v?hem kui 1%
y10.class_2<-y10.class%>%filter(counts<1)
y10.class_3<-y10.class_2%>%group_by(exp1)%>%summarise(counts=sum(counts))
#teen uue kolumni OTU, kus k?ikides ridades on "Ja teised"
y10.class_3<-y10.class_3%>%mutate(klass="Ja teised")
#liidan kaks tabelit
y10.class<-bind_rows(y10.class_3, y10.class_1)
#grupeerin katsete kaupa
y10.class<-y10.class%>%group_by(exp1)
#diagramm k6ik 5+5
d10hmk_1perc<-y10.class%>%filter(counts>1)%>%ggplot(aes(x=exp1, y=counts, fill=klass)) +
  geom_bar(position= "fill", stat= "identity",width = 0.7, aes(fill=klass), colour="black")+
  scale_fill_manual(values = palette, breaks=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "Ja teised"),
                    labels=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Beetaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "teised kõrgemad taksonid"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
d10hmk_1perc

##diagramm 5 KOOR
palette_5_KOOR=c("#FFC300", "#DC143C", "olivedrab3", "yellow", "purple", "black", "deepskyblue1", "magenta", "cornsilk")

y5.class_KOOR<-y10.class%>%filter(exp1 %in% c("kalk", "kark", "pek", "pork", "topk"))
d5hmk_1perc_KOOR<-y5.class_KOOR%>%filter(counts>1)%>%ggplot(aes(x=exp1, y=counts, fill=klass)) +
  geom_bar(position= "fill", stat= "identity",width = 0.7, aes(fill=klass), colour="black")+
  scale_fill_manual(values = palette_5_KOOR)+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
d5hmk_1perc_KOOR

####diagramm 5 SISU
palette_5_SISU=c("#FFC300", "#DC143C", "olivedrab3", "yellow", "purple", "deepskyblue1","lemonchiffon3", "magenta", "cornsilk")

y5.class_SISU<-y10.class%>%filter(exp1 %in% c("kals", "kars", "pes", "pors", "tops"))
d5hmk_1perc_SISU<-y5.class_SISU%>%filter(counts>1)%>%ggplot(aes(x=exp1, y=counts, fill=klass)) +
  geom_bar(position= "fill", stat= "identity",width = 0.7, aes(fill=klass), colour="black")+
  scale_fill_manual(values = palette_5_SISU, breaks=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Betaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "Ja teised"),
                    labels=c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidetes", "Beetaproteobacteria", "Gammaproteobacteria", "Saccharibacteria", "Tenericutes", "teised k?rgemad taksonid"))+
  scale_y_continuous(labels = percent_format())+
  coord_flip()+
  xlab(" ")+
  ylab(" ")+
  labs(fill = " ")+
  theme(legend.key.size=unit(.7, "cm"), axis.text.x = element_text(size = 14, face = "bold" ))
d5hmk_1perc_SISU


#############shannoni indeksid
library(vegan)
#andmekogu50
y50.11<-y50%>%select(exp, OTU, counts.summed)%>%dcast(exp~OTU)
y50.11$exp<-as.character(y50.11$exp)
y50.11[is.na(y50.11)]<-0
y50.11<- y50.11%>%separate(exp, into=c("a", "b"), sep= -2)
y50.11<- y50.11%>%separate(a, into=c("a", "location"), sep= -2)
y50.12<-y50.11[,-(1:3)]
y50.12<-cbind(ex=y50.11$a, loc=y50.11$location, y50.12)
y50.matrix<-y50.12[,-(1:2)]
shannon_y50_OTU<-diversity(y50.matrix, index="shannon", MARGIN=1)%>%as.data.frame()
shannon_y50_OTU<-cbind(ex=y50.12$ex, loc=y50.12$loc, shannon_y50_OTU)
names(shannon_y50_OTU)[3]<-paste("indeks") #annab veerule uue nime

#varieeruvuse-erinevate juurviljade koor ja sisu
sh50otu<-ggplot(data=shannon_y50_OTU, aes(x=factor(ex), y=indeks, fill=loc))+
  geom_boxplot()+
  scale_x_discrete(name=" ", labels=c("kaalikas", "kartul", "peet", "porgand", "maapirn"))+
  scale_fill_manual(values=c("coral", "cadetblue1"),name="",breaks=c("k", "s"),labels=c("koor", "sisu"))+
  scale_y_continuous("Shannoni indeks")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(, size = 18, face = "bold"))

sh50otu

y10.11<-y10%>%select(-Nr.exps)%>%dcast(exp1~OTU)
y10.11$exp1<-as.character(y10.11$exp1)
y10.11[is.na(y10.11)]<-0
y10.1<-y10.11%>%separate(exp1, into=c("ex", "location"), sep= -2)
y10.matrix<-y10.1[,-(1:2)]
shannon_y10_OTU<-diversity(y10.matrix, index="shannon", MARGIN=1)%>%as.data.frame()
shannon_y10_OTU<-cbind(ex=y10.1$ex, loc=y10.1$location, shannon_y10_OTU)
names(shannon_y10_OTU)[3]<-paste("indeks") #annab veerule uue nime
###koor ja sisu III tase 5 juurvilja koor ja sisu
sh10otu<-ggplot(data=shannon_y10_OTU, aes(x=factor(loc), y=indeks, fill=factor(loc)))+
  geom_boxplot()+
  scale_fill_manual(values=c("coral", "cadetblue1"),name="",breaks=c("k", "s"),labels=c("koor", "sisu"))+
  scale_x_discrete(name="", labels=c("koor", "sisu"))+
  scale_y_continuous("Shannoni indeks")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(, size = 16, face = "bold"))
sh10otu

write.csv2(shannon_y10_OTU, "sh10_OTU.csv")
write.csv2(shannon_y50_OTU, "sh50_OTU.csv")

#korrelatsioonid
library("corrgram")
source("sparCC.R")
#########SparCC
## count matrix x should be samples on the rows and OTUs on the colums,
## assuming dim(x) -> samples by OTUs

## count matrix x should be samples on the rows and OTUs on the colums,
## assuming dim(x) -> samples by OTUs

#sparcc <- function(x, max.iter=20, th=0.1, exiter=10)
## 10 katsega OTU-de oma, aluseks shannoni tabel y10.11
y10.matrix<-y10.11[,-1]
y10.matrix<-y10.matrix[,sapply(y10.matrix, function(x) any(x>5))]
spar<-sparcc(y10.matrix)
names.y10.matrix<-names(y10.matrix)
spar.cor<-spar$CORR
colnames(spar.cor) <- names.y10.matrix
rownames(spar.cor) <- names.y10.matrix
corrgram(spar.cor,order=TRUE, 
         main=" ", 
         lower.panel=panel.shade, 
         upper.panel=panel.pie, 
         diag.panel=NULL, 
         text.panel=panel.txt)



