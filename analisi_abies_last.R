################################################################################
# Questo codice è un esercizio su analisi discriminante e per problemi di verifica di ipotesi di classificazione


################################################################################################
# Fare un controllo sull'installazione dei seguenti pacchetti

# install.packages(c("papeR","rtables","postHoc","PMCMRplus","readr","psych","janitor")

# è possibile installazione manuale del pacchetto postHoc recente da un file postHoc1.4.tar.gz presente nella directory


################################################################################
# Impostare la directory di lavoro  e i semi per le funzioni randomiche 
# per assicurare la riproducibilità


setwd("C:\\aaa_lavori\\lav_waed_abies")

# setwd("")

set.seed(123)

##########################################################################################################
# caricamento dei pacchetti utilizzati nel codice

library(papeR)
library(rtables) # https://rookie.rbind.io/post/making-summary-tables-in-r
library(data.table) 
library(postHoc) # https://cran.r-project.org/web/packages/postHoc/vignettes/Post-hoc-analysis.html
library(PMCMRplus)

source("load_libraries_abies.R")

#######################################################################################################
# pulisco la console da riga di comando

cat("\014") 

################################################################################
# Carico 

dati_sel=read.xlsx("dati_abies_GC.xlsx",1) # dati relativi alle analisi GC MS

dati_seeds=read.xlsx("dati_abies_semi.xlsx",1) # dati relativi ai paraemtri morfologici dei semi


cat("\014") 

################################################################################
# eliminazione dei dati outlier e salvataggio su un file di archiviazione R formato rds 

# dati_sel=dati_sel[-c(165,167,162,20,45),] 

saveRDS(dati_sel,"dati_sel_abies.rds")

################################################################################
# Visualizzazione della matrice di lavoro GC MS e dei  nomi delle variabili


View(dati_sel) # vedere i dati 

names(dati_sel)  # vedere i nomi delle matrici 

# [1] "ID"                    "Species"               "date"                  "a-pinene"              "camphene"             
# [6] "b-pinene"              "myrcene"               "Limonene"              "β-Phellandrene"        "g-terpinene"          
# [11] "p-cymene"              "terpinolene"           "Limonene_oxide_cis"    "Pinocarvone"           "Unk_1"                
# [16] "bornylacetate"         "4-ol-terpinen"         "trans-Verbenol"        "verbenone"             "b-elemene"            
# [21] "sesquiterpeni.1"       "sesquiterpene.2"       "a-humulene"            "GermacreneD"           "sesquiterpene.3"      
# [26] "B.caryophillene-oxide" "1.Germacrene.D-4-ol"   "unk_.2"                "α-epi-Cadinol"         "Selina-6-en-4-ol"     
# [31] "α-Cubebene"            "a-copaene"             "b-copaene"             "TOT_Mono"              "TOT_mono&sesqui"    


nrow(dati_sel) #  170

#######################################################################################################
# Normalizzazione dei dati sul totale dei mono&sesqui perchè parto dal dato delle aree

dati_sel_rel=100*dati_sel[,4:31]/dati_sel$`TOT_mono&sesqui`

mat_final_rel=data.frame(species=dati_sel$Species,dati_sel_rel)

write.xlsx(mat_final_rel,"dati_relativi_totali.xlsx",overwrite = T) 
 
# solo dei mono e dei sesqui separati
dati_monosesqui=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,
                      100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`-dati_sel$TOT_Mono))

# solo dei mono e dei sesqui sommati per fare pesare meno i sesqui

dati_monosesquimiche=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,
                           100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`+dati_sel$TOT_Mono))

# la scelta della normalizzazione dipende dalla situazione e dagli obbiettivi degli analisi

###############################################################################################################################
# controllo di collinearità opzionale in questo caso non utilizzato

dati_sel_rel=dati_monosesquimiche

CN(dati_sel_rel) # 
multiCol(dati_sel_rel, graf = TRUE)
names(dati_sel_rel)[c(20,16,23)]

dati_sel_rel=dati_sel_rel[,-c(20,16,23)]
CN(dati_sel_rel) 
multiCol(dati_sel_rel, graf = TRUE)

##########################################################################

dati_sel_rel=dati_monosesquimiche # ignoro la multicollinearità ricaricando i dati normalizzati

##########################################################################
# Creazione delle matrici di analisi per PCA

X=dati_sel_rel
Y=dati_sel$Species

# X=asin(sqrt(X/100)) # eventuale trasfomazione dei dati normalizzati 


Xseed=dati_seeds[,4:8]
Yseed=dati_seeds$species

write.xlsx(list(data.frame(Y,X),data.frame(Yseed,Xseed)),"dati_PCA.xlsx")

cat("\014") 

#########################################################################################
# PCA explore variable and outlier detection

data_pca=data.frame(X,Species=Y)


abies_pca=ordinate(data_pca , cols = 1:25, model = ~ prcomp(., scale. = TRUE)) 

abies_pca %>%
  augment_ord() %>%
  ggbiplot(axis.type = "predictive",aes(alpha = .3)) +
  theme_bw() +
  geom_rows_point(aes(color=Species,shape=Species))+
  stat_rows_center(alpha = .8, fun.center = "mean",size = 5,aes(color = Species))+
  geom_origin() +
  ggtitle("Principal Component Analisys Abies data") +
  labs(color = "Species")

ggsave("PCA_biplot.png")


#########################################################################################
# data selection by non parametric ANOVA o Kruskal Wallis ( optional) 

a=col_kruskalwallis(X,Y) # uso il pacchetto MatrixTests

row.names(a)[which(a$pvalue<0.05)] # all coumpound discriminate!
row.names(a)[which(a$pvalue>0.05)] 

write.xlsx(data.frame(a,compound_terpene=row.names(a)),"test_kruskal.wallis.xlsx")

##########################################################################################################
# Faccio la Posthoc analysis

res_posthoc=list()
res_posthoc_dunn=list()
res_GG=list()

#######################################################

j=1  

for ( i in 2:29) { # run on the columns starting fron numeric values
  
formula_t=paste(names(mat_final_rel)[i], "~ species + 0")

MM=glm(formula_t,  data = mat_final_rel)

GG <- posthoc(MM)


dunn_df=as.data.frame(gsub("\\s+", " ", paste0(names(mat_final_rel[i])," ",capture.output(summary(kwAllPairsDunnTest(mat_final_rel[,i], factor(mat_final_rel$species),p.adjust.method = "bonferroni")))[2:4])))

names(dunn_df)=c("Compound Pair Z p_values sign")
res_posthoc_dunn[[j]]=dunn_df
res_posthoc[[j]]=data.frame(terpene_compound=names(mat_final_rel[i]),
                            GG$CI,groups=GG$Grouping,
                            test_sp=gsub("species","",row.names(GG$PvaluesMatrix)),
                            GG$PvaluesMatrix,
                            ks_pval=KruskalWallisAllPvalues(mat_final_rel[,i], factor(mat_final_rel$species))
                            )
res_GG[[j]]=GG

j=j+1

}

res_df=do.call("rbind",res_posthoc)
res_df_dunn=do.call("rbind",res_posthoc_dunn)

write.xlsx(list(dunn_test=res_df_dunn,t_test=res_df),"posthoc_kruskal.wallis.xlsx")

##############################################################################################
# summary tables

my_desc=psych::describeBy(mat_final_rel[,-1], mat_final_rel$species)

readr::write_excel_csv2(rbindlist(my_desc, idcol = 'species'), file = 'summary_terpene_by_species.xlsx')


######################################################################################
# boxplot main discriminant compounds

dir.create("boxplot")

setwd("boxplot")

dati_sel_box=dati_sel_rel
dati_sel_box$Species=factor(dati_sel$Species)
dati_sel_box=janitor::clean_names(dati_sel_box) # pulisco i nomi

for ( i in 1:25) {

ggboxplot(dati_sel_box,"species",names(dati_sel_box)[i],fill="red") +ylim(0,50)
  
ggsave(paste0("boxplot_",names(dati_sel_box)[i],".png"))

}

setwd("..")



# my_plot_list=list(a,b,c)
# ggexport(
#   plotlist =my_plot_list, filename = "composti_selected.pdf", 
#   ncol = 1, nrow = 4
# )


####################################################################################
# Ordination Supervised methods : LDA & & KNN

# LDA
# Prior probabilities of groups: the proportion of training observations in each group. For example, there are 31% of the training observations in the setosa group
# Group means: group center of gravity. Shows the mean of each variable in each group.
# Coefficients of linear discriminants: Shows the linear combination of predictor variables that are used to form the LDA d

#########################################################################

training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)

dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model <- lda(Y~., data =dati_train)

predictions_lda <- model %>% predict(dati_test)
model_accuracy_lda=mean(predictions_lda$class==dati_test$Y)
res_lda_abiesGC=confusionMatrix(predictions_lda$class,factor(dati_test$Y))

###########################################################################################################################################
# plots LDA

abies_lda <- lda_ord(X, Y, axes.scale = "standardized")
abies_lda %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  mutate_rows(
    species = grouping,
    data = ifelse(.element == "active", "centroid", "case")
  ) %>%
  ggbiplot() +
  theme_bw() +
  geom_rows_point(aes(
    color = grouping,
    size = data, 
    alpha = data
  )) +
  ggtitle("Standardized coefficient biplot of mediterranean Abies spp. ") +
  expand_limits(y = c(-3, 5))+
  labs(color = "Species")
ggsave("LDA_biplot.png")

################################################################################
# knn classifiers

test_pred <- class::knn(
  train = dati_train[,1:25], 
  test = dati_test[,1:25],
  cl = dati_train$Y, 
  k=10
)
res_knn_abiesGC=confusionMatrix(test_pred,factor(dati_test$Y))

################################################################################
# Ora sto lavorando sui dati dei semi


col_oneway_welch(Xseed, Yseed) # utilizzo il pacchetto MatrixTests per fare l'ANOVA a 1 via 

summarySE(Xseed,"length","Yseed")
summarySE(Xseed,"mc","Yseed") 
summarySE(Xseed,"weight100","Yseed") 
summarySE(Xseed,"length","Yseed") 
summarySE(Xseed,"width","Yseed")

X=Xseed
Y=Yseed

################################################################################
# PCA

ord <- PCA(Xseed, graph = FALSE)
ggord(ord, Yseed,arrow=NULL,txt=NULL)

################################################################################
# linear DA 

heplots::boxM(X, Y)

training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)
dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model <- lda(Y~., data =dati_train)
predictions_lda <- model %>% predict(dati_test)
model_accuracy_lda=confusionMatrix(predictions_lda$class,factor(dati_test$Y))

###########################################################################################################################################
# plots DA seeds

ggord(model, Y[training.samples],
      arrow=NULL,
      txt=NULL,
      xlims=c(-8,8),
      ylims =c(-5,7),
      obslab =F)

#####################################################################
# References

# https://cmdlinetips.com/2020/12/canonical-correlation-analysis-in-r/
# https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/
# https://vitalflux.com/pca-vs-lda-differences-plots-examples/
# https://towardsai.net/p/data-science/lda-vs-pca
# https://stats.stackexchange.com/questions/23353/pca-lda-cca-and-pls
# https://www.geeksforgeeks.org/classifying-data-using-support-vector-machinessvms-in-r/
# https://mdatools.com/docs/pca.html
# https://rdrr.io/cran/mixOmics/man/plsda.html
# http://mixomics.org/methods/spls-da/
