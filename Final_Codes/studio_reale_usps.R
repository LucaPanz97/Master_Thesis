### COVGLASSO_TCLUST ###

setwd("C:/Users/panze/Desktop/TESI/Personal_Files")

library("rhdf5")
h5ls("C:/Users/panze/Desktop/TESI/Personal_Files/usps.h5")
train <- h5read("C:/Users/panze/Desktop/TESI/Personal_Files/usps.h5", "train")
usps_trasposto <- train$data
usps <- t(usps_trasposto)
label <- train$target


### 0,1 e 4 ###

i0 <- which(label==0)
i1 <- which(label==1)
i4 <- which(label==4)
iout <- which(label!=0 & label!=1 & label!=4)

usps0 <- usps[i0,]
usps1 <- usps[i1,]
usps4 <- usps[i4,]
uspsout <- usps[iout,]

selected_rows_0 <- sample(nrow(usps0), size = 50, replace = FALSE)
usps0_filtered <- 2*usps0[selected_rows_0, ]

selected_rows_1 <- sample(nrow(usps1), size = 50, replace = FALSE)
usps1_filtered <- 2*usps1[selected_rows_1, ]

selected_rows_4 <- sample(nrow(usps4), size = 50, replace = FALSE)
usps4_filtered <- 2*usps4[selected_rows_4, ]

selected_rows_out <- sample(nrow(uspsout), size = 5, replace = FALSE)
uspsout_filtered <- 2*uspsout[selected_rows_out, ]

#vedere che cifre sono gli outliers che abbiamo selezionato
labelout <- label[which(label!=0 & label!=1 & label!=4)]
labeloutsel <- labelout[selected_rows_out]
labeloutsel

subfin <- rbind(usps0_filtered, usps1_filtered, usps4_filtered, uspsout_filtered)
group <- factor(c(rep(0,50),rep(1,50),rep(4,50),rep("out",5)))

#scelgo di eliminare le covariate con varianza vicina allo zero 
#variabili_basse_varianze <- caret::nearZeroVar(subfin) #34 variabili, troppo poche da eliminare

#scelgo di tenere solo le covariate con varianza sopra una certa soglia
varianze_colonne <- apply(subfin, 2, var)
subfin_filtered <- subfin[, varianze_colonne > 0.55]
p <- dim(subfin_filtered)[2]
plot(subfin_filtered[,26:27], col = group)

tclust <- .covglasso_in_tclust_newinit(subfin_filtered, k = 3, alpha = 5/155, nstart = 20, iter.max = 8, equal.weights = FALSE, zero.tol = 1e-16, lambda = 10, P = matrix(1, nrow = p, ncol = p))

#vedere quanto matcha 
matching_counts <- table(group, tclust$cluster)
matching_counts
# accuracy = (5 + 48 + 50 + 48) / 155 = 151 / 155 = 0.97

#Calcolo dell'Adjusted Rand Index (ARI)
library(mclust)
ari <- adjustedRandIndex(group, tclust$cluster)
print(ari)
# ari = 0.92

#Salvare il modello su disco
saveRDS(tclust, file = "covglassotclust014.rds")
#Caricare il modello precedentemente salvato
tclust <- readRDS(file = "covglassotclust014.rds")


### 3,5 e 8 ###

i3 <- which(label==3)
i5 <- which(label==5)
i8 <- which(label==8)
iout <- which(label!=3 & label!=5 & label!=8)

usps3 <- usps[i3,]
usps5 <- usps[i5,]
usps8 <- usps[i8,]
uspsout <- usps[iout,]

#data("usps358")
#notiamo che usps358 è effettivamente il subset di train considerando solo le cifre 3,5 e 8 (il numero di data points 1756 coincide)
#tuttavia i valori delle variabili in usps358 sono compresi tra 0 a 2 invece che tra 0 e 1 come nel train
#e notiamo che sono semplicemente i valori del subset di train moltiplicati per 2, forse per allargare il dominio delle variabili e avere così data points leggermente più distanziati

selected_rows_3 <- sample(nrow(usps3), size = 50, replace = FALSE)
usps3_filtered <- 2*usps3[selected_rows_3, ]

selected_rows_5 <- sample(nrow(usps5), size = 50, replace = FALSE)
usps5_filtered <- 2*usps5[selected_rows_5, ]

selected_rows_8 <- sample(nrow(usps8), size = 50, replace = FALSE)
usps8_filtered <- 2*usps8[selected_rows_8, ]

selected_rows_out <- sample(nrow(uspsout), size = 5, replace = FALSE)
uspsout_filtered <- 2*uspsout[selected_rows_out, ]

#vedere che cifre sono gli outliers che abbiamo selezionato
labelout <- label[which(label!=3 & label!=5 & label!=8)]
labeloutsel <- labelout[selected_rows_out]
labeloutsel

subfin <- rbind(usps3_filtered, usps5_filtered, usps8_filtered, uspsout_filtered)
group <- factor(c(rep(3,50),rep(5,50),rep(8,50),rep("out",5)))

#scelgo di eliminare le covariate con varianza vicina allo zero 
#variabili_basse_varianze <- caret::nearZeroVar(subfin) #41 variabili, troppo poche da eliminare

#scelgo di tenere solo le covariate con varianza sopra una certa soglia
varianze_colonne <- apply(subfin, 2, var)
subfin_filtered <- subfin[, varianze_colonne > 0.50]
p <- dim(subfin_filtered)[2]
plot(subfin_filtered[,26:27], col = group)

tclust <- .covglasso_in_tclust_newinit(subfin_filtered, k = 3, alpha = 5/155, nstart = 20, iter.max = 10, equal.weights = FALSE, zero.tol = 1e-16, lambda = 10, P = matrix(1, nrow = p, ncol = p))

#vedere quanto matcha 
matching_counts <- table(group, tclust$cluster)
matching_counts
# accuracy = (5 + 43 + 23 + 21) / 155 = 92 / 155 = 0.60

#Calcolo dell'Adjusted Rand Index (ARI)
library(mclust)
ari <- adjustedRandIndex(group, tclust$cluster)
print(ari)
# ari = 0.25

#Salvare il modello su disco
saveRDS(tclust, file = "covglassotclust358.rds")
#Caricare il modello precedentemente salvato
tclust <- readRDS(file = "covglassotclust358.rds")


######################################################################################################################################################

### LEDOITWOLF_TCLUST ###

setwd("C:/Users/panze/Desktop/TESI/Personal_Files")

library("rhdf5")
h5ls("C:/Users/panze/Desktop/TESI/Personal_Files/usps.h5")
train <- h5read("C:/Users/panze/Desktop/TESI/Personal_Files/usps.h5", "train")
usps_trasposto <- train$data
usps <- t(usps_trasposto)
label <- train$target


### 0,1 e 4 ###

i0 <- which(label==0)
i1 <- which(label==1)
i4 <- which(label==4)
iout <- which(label!=0 & label!=1 & label!=4)

usps0 <- usps[i0,]
usps1 <- usps[i1,]
usps4 <- usps[i4,]
uspsout <- usps[iout,]

selected_rows_0 <- sample(nrow(usps0), size = 50, replace = FALSE)
usps0_filtered <- 5*usps0[selected_rows_0, ]

selected_rows_1 <- sample(nrow(usps1), size = 50, replace = FALSE)
usps1_filtered <- 5*usps1[selected_rows_1, ]

selected_rows_4 <- sample(nrow(usps4), size = 50, replace = FALSE)
usps4_filtered <- 5*usps4[selected_rows_4, ]

selected_rows_out <- sample(nrow(uspsout), size = 5, replace = FALSE)
uspsout_filtered <- 5*uspsout[selected_rows_out, ]

#5* è la trasformazione con cui otteniamo i risultati migliori nel problema di clustering con cifre 0,1,4 utilizzando il LedoitWolf_tclust 

#vedere che cifre sono gli outliers che abbiamo selezionato
labelout <- label[which(label!=0 & label!=1 & label!=4)]
labeloutsel <- labelout[selected_rows_out]
labeloutsel

subfin <- rbind(usps0_filtered, usps1_filtered, usps4_filtered, uspsout_filtered)
group <- factor(c(rep(0,50),rep(1,50),rep(4,50),rep("out",5)))

#scelgo di eliminare le covariate con varianza vicina allo zero 
#variabili_basse_varianze <- caret::nearZeroVar(subfin) #34 variabili, troppo poche da eliminare

#scelgo di tenere solo le covariate con varianza sopra una certa soglia
varianze_colonne <- apply(subfin, 2, var)
subfin_filtered <- subfin[, varianze_colonne > 3.5]
p <- dim(subfin_filtered)[2]
plot(subfin_filtered[,26:27], col = group)

tclust <- .LedoitWolf_in_tclust_newinit(subfin_filtered, k = 3, alpha = 5/155, nstart = 150, iter.max = 20, equal.weights = FALSE, zero.tol = 1e-16)

#vedere quanto matcha 
matching_counts <- table(group, tclust$cluster)
matching_counts
# accuracy = (5 + 45 + 50 + 40) / 155 = 140 / 155 = 0.90

#Calcolo dell'Adjusted Rand Index (ARI)
library(mclust)
ari <- adjustedRandIndex(group, tclust$cluster)
print(ari)
# ari = 0.73

#Salvare il modello su disco
saveRDS(tclust, file = "LedoitWolftclust014.rds")
#Caricare il modello precedentemente salvato
tclust <- readRDS(file = "LedoitWolftclust014.rds")


### 3,5 e 8 ###

i3 <- which(label==3)
i5 <- which(label==5)
i8 <- which(label==8)
iout <- which(label!=3 & label!=5 & label!=8)

usps3 <- usps[i3,]
usps5 <- usps[i5,]
usps8 <- usps[i8,]
uspsout <- usps[iout,]

selected_rows_3 <- sample(nrow(usps3), size = 50, replace = FALSE)
usps3_filtered <- 5*usps3[selected_rows_3, ]

selected_rows_5 <- sample(nrow(usps5), size = 50, replace = FALSE)
usps5_filtered <- 5*usps5[selected_rows_5, ]

selected_rows_8 <- sample(nrow(usps8), size = 50, replace = FALSE)
usps8_filtered <- 5*usps8[selected_rows_8, ]

selected_rows_out <- sample(nrow(uspsout), size = 5, replace = FALSE)
uspsout_filtered <- 5*uspsout[selected_rows_out, ]

#vedere che cifre sono gli outliers che abbiamo selezionato
labelout <- label[which(label!=3 & label!=5 & label!=8)]
labeloutsel <- labelout[selected_rows_out]
labeloutsel

subfin <- rbind(usps3_filtered, usps5_filtered, usps8_filtered, uspsout_filtered)
group <- factor(c(rep(3,50),rep(5,50),rep(8,50),rep("out",5)))

#scelgo di eliminare le covariate con varianza vicina allo zero 
#variabili_basse_varianze <- caret::nearZeroVar(subfin) #41 variabili, troppo poche da eliminare

#scelgo di tenere solo le covariate con varianza sopra una certa soglia
varianze_colonne <- apply(subfin, 2, var)
subfin_filtered <- subfin[, varianze_colonne > 3.2]
p <- dim(subfin_filtered)[2]
plot(subfin_filtered[,26:27], col = group)

tclust <- .LedoitWolf_in_tclust_newinit(subfin_filtered, k = 3, alpha = 5/155, nstart = 150, iter.max = 20, equal.weights = FALSE, zero.tol = 1e-16)

#vedere quanto matcha 
matching_counts <- table(group, tclust$cluster)
matching_counts
# accuracy = (5 + 29 + 29 + 30) / 155 = 93 / 155 = 0.60

#Calcolo dell'Adjusted Rand Index (ARI)
library(mclust)
ari <- adjustedRandIndex(group, tclust$cluster)
print(ari)
# ari = 0.17

#Salvare il modello su disco
saveRDS(tclust, file = "LedoitWolftclust358.rds")
#Caricare il modello precedentemente salvato
tclust <- readRDS(file = "LedoitWolftclust358.rds")


### Altre misure per valutare la bontà del clustering ###

#indice di validità di Silhouette: misura quanto bene ogni oggetto si adatta al proprio cluster rispetto agli altri cluster
library(cluster)
silhouette_values <- silhouette(tclust$cluster)
mean_silhouette <- mean(silhouette_values[, "silhouette.width"])

#coefficiente di correlazione di Pearson: misura la correlazione tra la distanza euclidea dei punti nello spazio originale e la distanza tra i punti nello spazio di proiezione
library(fpc)
correlation <- cluster.stats(subfin_filtered, tclust$cluster)$pearson

#indice di Dunn: misura la distanza minima tra i centroidi di ogni coppia di cluster rispetto alla massima distanza tra punti all'interno di ogni cluster
library(fpc)
dunn_index <- dunn(subfin_filtered, tclust$cluster)

#indice di validità di Xie-Beni: misura il rapporto tra la somma delle distanze quadrate tra i punti e i centroidi dei cluster e la distanza tra i punti e i centroidi dei cluster più vicini 
library(fpc)
xie_beni_index <- index.G1(subfin_filtered, tclust$cluster)





