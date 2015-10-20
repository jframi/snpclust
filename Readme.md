# Installation

Installer le package devtools si vous ne l'avez pas déjà

```{r}
install.packages(devtools)
```
Charger le package devtools puis installer le package snpclust
```{r}
library(devtools)
install_github("https://github.com/jframi/snpclust")
```
Puis charger le package snpclust

```{r}
library(snpclust)
```


# Usage

```{r}
sonde<-read.table("Sonde1.txt",header=T,sep="\t",dec=",")
```

## Clustering avec cmeans

clust.cmeans(sonde,Rmin = 1.5,seuil = .6)

## Clustering avec lga

clust.lga(sonde,niter=1000,scale=F,Rmin = 1.5)

## Clustering manuel

On execute une à une les commandes suivantes pour définir chaque groupe d'alleles:
Après avoir déssiné un polygone, on clique une valeur de X negative pour quitter

```{r}
sonde<-man.reclust(sonde,what="Allele_X",update.all = T)
sonde<-man.reclust(sonde,what="Allele_Y",update.all = T)
sonde<-man.reclust(sonde,what="Both_Alleles",update.all = T)
sonde<-man.reclust(sonde,what="Negative",update.all = T)
sonde<-man.reclust(sonde,what="Unknown",update.all = F)
```
Si `update.all=T` tous les individus du cluster considéré du sont d'abord passés en Unknown
Si `update.all=F` les individus désignés par la polygone sont ajoutés au cluster

# A Faire

Améliorer clust.cmeans et clust.lga
Rajouter l'update de la colonne Call dans clust.cmeans et clust.lga



