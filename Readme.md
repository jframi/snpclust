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

```{r}
  clust.cmeans(sonde,Rmin = 1.5,seuil = .6)
```

## Clustering avec lga

```{r}
clust.lga(sonde,niter=1000,scale=F,Rmin = 1.5)
```

## Clustering manuel

On execute une à une les commandes suivantes pour définir chaque groupe d'alleles.  
Après avoir dessiné un polygone, on clique une valeur de X negative pour quitter

```{r}
sonde<-man.reclust(sonde,what="Allele_X",update.all = T)
sonde<-man.reclust(sonde,what="Allele_Y",update.all = T)
sonde<-man.reclust(sonde,what="Both_Alleles",update.all = T)
sonde<-man.reclust(sonde,what="Negative",update.all = T)
sonde<-man.reclust(sonde,what="Unknown",update.all = F)
```
Si `update.all=T` tous les individus du cluster considéré du sont d'abord passés en Unknown puis ceux désignés par le polygone sont codés avec l'allele concerné.  
Si `update.all=F` les individus désignés par le polygone sont codés avec l'allele concerné sans que les autres soient modifiés. 

# A Faire

- Améliorer clust.cmeans et clust.lga   
- Rajouter l'update de la colonne Call dans clust.cmeans et clust.lga  
- __Faire le clustering cmeans en R/Theta__  


