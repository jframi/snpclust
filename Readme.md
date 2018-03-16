# Version 0.2 (01/09/2017)

## Installation

Install devtools package if not already done

```{r}
install.packages(devtools)
```
Load devtools then install snpclust package

```{r}
library(devtools)
install_github("jframi/snpclust")
```
Then load snpclust

```{r}
library(snpclust)
```


## Usage


### Shiny app for manual clustering

After installing the package, the shiny app can be run using:

```{r}
runmanclust()
```


### Other functions that may be useful for manual clustering

```{r}
sonde<-read.table("Sonde1.txt",header=T,sep="\t",dec=",")
```

#### Clustering with cmeans

```{r}
clust.cmeans(sonde,Rmin = 1.5,seuil = .6)
```

#### Clustering with lga

```{r}
clust.lga(sonde,niter=1000,scale=F,Rmin = 1.5)
```

#### Manual clustering (should use shiny app since v0.2)

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

## To do

- Improve clust.cmeans et clust.lga   
- Add update of Call column in clust.cmeans and clust.lga  
- __Use R/Theta for cmeans clustering__  


