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

## Clustering avec lga

## Clustering manuel

