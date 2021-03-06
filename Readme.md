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
Choose an input file that contain fluorescence values and calls for one or several SNPs and one or several DNA plates.  
![screenshot01](figure/ss01.png)  
Adjust the file format to read correctly the input file.  
![screenshot01](figure/ss02.png)  
Match the columns of your file to the column that the app is expecting  
<img src="figure/ss03.png" width="300">  
Then select a SNP and possibly a Plate to display a scatter plot of fluorescence values.  
![screenshot01](figure/ss04.png)  
You can switch to Show New Call view, select individual genotypes in the plot using the lasso tool, and recall selected points.  
![screenshot01](figure/ss05.png)  


