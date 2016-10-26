# omics_scripts
My analysis scripts for Genomics project

它可以输入中文吗？似乎是可以的^_^
有时中文注释可以用来加密，哈哈

the scripts are written for the genomics project, Including the following steps, or more.

##1. Detect somatic mutation and indels
##2. Find driver gene
##3. mutation signatures and their associated clinical and genetic variables
##4. clinical subtype classification
##5. anymore?


##6. plotting for genomics paper
R scripts were written to plot similar state of art figures as similar as the ones published in high-level journal, for example, Nature, Nature Genetics, etc. Many of the scripts were under the `ggplot2` plotting language framework, and combined and adjusted by `grid`, `gridExtra` and `gtable` packages. 

### Before start
As my strong recommendation, one could be falimiar with the `dplyr` and `reshape2` packages because they are very useful, readable and elegant for data manipulation and debuging. All packages can be installed locally by use `install.packages()` function. e.g.

```R
install.packages("ggplot2")
```

Once finished installing these packages, enjoy my journal of plotting !

#### 1. oncoplot
The oncoplot is a rainfall like combination of several figures, including heatmap and barplot. The main part is a heatmap to show the occurence of mutations and their types for each gene on each sample. the top and the right(left) are the marginal barplot for sample and gene, respectively.

One Nature published example is like this:
![image](https://github.com/Hiuyu/omics_scripts/blob/oncoplot/oncoplot/example_hope_to_plot.jpg?raw=true)
