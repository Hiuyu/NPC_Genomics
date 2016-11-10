## demo for heatmap with barplot on top and right panels, 
## also called the rainfall plot or oncoplot
## the basic steps are as follows:
##    1. create each plots seperately
##    2. save the legend of one plot as a seperate grob
##    3. strip the legends from the plots
##    4. arrange plots by gtable or arrangeGrob 
##    5. draw the four plots and the legend below or above (or any other) them

## 2016-10-26
## problems remained to resolved:
##    1. increase spacing between legend keys
##    2. orientation of gene names
##    3. let y-axis of the right plot show on top
##    4. how to add a left vertical plot ?
##    
## 2016-10-28
## solution to problem 3 and 4:
##    use scale_y_reverse() to make the Y axis reverse. 
## 2016-10-30
## To apply it into practice
##    1. assume the input is a tabular format, in which columns with name sample, gene and mutation_type, are required.
##    2. write a convert function to recode the mutation type given a dictionary
##    3. write a function to generate analysis-ready data (onco): make_oncoMatrix()
##    4. if there are more than one mutation in a gene in the same sample, assign this gene a new category name "multi-hit".
## 2016-10-31
## Fixed some bugs and add new arguments
##    1. allow to choose whether to include silent or unknown mutation
##    2. fixed a bug that make the sample order not correct.
##    3. allow to input presorted sample order and gene order. could be a subset of genes or samples which are of interest to plot only.
##    4. allow removed samples with no mutation in genes included, when plotting the heatmap.
## 2016-11-01
## Fixed bugs and add new features
##    1. if included subsets of genes, the top bar will still count mutations on all genes genes.
##    2. if included subsets of samples, no genes were dropped and the top 30 were included by default
##    3. add the left plot for mutsig log(P) or pathway assignment, by argument: gene.annotation.plot. required gene column
##    4. allow to add optional sample annotation informations by argument: sample.annotation.plot. required sample column
## 2016-11-05
## add new argument
##    1. add sample annotations below the heatmap, with the legend on the bottomright cell. 
##    2. The argument sample.annotation.table is a data.frame with columns: sample, annotation1, annotation2, .....
##

oncoplot <- function(data, # main data, a data.frame
                     gene.annotation.plot, # topleft plot ggplot2 object
                     sample.annotation.table, # bottom plot for sample annotation, a data.frame
                     included.gene.list, included.sample.list, 
                     is.sort.gene = TRUE,
                     is.sort.sample = TRUE,
                     include.silent = TRUE,
                     gene.frac.threshold = 0.03,
                     gene.n.threshold = 5,
                     include.gene.n = 30, #  top 30 to included
                     is.drop.gene = TRUE,
                     include.unknown = FALSE){ # just for test

require(ggplot2)
require(grid)
require(gridExtra)
require(gtable)
require(reshape2)

# a demo mutation data
# set dictionary of mutation type
dictType=c("missense", 
           "nonsense",
           "silent",
           "frameshift",
           "non-frameshift",
           "multi-hit"
           )
names(dictType)=1:6
# define the dictionary for annovar type mutation category
dict.annovar.ExonicFunc = c(
  "nonsynonymous_SNV" = "missense",
  "synonymous_SNV"  = "silent",       
  "stopgain" = "nonsense",
  "stoploss" = "nonsense",
  "unknown" = "unknown",
  "frameshift_deletion" = "frameshift",
  "frameshift_insertion" = "frameshift",
  "nonframeshift_deletion" = "non-frameshift",
  "nonframeshift_insertion" = "non-frameshift",
  "non-exonic" = NA
)

# set custom theme()
mythm=theme(
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x=element_blank(),
  axis.ticks=element_blank(),
  panel.grid=element_blank(),
  panel.background=element_rect(fill="white",colour = NA),
  panel.border=element_rect(fill=NA,colour="grey50")
)

# how to increase the spacing of legend key ?
legend_theme <- theme(
  legend.position="bottom",
  legend.direction="horizontal",
  legend.text=element_text(size=11),
  legend.title=element_blank(),
  legend.key.height=unit(0.4,"cm"),
  legend.key.width=unit(0.3,"cm")
)

# set manual color
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## function to convert mutation for a given dictonary
data$mutation = dict.annovar.ExonicFunc[data$mutation]

## generate oncoMatrix
make_oncoMatrix <- function(data = data){
  require(reshape2)
  # the fun.aggregate function for dcast, used to aggregating multiple matched entries in one cell
  .func.aggregate <- function(hits) {
    # for simplicity, all were recoded as "multi-hit"
    if(length(hits) == 1){
      return(hits)
    } else if (length(hits) > 1) {
      return("multi-hit")
    } else {
      return("")
    }
  }
  # convert to matrix
  data1 = dcast(data, sample ~ gene, value.var = "mutation", fun.aggregate = .func.aggregate)
  # convert back to tabular format
  onco=melt(data1, id.vars = "sample", variable.name = "gene", value.name = "mutation")
  onco$mutation[onco$mutation == ""] = NA # convert non-mutation cell as NA, required for heatmap
  onco$mutation=factor(onco$mutation,levels = c(dictType)) # assign label
  onco$gene = as.character(onco$gene)
  return(onco)
}

## remove unknown and silent or not ?
if ( ! include.unknown) {
  data = data %>% filter(mutation != "unknown")
}
if ( ! include.silent) {
  data = data %>% filter(mutation != "synonymous_SNV")
}

## convert to oncoMatrix
onco = make_oncoMatrix(data)


## how to decide gene order ?
# or just provide a pre-sort gene list?
oncoplot_gene_order = function(n_threshold ,f_threshold) {
  order_gene = onco %>% 
    filter(!is.na(mutation)) %>% 
    group_by(gene) %>% 
    summarise(n=n()) %>% 
    arrange(-n)
  order_gene=as.data.frame(order_gene)
  order_gene$f = order_gene$n / length(unique(onco$sample))
  # if not set, 1% was chosen
  if(missing(f_threshold)){
    f_threshold = 0.01
  }
  # if not set, 5 was chosen
  if(missing(n_threshold)){
    n_threshold = 5
  }
  # if samples size <= 5, do not drop genes
  if(length(unique(onco$sample)) < 6) {
    f_threshold = 0
    n_threshold = 0
  }
  # select only subset of genes
  if(is.drop.gene) { # select only top genes
    order_gene = order_gene %>%
      filter(n >= n_threshold & f >= f_threshold) %>% # any of the threshold passed
      select(gene)
  }
  gene_included = as.character(order_gene$gene)
  return(gene_included)
}

## how to sort sample order ?
# or just provide a pre-sort sample list?
# or sort by frequency ?
oncoplot_sample_order = function(gene_order, weight) {
  if(missing(weight)){
    weight=rep(1,nlevels(onco$mutation)) # a toy set of weight
    names(weight)=levels(onco$mutation)
  }
  gene.value=length(gene_order):1
  names(gene.value)=gene_order
  sample_list = unique(onco$sample) # get sample list
  # compute score for each sample
  sample_score = onco %>%
      filter(!is.na(mutation)) %>%
      filter(gene %in% gene_order) %>%
      mutate(m.score=weight[mutation], 
             g.score=gene.value[gene], 
             score=m.score * g.score
      ) %>% group_by(sample) %>% 
    summarise(sum.score = max(score)) %>%
    arrange(-sum.score)
  sample_score = as.data.frame(sample_score)
  return(sample_score$sample)
}

## perform subsetting if defined included gene or sample sets
if (!missing(included.sample.list)) {
  onco = onco %>% filter(sample %in% included.sample.list)
}else{
  included.sample.list = unique(onco$sample)
}

## counting sample mutations for top barplot
tmp.z = with(onco,table(sample,mutation))
sample.mutation.count = as.data.frame(tmp.z)  # tabulate sample vs mutation 


if (!missing(included.gene.list)) {
  onco = onco %>% filter(gene %in% included.gene.list)
  is.drop.gene = FALSE # if defined included subset, we should not drop genes anymore
}else{
  included.gene.list = unique(onco$gene)
}


## perform sorting if defined is.sort.gene or is.sort.gene
gene_order = NULL
sample_order = NULL
# define the order of gene
if (is.sort.gene) {
  if (is.drop.gene) {
    gene_order = oncoplot_gene_order(f_threshold = gene.frac.threshold, n_threshold = gene.n.threshold)
  } else {
    gene_order = oncoplot_gene_order(f_threshold = 0, n_threshold = 0)
  }
} else {
  
  gene_order = included.gene.list
}
# define the order of sample
if (is.sort.sample) {
  sample_order = oncoplot_sample_order(gene_order)
} else {
  sample_order = included.sample.list
}

# reorder and subset onco.
# arrange sample order by factor levels
onco$sample = factor(onco$sample, levels = sample_order)
sample.mutation.count$sample = factor(sample.mutation.count$sample, levels=sample_order)

# arrange gene order
onco$gene = factor(onco$gene, levels = rev(gene_order)) # to make high order on top of heatmap

## subsetting onco if only subset of genes or samples are included
onco = onco %>% filter(!is.na(sample)) %>% filter(!is.na(gene))
sample.mutation.count = sample.mutation.count %>% filter(!is.na(sample))

# draw heatmap, the body part
Hp <- ggplot(onco, aes(x=sample,y=gene,fill=mutation)) + 
  geom_tile(width=1,height=1,colour="grey70") +
  scale_fill_manual(values = cbPalette, na.value="white") +
  mythm + 
  theme(panel.border=element_rect(size=1, color="black"),
        axis.text.y=element_text(face="italic", size=11),
        legend.position="none"
  )

# draw sample status, the top part
Tp <- ggplot(sample.mutation.count, aes(x=sample, y=Freq,fill=mutation)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=cbPalette) + # set fill color
  scale_y_continuous(breaks=seq(0,5000,500),expand = c(0, 0)) +  # set y axis
  mythm + theme(legend.position="none", 
        rect=element_blank(), 
        axis.line.y=element_line(size=0.8), 
        axis.line.x=element_line(size=0.5),
        axis.ticks.y=element_line(size=0.3,color="black"),
        # how to add tick on Y-axis ?
       # axis.ticks.length=unit(1,"cm"),
        axis.text.y=element_text(face="plain", size=10)
        )

# draw gene status, the right part 
# how to add a top y-axis line and tick marker? use secondary axis?
Rp <- ggplot(subset(onco,!is.na(mutation)), aes(x=gene, fill=mutation)) + geom_bar() +
  coord_flip() +  
  scale_fill_manual(values=cbPalette) + # set fill color
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0, 0)) +  # set y axis
  mythm + theme(legend.position="none", 
                rect=element_blank(),
                axis.line.y=element_line(size=0.5),
                axis.text.y=element_blank(),
                axis.text.x=element_text(size=10),
                axis.line.x=element_line(size=0.8)
                )


# draw an empty figure, for debuging
empty <- ggplot() + 
  geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

# combine figures
# grob all elements
#leg <- extract_ggplot2_legend(Hp)
gH=ggplotGrob(Hp + legend_theme) # grob the heatmap part, main part of the oncoplot
leg=gH$grobs[which(sapply(gH$grobs,function(x) x$name) == "guide-box")] # extract the legend grob
gH=ggplotGrob(Hp) # grob it again
gT=ggplotGrob(Tp) 
gR=ggplotGrob(Rp)

# let the top and heatmap with same width, right and heatmap with same height
maxWidth=grid::unit.pmax(gT$width[2:5],gH$widths[2:5])
maxHeight=grid::unit.pmax(gR$heights[2:5],gH$heights[2:5])
# let the heights and widths consistent
gT$widths[2:5] <- gH$widths[2:5] <- as.list(maxWidth)
gR$heights[2:5] <- gH$heights[2:5] <- as.list(maxHeight)

# arrange the grobs
# how to arrange the heights sizes ? 
gC <- arrangeGrob(gT,empty,gH,gR,ncol=2,nrow=2, widths=c(7,1), heights=c(1,7))
gC <- gtable_add_rows(gC, heights = unit(3,"line"), pos=0)
gC <- gtable_add_grob(gC, leg, t=1, b=1, l=1, r=1)


#### if there are additional plots, add it to the left and the bottom
if(!missing(gene.annotation.plot)) {
  ## arrange the genes of the left part plot 
  tmp.df = gene.annotation.plot$data 
  tmp.df = gene.annotation.plot$data %>% 
    filter(gene %in% gene_order)
  tmp.gene = unique(tmp.df$gene)
  stopifnot(all(tmp.gene %in% gene_order)) # if not all genes in both data
  tmp.df$gene = factor(tmp.df$gene, levels = rev(gene_order))
  gene.annotation.plot$data = tmp.df # updated new plot data
  # add its grob to the combined plot
  gL = ggplotGrob(gene.annotation.plot + ylab("-log(q-value)") + mythm + 
                    theme(legend.position="none", 
                          rect=element_blank(),
                          axis.line.x=element_line(size=0.5),
                          axis.text.y=element_blank(),
                          axis.text.x=element_text(size=10),
                          axis.line.y=element_blank(),
                          axis.title.x=element_text(size=10)
                    )
  )
  gL$heights[2:5] <- gH$heights[2:5]
  gC <- gtable_add_cols(gC, widths = unit(1,"null"), pos=0)
  gC <- gtable_add_grob(gC, gL, t=3, b=3, l=1, r=1)  
}


## arrange the sample of the bottom part plot
if(!missing(sample.annotation.table)) {
  if(class(sample.annotation.table) == "matrix"){
    sample.annotation.table = as.data.frame(sample.annotation.table)
  }
  # check whether the sample.annotation.table has the sample column
  if(!any(colnames(sample.annotation.table) %in% "sample")){
    stop("the sample.annotation.table argument requires 'sample' column")
  }
  # check whether the samples in sample.annotation.table contains the sample_order
  if(!all(sample_order %in% sample.annotation.table$sample)){
    stop("the sample.annotation.table should cover all samples of the mutation data")
  }
  tmp.sa = sample.annotation.table %>% filter(sample %in% sample_order)
  tmp.sa$sample = factor(tmp.sa$sample, levels=sample_order)
  # iteratively add sample annotation
  annot.title = colnames(tmp.sa)[colnames(tmp.sa) != "sample"]
  for(i in 1:length(annot.title)) {
    tmp.sa.df = melt(tmp.sa[,c("sample",annot.title[i])], id.vars = "sample", variable.name = "annotation", value.name = "value")
    # heatmap type plot
    Bp <- ggplot(tmp.sa.df, aes(x=sample, y = annotation, fill=value)) +
      geom_tile(width=1,height=1,colour="grey70") +
      scale_fill_manual(values = cbPalette, na.value="white") +
      mythm + 
      guides(fill=guide_legend(
        title=annot.title[i], 
        keyheight = unit(0.4,"cm"),
        keywidth = unit(0.3,"cm"),
        title.position="left",
        direction="horizontal",
        ncol=3,byrow=T
      )) +
      theme(panel.border=element_rect(size=1, color="black"),
            axis.text.y=element_text(face="bold", size=11)
      )
    gB=ggplotGrob(Bp) # grob the heatmap part, main part of the oncoplot
    legB=gB$grobs[which(sapply(gB$grobs,function(x) x$name) == "guide-box")] # extract the legend grob
    gB=ggplotGrob(Bp + guides(fill=FALSE)) # grob it again
    gB$widths[2:5] <- gH$widths[2:5]
    gC <- gtable_add_rows(gC, heights = unit(1,"null"), pos=-1)
    gC <- gtable_add_grob(gC, gB, t=3+i, b=3+i, l=2, r=2) # add plot
    gC <- gtable_add_grob(gC,legB, t=3+i, b=3+i,l=3, r=3)
  }
}

#jpeg("oncoplot.jpeg",width=1200,height=2000)
grid.newpage()
grid.draw(gC)
#dev.off()

return(list(sample_order = sample_order,
            gene_order = gene_order,
            gtable.object = gC
            )
       )
}



