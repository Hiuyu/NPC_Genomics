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

oncoplot <- function(data){ # just for test

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
names(dictType)=1:5
# define the dictionary for annovar type mutation category
dict.annovar.ExonicFunc = c(
  "nonsynonymous_SNV" = "missense",
  "synonymous_SNV"  = "silent",       
  "stopgain" = "nonsense",
  "stoploss" = "nonsense",
  "unknown" = NA,
  "frameshift_deletion" = "frameshift",
  "frameshift_insertion" = "frameshift",
  "nonframeshift_deletion" = "non-frameshift",
  "nonframeshift_insertion" = "non-frameshift"
)

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
  # select only subset of genes
  order_gene_ocur = order_gene %>% 
    filter(n >= n_threshold & f >= f_threshold) %>% # any of the threshold passed
    arrange(-f) %>%
    select(gene) 
  return(as.character(order_gene_ocur$gene))
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

# perform sorting
gene_order = oncoplot_gene_order(f_threshold = 0.02)
sample_order = oncoplot_sample_order(gene_order)

# reorder and subset onco.
onco$sample = factor(onco$sample, levels = sample_order)
onco.old = onco
onco = onco.old %>% filter(gene %in% gene_order)
onco$gene = factor(onco$gene, levels = rev(gene_order)) # to make high order on top of heatmap

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

# draw heatmap, the body part
Hp <- ggplot(onco, aes(x=sample,y=gene,fill=mutation)) + 
  geom_tile(width=1,height=1.2,colour="grey70") +
  scale_fill_manual(values = cbPalette, na.value="white") +
  mythm + 
  theme(panel.border=element_rect(size=1, color="black"),
        axis.text.y=element_text(face="italic", size=11),
        legend.position="none"
  )

# draw sample status, the top part
Tp <- ggplot(subset(onco.old,!is.na(mutation)), aes(x=sample, fill=mutation)) + geom_bar() +
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
maxHeight=grid::unit.pmax(gT$heights[2:5],gH$heights[2:5])
# let the heights and widths consistent
gT$widths[2:5] <- gH$widths[2:5] <- as.list(maxWidth)
gT$heights[2:5] <- gH$heights[2:5] <- as.list(maxHeight)

# arrange the grobs
# how to arrange the heights sizes ? 
gC <- arrangeGrob(gT,empty,gH,gR,ncol=2,nrow=2, widths=c(7,1), heights=c(1,7))
gC <- gtable_add_rows(gC, heights = unit(3,"line"), pos=-1)
gC <- gtable_add_grob(gC, leg, t=3, b=3, l=1, r=1)


#jpeg("oncoplot.jpeg",width=1200,height=2000)
grid.newpage()
grid.draw(gC)
#dev.off()

}



