## demo for heatmap with barplot on top and right panels, 
## also called the rainfall plot or oncoplot
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(reshape2)

# a demo mutation data
# set dictionary of mutation type
dictType=c("missense", 
           "nonsense",
           "silent",
           "frameshift",
           "nonframeshift"
           )
names(dictType)=1:5

# simulate a mutation data
# simulate a oncomatrix
df=data.frame(matrix(
  sample(0:5,100*38,replace = T,prob=c(0.9,0.05,0.02,0.003,0.07,0.08)),
  nr=100,nc=38)
)
colnames(df)=paste0("g",1:ncol(df))
df$sample=sample(paste0("s",1:nrow(df)))

# convert it to data.frame format
onco=melt(df,id.vars = "sample", variable.name = "gene", value.name = "mutation")
onco$mutation[onco$mutation==0]=NA # remove non-mutation rows
onco$mutation=factor(onco$mutation,labels = c(dictType)) # assign label

## how to decide gene order ?
# or just provide a pre-sort gene list?



## how to sort sample order ?
# or just provide a pre-sort sample list?
# or sort by frequency ?



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

# set manual color
mutColor=c("lightblue","green","yellow","pink","red")
names(mutColor)=c(dictType)

# draw heatmap, the body part
Hp <- ggplot(onco,aes(x=sample,y=gene,fill=mutation)) + 
  geom_tile(width=0.9,height=1) +
  scale_fill_manual(values = mutColor, na.value="grey95") +
  mythm + theme(legend.position="bottom",
                legend.direction="horizontal", 
                legend.title=element_blank(),
                legend.key.size=unit(0.01,"npc"),
                panel.border=element_rect(size=1.1, color="black"),
                axis.text.y=element_text(family="serif", face="bold.italic", size=13)
                )

# extract the legend
# a function to extract legend
extract_ggplot2_legend<-function(ggplot_obj) {
  tmp <- ggplot_gtable(ggplot_build(ggplot_obj))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

leg <- extract_ggplot2_legend(Hp)
Hp <- Hp + theme(legend.position="none") # remove legend after extracted

# draw sample status, the top part
Tp <- ggplot(subset(onco,!is.na(mutation)), aes(x=sample, fill=mutation)) + geom_bar() +
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
  mythm + theme(legend.position="none", 
        rect=element_blank(), 
        axis.line.y=element_line(size=1.1), 
        axis.line.x=element_line(size=1.1),
        axis.ticks.y=element_line(size=1.1,color="black"),
        # how to add tick on Y-axis ?
       # axis.ticks.length=unit(1,"cm"),
        axis.text.y=element_text(family="serif", face="bold", size=13)
        )

# draw gene status, the right part 
# how to add a top y-axis line and tick marker? use secondary axis?
Rp <- ggplot(subset(onco,!is.na(mutation)), aes(x=gene, fill=mutation)) + geom_bar() +
  coord_flip() +  # note that the cols order should be reverse
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
  mythm + theme(legend.position="none", 
                rect=element_blank(),
                axis.line.x=element_blank()
                axis.text.y=element_blank(),
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
gT=ggplotGrob(Tp) # grob all elements
gH=ggplotGrob(Hp)
gR=ggplotGrob(Rp)

# let the top and heatmap with same width, right and heatmap with same height
maxWidth=grid::unit.pmax(gT$width[2:5],gH$widths[2:5])
maxHeight=grid::unit.pmax(gT$heights[2:5],gH$heights[2:5])
# let the heights and widths consistent
gT$widths[2:5] <- gH$widths[2:5] <- as.list(maxWidth)
gT$heights[2:5] <- gH$heights[2:5] <- as.list(maxHeight)

# arrange the grobs
gC <- arrangeGrob(gT,empty,gH,gR,ncol=2,nrow=2, widths=c(4,1), heights=c(1,4))
gC <- gtable_add_rows(gC, heights = unit(0.1,"npc"), pos=-1)
gC <- gtable_add_grob(gC, leg$grobs, t=3, b=3, l=1.5, r=1.5)

grid.newpage()
grid.draw(gC)






