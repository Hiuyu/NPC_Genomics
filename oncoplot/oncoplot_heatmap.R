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

# how to increase the spacing of legend key ?
legend_theme <- theme(
  legend.position="bottom",
  legend.direction="horizontal",
  legend.text=element_text(size=10),
  legend.title=element_blank(),
  legend.key.size=unit(0.4,"cm")
)

# set manual color
mutColor=c("lightblue","green","yellow","pink","red")
names(mutColor)=c(dictType)

# draw heatmap, the body part
Hp <- ggplot(onco,aes(x=sample,y=gene,fill=mutation)) + 
  geom_tile(width=0.9,height=1,colour="white") +
  scale_fill_manual(values = mutColor, na.value="grey95") +
  mythm + 
  theme(panel.border=element_rect(size=1, color="black"),
        axis.text.y=element_text(face="italic", size=11),
        legend.position="none"
  )

# draw sample status, the top part
Tp <- ggplot(subset(onco,!is.na(mutation)), aes(x=sample, fill=mutation)) + geom_bar() +
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
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
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
  mythm + theme(legend.position="none", 
                rect=element_blank(),
                axis.line.x=element_blank(),
                axis.text.y=element_blank()
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
gC <- gtable_add_rows(gC, heights = unit(0.1,"npc"), pos=-1)
gC <- gtable_add_grob(gC, leg, t=3, b=3, l=1.5, r=1.5)


jpeg("C:/Users/LibJu/workspace/Github-project_scripts-repository/NPC_Genomics/oncoplot/demo_1.jpeg",width=600,height=800)
grid.newpage()
grid.draw(gC)
dev.off()





