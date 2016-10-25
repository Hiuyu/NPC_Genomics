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




## how to decide sample order ?




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
  geom_tile(width=0.9,height=0.9) +
  scale_fill_manual(values = mutColor, na.value="grey95") +
  mythm

# draw sample status, the top part
Tp <- ggplot(subset(onco,!is.na(mutation)), aes(x=sample, fill=mutation)) + geom_bar() +
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
  mythm + 
  theme(legend.position="none", rect=element_blank(), axis.line.y=element_line(), axis.line.x=element_line())

# draw gene status, the right part 
Rp <- ggplot(subset(onco,!is.na(mutation)), aes(x=gene, fill=mutation)) + geom_bar() +
  coord_flip() +  # note that the cols order should be reverse
  scale_fill_manual(values=mutColor) + # set fill color
  scale_y_continuous(breaks=c(0,10,30),limits = c(0,30),expand = c(0, 0)) +  # set y axis
  mythm + 
  theme(legend.position="none", rect=element_blank(), axis.line.x=element_line(), axis.text.y=element_blank())











