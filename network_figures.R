install.packages("ggraph")
library(ggraph)
library(igraph)
bac.otu.0_2.Control_postive.modularity <- edge.betweenness.community(bac.otu.0_2.Control_postive.1, weights = NULL)
modularity(bac.otu.0_2.Control_postive.modularity)


tiff('bac.otu.0_2.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.Control_postive.1)$width/5),
                                    vertex.color=V(bac.otu.0_2.Control_postive.1)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(bac.otu.0_2.Control_postive.1)$name,vertex.label.color="black",
                                    vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#b73616"),
                                    vertex.size=5,
                                    layout=layout_nicely(bac.otu.0_2.Control_postive.1),main="Bacteria(+) 0_2 Control")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.Control_postive.modularity),2)))
dev.off()

tiff('1bac.otu.0_2.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.Control_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.Control_postive.1),main="Bacteria(+) 0_2 Control")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.Control_postive.modularity),2)))
dev.off()

bac.otu.0_2.Control_negative.modularity <- edge.betweenness.community(bac.otu.0_2.Control_negative.1,weights = NULL)
modularity(bac.otu.0_2.Control_negative.modularity)

tiff('bac.otu.0_2.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.Control_negative.modularity,bac.otu.0_2.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.Control_negative.1),main="Bacteria(-) 0_2 Control")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.Control_negative.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.Control_negative.1),main="Bacteria(-) 0_2 Control")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.Control_negative.modularity),2)))

dev.off()
bac.otu.0_2.W_positive.modularity <- edge.betweenness.community(bac.otu.0_2.W_postive.1,weights = NULL)
modularity(bac.otu.0_2.W_positive.modularity)

tiff('bac.otu.0_2.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.W_positive.modularity,bac.otu.0_2.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.W_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.W_postive.1),main="Bacteria(+) 0_2 W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.W_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.W_postive.1),main="Bacteria(+) 0_2 W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.W_positive.modularity),2)))

dev.off()

bac.otu.0_2.W_negative.modularity <- edge.betweenness.community(bac.otu.0_2.W_negative.1,weights = NULL)
modularity(bac.otu.0_2.W_negative.modularity)

tiff('bac.otu.0_2.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.W_negative.modularity,bac.otu.0_2.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.W_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.W_negative.1),main="Bacteria(-) 0_2 W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.W_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.W_negative.1),main="Bacteria(-) 0_2 W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.W_negative.modularity),2)))

dev.off()
bac.otu.0_2.T2W_positive.modularity <- edge.betweenness.community(bac.otu.0_2.T2W_postive.1,weights = NULL)
modularity(bac.otu.0_2.T2W_positive.modularity)

tiff('bac.otu.0_2.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2W_positive.modularity,bac.otu.0_2.T2W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2W_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.T2W_postive.1),main="Bacteria(+) 0_2 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2W_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.T2W_postive.1),main="Bacteria(+) 0_2 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2W_positive.modularity),2)))

dev.off()

bac.otu.0_2.T2W_negative.modularity <- edge.betweenness.community(bac.otu.0_2.T2W_negative.1,weights = NULL)
modularity(bac.otu.0_2.T2W_negative.modularity)

tiff('bac.otu.0_2.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2W_negative.modularity,bac.otu.0_2.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.T2W_negative.1),main="Bacteria(-) 0_2 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.T2W_negative.1),main="Bacteria(-) 0_2 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2W_negative.modularity),2)))

dev.off()

bac.otu.0_2.T2_positive.modularity <- edge.betweenness.community(bac.otu.0_2.T2_postive.1,weights = NULL)
modularity(bac.otu.0_2.T2_positive.modularity)

tiff('bac.otu.0_2.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2_positive.modularity,bac.otu.0_2.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.T2_postive.1),main="Bacteria(+) 0_2 T2")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2_positive.modularity),2)))

dev.off()

tiff('1bac.otu.0_2.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.0_2.T2_postive.1),main="Bacteria(+) 0_2 T2")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2_positive.modularity),2)))

dev.off()


bac.otu.0_2.T2_negative.modularity <- edge.betweenness.community(bac.otu.0_2.T2_negative.1,weights = NULL)
modularity(bac.otu.0_2.T2_negative.modularity)

tiff('bac.otu.0_2.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2_negative.modularity,bac.otu.0_2.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.T2_negative.1),main="Bacteria(-) 0_2 T2")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2_negative.modularity),2)))

dev.off()
tiff('1bac.otu.0_2.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.0_2.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.0_2.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.0_2.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.0_2.T2_negative.1),main="Bacteria(-) 0_2 T2")
legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.0_2.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2_negative.modularity),2)))

dev.off()
# par(mfrow=c(2,1))
# plot(bac.otu.0_2.Control_postive.modularity,bac.otu.0_2.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.Control_postive.1)$width/5),
#      vertex.color=V(bac.otu.0_2.Control_postive.1)$vertex_color,vertex.frame.color="black",
#      vertex.label=NA,vertex.label.color="black",
#      vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
#      edge.color=c("#b73616"),
#      vertex.size=V(bac.otu.0_2.Control_postive.1)$size*500,
#      layout=layout_nicely(bac.otu.0_2.Control_postive.1),main="Bacteria(+) 0_2 Control")
# legend(x=1.2, y=1, unique(V(bac.otu.0_2.Control_postive.1)$Phylum), pch=21,
#        col="#777777", pt.bg=unique(V(bac.otu.0_2.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
# text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.Control_postive.modularity),2)))
# 
# plot(bac.otu.0_2.T2_negative.modularity,bac.otu.0_2.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.0_2.T2_negative.1)$width/5),
#      vertex.color=V(bac.otu.0_2.T2_negative.1)$vertex_color,vertex.frame.color="black",
#      vertex.label=NA,vertex.label.color="black",
#      vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
#      edge.color=c("#0b4d2b"),
#      vertex.size=V(bac.otu.0_2.T2_negative.1)$size*500,
#      layout=layout_with_dh(bac.otu.0_2.T2_negative.1),main="Bacteria(-) 0_2 T2")
# legend(x=1.2, y=1, unique(V(bac.otu.0_2.T2_negative.1)$Phylum), pch=21,
#        col="#777777", pt.bg=unique(V(bac.otu.0_2.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
# text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.0_2.T2_negative.modularity),2)))




################################2-10cm################################
bac.otu.2_10.Control_postive.modularity <- edge.betweenness.community(bac.otu.2_10.Control_postive.1, weights = NULL)
modularity(bac.otu.2_10.Control_postive.modularity)


tiff('bac.otu.2_10.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.Control_postive.modularity,bac.otu.2_10.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.Control_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.Control_postive.1),main="Bacteria(+) 2_10 Control")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.Control_postive.modularity),2)))
dev.off()

tiff('1bac.otu.2_10.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.Control_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.Control_postive.1),main="Bacteria(+) 2_10 Control")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.Control_postive.modularity),2)))
dev.off()


bac.otu.2_10.Control_negative.modularity <- edge.betweenness.community(bac.otu.2_10.Control_negative.1,weights = NULL)
modularity(bac.otu.2_10.Control_negative.modularity)

tiff('bac.otu.2_10.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.Control_negative.modularity,bac.otu.2_10.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.Control_negative.1),main="Bacteria(-) 2_10 Control")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.Control_negative.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.Control_negative.1),main="Bacteria(-) 2_10 Control")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.Control_negative.modularity),2)))

dev.off()

bac.otu.2_10.W_positive.modularity <- edge.betweenness.community(bac.otu.2_10.W_postive.1,weights = NULL)
modularity(bac.otu.2_10.W_positive.modularity)

tiff('bac.otu.2_10.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.W_positive.modularity,bac.otu.2_10.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.W_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.W_postive.1),main="Bacteria(+) 2_10 W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.W_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.W_postive.1),main="Bacteria(+) 2_10 W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.W_positive.modularity),2)))

dev.off()

bac.otu.2_10.W_negative.modularity <- edge.betweenness.community(bac.otu.2_10.W_negative.1,weights = NULL)
modularity(bac.otu.2_10.W_negative.modularity)

tiff('bac.otu.2_10.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.W_negative.modularity,bac.otu.2_10.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.W_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.W_negative.1),main="Bacteria(-) 2_10 W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.W_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.W_negative.1),main="Bacteria(-) 2_10 W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.W_negative.modularity),2)))

dev.off()
bac.otu.2_10.T2W_positive.modularity <- edge.betweenness.community(bac.otu.2_10.T2W_postive.1,weights = NULL)
modularity(bac.otu.2_10.T2W_positive.modularity)

tiff('bac.otu.2_10.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2W_positive.modularity,bac.otu.2_10.T2W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2W_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.T2W_postive.1),main="Bacteria(+) 2_10 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2W_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.T2W_postive.1),main="Bacteria(+) 2_10 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2W_positive.modularity),2)))

dev.off()


bac.otu.2_10.T2W_negative.modularity <- edge.betweenness.community(bac.otu.2_10.T2W_negative.1,weights = NULL)
modularity(bac.otu.2_10.T2W_negative.modularity)

tiff('bac.otu.2_10.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2W_negative.modularity,bac.otu.2_10.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.T2W_negative.1),main="Bacteria(-) 2_10 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.T2W_negative.1),main="Bacteria(-) 2_10 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2W_negative.modularity),2)))

dev.off()

bac.otu.2_10.T2_positive.modularity <- edge.betweenness.community(bac.otu.2_10.T2_postive.1,weights = NULL)
modularity(bac.otu.2_10.T2_positive.modularity)

tiff('bac.otu.2_10.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2_positive.modularity,bac.otu.2_10.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.T2_postive.1),main="Bacteria(+) 2_10 T2")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2_positive.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.2_10.T2_postive.1),main="Bacteria(+) 2_10 T2")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2_positive.modularity),2)))

dev.off()
bac.otu.2_10.T2_negative.modularity <- edge.betweenness.community(bac.otu.2_10.T2_negative.1,weights = NULL)
modularity(bac.otu.2_10.T2_negative.modularity)

tiff('bac.otu.2_10.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2_negative.modularity,bac.otu.2_10.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.T2_negative.1),main="Bacteria(-) 2_10 T2")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2_negative.modularity),2)))

dev.off()

tiff('1bac.otu.2_10.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.2_10.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.2_10.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.2_10.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.2_10.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.2_10.T2_negative.1),main="Bacteria(-) 2_10 T2")
legend(x=1.2, y=1, unique(V(bac.otu.2_10.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.2_10.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.2_10.T2_negative.modularity),2)))

dev.off()
####################20-30
bac.otu.20_30.Control_postive.modularity <- edge.betweenness.community(bac.otu.20_30.Control_postive.1, weights = NULL)
modularity(bac.otu.20_30.Control_postive.modularity)


tiff('bac.otu.20_30.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.Control_postive.modularity,bac.otu.20_30.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.Control_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.Control_postive.1),main="Bacteria(+) 20_30 Control")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.Control_postive.modularity),2)))
dev.off()

tiff('1bac.otu.20_30.Control_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.Control_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.Control_postive.1),main="Bacteria(+) 20_30 Control")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.Control_postive.modularity),2)))
dev.off()

bac.otu.20_30.Control_negative.modularity <- edge.betweenness.community(bac.otu.20_30.Control_negative.1,weights = NULL)
modularity(bac.otu.20_30.Control_negative.modularity)

tiff('bac.otu.20_30.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.Control_negative.modularity,bac.otu.20_30.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.Control_negative.1),main="Bacteria(-) 20_30 Control")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.Control_negative.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.Control_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.Control_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.Control_negative.1),main="Bacteria(-) 20_30 Control")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.Control_negative.modularity),2)))

dev.off()
bac.otu.20_30.W_positive.modularity <- edge.betweenness.community(bac.otu.20_30.W_postive.1,weights = NULL)
modularity(bac.otu.20_30.W_positive.modularity)

tiff('bac.otu.20_30.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.W_positive.modularity,bac.otu.20_30.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.W_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.W_postive.1),main="Bacteria(+) 20_30 W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.W_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.W_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.W_postive.1),main="Bacteria(+) 20_30 W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.W_positive.modularity),2)))

dev.off()

bac.otu.20_30.W_negative.modularity <- edge.betweenness.community(bac.otu.20_30.W_negative.1,weights = NULL)
modularity(bac.otu.20_30.W_negative.modularity)

tiff('bac.otu.20_30.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.W_negative.modularity,bac.otu.20_30.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.W_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.W_negative.1),main="Bacteria(-) 20_30 W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.W_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.W_negative.1),main="Bacteria(-) 20_30 W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.W_negative.modularity),2)))

dev.off()
bac.otu.20_30.T2W_positive.modularity <- edge.betweenness.community(bac.otu.20_30.T2W_postive.2,weights = NULL)
modularity(bac.otu.20_30.T2W_positive.modularity)

tiff('bac.otu.20_30.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2W_positive.modularity,bac.otu.20_30.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2W_postive.2)$width/5),
     vertex.color=V(bac.otu.20_30.T2W_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2W_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.T2W_postive.2),main="Bacteria(+) 20_30 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(bac.otu.20_30.T2W_positive.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.T2W_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2W_postive.2)$width/5),
     vertex.color=V(bac.otu.20_30.T2W_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2W_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.T2W_postive.2),main="Bacteria(+) 20_30 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(bac.otu.20_30.T2W_positive.modularity),2)))

dev.off()

bac.otu.20_30.T2W_negative.modularity <- edge.betweenness.community(bac.otu.20_30.T2W_negative.1,weights = NULL)
modularity(bac.otu.20_30.T2W_negative.modularity)

tiff('bac.otu.20_30.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2W_negative.modularity,bac.otu.20_30.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.T2W_negative.1),main="Bacteria(-) 20_30 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2W_negative.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.T2W_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2W_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.T2W_negative.1),main="Bacteria(-) 20_30 T2W")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2W_negative.modularity),2)))

dev.off()

bac.otu.20_30.T2_positive.modularity <- edge.betweenness.community(bac.otu.20_30.T2_postive.1,weights = NULL)
modularity(bac.otu.20_30.T2_positive.modularity)

tiff('bac.otu.20_30.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2_positive.modularity,bac.otu.20_30.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.T2_postive.1),main="Bacteria(+) 20_30 T2")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2_positive.modularity),2)))

dev.off()
tiff('1bac.otu.20_30.T2_postive.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2_postive.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(bac.otu.20_30.T2_postive.1),main="Bacteria(+) 20_30 T2")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2_positive.modularity),2)))

dev.off()

bac.otu.20_30.T2_negative.modularity <- edge.betweenness.community(bac.otu.20_30.T2_negative.1,weights = NULL)
modularity(bac.otu.20_30.T2_negative.modularity)

tiff('bac.otu.20_30.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2_negative.modularity,bac.otu.20_30.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.T2_negative.1),main="Bacteria(-) 20_30 T2")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2_negative.modularity),2)))

dev.off()

tiff('1bac.otu.20_30.T2_negative.tiff', units="in", width=10, height=8, res=300)
plot(bac.otu.20_30.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(bac.otu.20_30.T2_negative.1)$width/5),
     vertex.color=V(bac.otu.20_30.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(bac.otu.20_30.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.5,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(bac.otu.20_30.T2_negative.1),main="Bacteria(-) 20_30 T2")
legend(x=1.2, y=1, unique(V(bac.otu.20_30.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(bac.otu.20_30.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(bac.otu.20_30.T2_negative.modularity),2)))

dev.off()