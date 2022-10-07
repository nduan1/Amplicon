install.packages("ggraph")
library(ggraph)
library(igraph)
its.otu.0_2.Control_postive.modularity <- edge.betweenness.community(its.otu.0_2.Control_postive.1, weights = NULL)
modularity(its.otu.0_2.Control_postive.modularity)

its.otu.0_2.Control_postive.plot <-ggraph(its.otu.0_2.Control_postive.1, layout="linear") +
  geom_edge_fan(color="gray50", width=0.8, alpha=0.5) +
  geom_node_point(color=V(its.otu.0_2.Control_postive.1)$vertex_color, aes(size = V(its.otu.0_2.Control_postive.1)$size*500)) +
  theme_void()

tiff('its.otu.0_2.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.Control_postive.modularity,its.otu.0_2.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.Control_postive.1)$width/5),
                                         vertex.color=V(its.otu.0_2.Control_postive.1)$vertex_color,vertex.frame.color="black",
                                         vertex.label=V(its.otu.0_2.Control_postive.1)$name,vertex.label.color="black",
                                         vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                         edge.color=c("#b73616"),
                                         vertex.size=degree(its.otu.0_2.Control_postive.1),
                                         layout=layout_nicely(its.otu.0_2.Control_postive.1),main="Fungi(+) 0_2 Control")
legend(x=1.2, y=1, unique(V(its.otu.0_2.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.Control_postive.modularity),2)))
dev.off()

tiff('1its.otu.0_2.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.Control_postive.1)$width/5),
     vertex.color=V(its.otu.0_2.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.0_2.Control_postive.1),main="Fungi(+) 0_2 Control")
legend(x=1.2, y=1, unique(V(its.otu.0_2.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.Control_postive.modularity),2)))
dev.off()

its.otu.0_2.Control_negative.modularity <- edge.betweenness.community(its.otu.0_2.Control_negative.1,weights = NULL)
modularity(its.otu.0_2.Control_negative.modularity)

tiff('its.otu.0_2.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.Control_negative.modularity,its.otu.0_2.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.Control_negative.1)$width/5),
                                         vertex.color=V(its.otu.0_2.Control_negative.1)$vertex_color,vertex.frame.color="black",
                                         vertex.label=V(its.otu.0_2.Control_negative.1)$name,vertex.label.color="black",
                                         vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                         edge.color=c("#0b4d2b"),
                                         vertex.size=degree(its.otu.0_2.Control_negative.1),
                                         layout=layout_with_dh(its.otu.0_2.Control_negative.1),main="Fungi(-) 0_2 Control")
legend(x=1.2, y=1, unique(V(its.otu.0_2.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.Control_negative.modularity),2)))

dev.off()

tiff('1its.otu.0_2.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.Control_negative.1)$width/5),
     vertex.color=V(its.otu.0_2.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.0_2.Control_negative.1),main="Fungi(-) 0_2 Control")
legend(x=1.2, y=1, unique(V(its.otu.0_2.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.Control_negative.modularity),2)))

dev.off()
its.otu.0_2.W_positive.modularity <- edge.betweenness.community(its.otu.0_2.W_postive.1,weights = NULL)
modularity(its.otu.0_2.W_positive.modularity)

tiff('its.otu.0_2.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.W_positive.modularity,its.otu.0_2.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.W_postive.1)$width/5),
                                   vertex.color=V(its.otu.0_2.W_postive.1)$vertex_color,vertex.frame.color="black",
                                   vertex.label=V(its.otu.0_2.W_postive.1)$name,vertex.label.color="black",
                                   vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                   edge.color=c("#b73616"),
                                   vertex.size=degree(its.otu.0_2.W_postive.1),
                                   layout=layout_nicely(its.otu.0_2.W_postive.1),main="Fungi(+) 0_2 W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.W_positive.modularity),2)))

dev.off()

tiff('1its.otu.0_2.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.W_postive.1)$width/5),
     vertex.color=V(its.otu.0_2.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.0_2.W_postive.1),main="Fungi(+) 0_2 W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.W_positive.modularity),2)))

dev.off()

its.otu.0_2.W_negative.modularity <- edge.betweenness.community(its.otu.0_2.W_negative.1,weights = NULL)
modularity(its.otu.0_2.W_negative.modularity)

tiff('its.otu.0_2.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.W_negative.modularity,its.otu.0_2.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.W_negative.1)$width/5),
                                   vertex.color=V(its.otu.0_2.W_negative.1)$vertex_color,vertex.frame.color="black",
                                   vertex.label=V(its.otu.0_2.W_negative.1)$name,vertex.label.color="black",
                                   vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                   edge.color=c("#0b4d2b"),
                                   vertex.size=degree(its.otu.0_2.W_negative.1),
                                   layout=layout_with_dh(its.otu.0_2.W_negative.1),main="Fungi(-) 0_2 W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.W_negative.modularity),2)))

dev.off()

tiff('1its.otu.0_2.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.W_negative.1)$width/5),
     vertex.color=V(its.otu.0_2.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.0_2.W_negative.1),main="Fungi(-) 0_2 W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.W_negative.modularity),2)))

dev.off()
its.otu.0_2.T2W_positive.modularity <- edge.betweenness.community(its.otu.0_2.T2W_postive.2,weights = NULL)
modularity(its.otu.0_2.T2W_positive.modularity)

tiff('its.otu.0_2.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2W_positive.modularity,its.otu.0_2.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2W_postive.2)$width/5),
                                   vertex.color=V(its.otu.0_2.T2W_postive.2)$vertex_color,vertex.frame.color="black",
                                   vertex.label=V(its.otu.0_2.T2W_postive.2)$name,vertex.label.color="black",
                                   vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                   edge.color=c("#b73616"),
                                   vertex.size=degree(its.otu.0_2.T2W_postive.2),
                                   layout=layout_nicely(its.otu.0_2.T2W_postive.2),main="Fungi(+) 0_2 T2W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.0_2.T2W_positive.modularity),2)))

dev.off()

tiff('1its.otu.0_2.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2W_postive.2)$width/5),
     vertex.color=V(its.otu.0_2.T2W_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.T2W_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.0_2.T2W_postive.2),main="Fungi(+) 0_2 T2W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.0_2.T2W_positive.modularity),2)))

dev.off()
its.otu.0_2.T2W_negative.modularity <- edge.betweenness.community(its.otu.0_2.T2W_negative.1,weights = NULL)
modularity(its.otu.0_2.T2W_negative.modularity)

tiff('its.otu.0_2.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2W_negative.modularity,its.otu.0_2.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2W_negative.1)$width/5),
                                      vertex.color=V(its.otu.0_2.T2W_negative.1)$vertex_color,vertex.frame.color="black",
                                      vertex.label=V(its.otu.0_2.T2W_negative.1)$name,vertex.label.color="black",
                                      vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                      edge.color=c("#0b4d2b"),
                                      vertex.size=degree(its.otu.0_2.T2W_negative.1),
                                      layout=layout_with_dh(its.otu.0_2.T2W_negative.1),main="Fungi(-) 0_2 T2W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2W_negative.modularity),2)))

dev.off()

tiff('1its.otu.0_2.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2W_negative.1)$width/5),
     vertex.color=V(its.otu.0_2.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.0_2.T2W_negative.1),main="Fungi(-) 0_2 T2W")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2W_negative.modularity),2)))

dev.off()

its.otu.0_2.T2_positive.modularity <- edge.betweenness.community(its.otu.0_2.T2_postive.1,weights = NULL)
modularity(its.otu.0_2.T2_positive.modularity)

tiff('its.otu.0_2.T2_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2_positive.modularity,its.otu.0_2.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2_postive.1)$width/5),
                                   vertex.color=V(its.otu.0_2.T2_postive.1)$vertex_color,vertex.frame.color="black",
                                   vertex.label=V(its.otu.0_2.T2_postive.1)$name,vertex.label.color="black",
                                   vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                   edge.color=c("#b73616"),
                                   vertex.size=degree(its.otu.0_2.T2_postive.1),
                                   layout=layout_nicely(its.otu.0_2.T2_postive.1),main="Fungi(+) 0_2 T2")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2_positive.modularity),2)))

dev.off()

tiff('1its.otu.0_2.T2_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2_postive.1)$width/5),
     vertex.color=V(its.otu.0_2.T2_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.T2_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.0_2.T2_postive.1),main="Fungi(+) 0_2 T2")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2_positive.modularity),2)))

dev.off()

its.otu.0_2.T2_negative.modularity <- edge.betweenness.community(its.otu.0_2.T2_negative.1,weights = NULL)
modularity(its.otu.0_2.T2_negative.modularity)

tiff('its.otu.0_2.T2_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2_negative.modularity,its.otu.0_2.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2_negative.1)$width/5),
                                    vertex.color=V(its.otu.0_2.T2_negative.1)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(its.otu.0_2.T2_negative.1)$name,vertex.label.color="black",
                                    vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#0b4d2b"),
                                    vertex.size=degree(its.otu.0_2.T2_negative.1),
                                    layout=layout_with_dh(its.otu.0_2.T2_negative.1),main="Fungi(-) 0_2 T2")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2_negative.modularity),2)))

dev.off()

tiff('1its.otu.0_2.T2_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.0_2.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.0_2.T2_negative.1)$width/5),
     vertex.color=V(its.otu.0_2.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.0_2.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.0_2.T2_negative.1),main="Fungi(-) 0_2 T2")
legend(x=1.2, y=1, unique(V(its.otu.0_2.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.0_2.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.0_2.T2_negative.modularity),2)))

dev.off()
library(ggpubr)
#put it together
tiff('its.otu.0_2.network.tiff', units="in", width=15, height=7, res=300)
ggarrange(its.otu.0_2.Control_postive.plot,its.otu.0_2.W_postive.plot,its.otu.0_2.T2W_postive.plot,its.otu.0_2.T2_postive.plot,
          its.otu.0_2.Control_negative.plot,its.otu.0_2.W_negative.plot,its.otu.0_2.T2W_negative.plot,its.otu.0_2.T2_negative.plot,
          labels = c("A", "B", "C","D","E","F","G","H"),ncol = 4, nrow = 2,common.legend = TRUE,legend = "bottom",heights = c(1,1))
dev.off()

################################2-10cm################################
its.otu.2_10.Control_postive.modularity <- edge.betweenness.community(its.otu.2_10.Control_postive.1, weights = NULL)
modularity(its.otu.2_10.Control_postive.modularity)

its.otu.2_10.Control_postive.plot <-ggraph(its.otu.2_10.Control_postive.1, layout="linear") +
  geom_edge_fan(color="gray50", width=0.8, alpha=0.5) +
  geom_node_point(color=V(its.otu.2_10.Control_postive.1)$vertex_color, aes(size = V(its.otu.2_10.Control_postive.1)$size*500)) +
  theme_void()

tiff('its.otu.2_10.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.Control_postive.modularity,its.otu.2_10.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.Control_postive.1)$width/5),
                                          vertex.color=V(its.otu.2_10.Control_postive.1)$vertex_color,vertex.frame.color="black",
                                          vertex.label=V(its.otu.2_10.Control_postive.1)$name,vertex.label.color="black",
                                          vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                          edge.color=c("#b73616"),
                                          vertex.size=degree(its.otu.2_10.Control_postive.1),
                                          layout=layout_nicely(its.otu.2_10.Control_postive.1),main="Fungi(+) 2_10 Control")
legend(x=1.2, y=1, unique(V(its.otu.2_10.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.Control_postive.modularity),2)))
dev.off()

tiff('1its.otu.2_10.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.Control_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.Control_postive.1)$width/5),
     vertex.color=V(its.otu.2_10.Control_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.Control_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.2_10.Control_postive.1),main="Fungi(+) 2_10 Control")
legend(x=1.2, y=1, unique(V(its.otu.2_10.Control_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.Control_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.Control_postive.modularity),2)))
dev.off()
its.otu.2_10.Control_negative.modularity <- edge.betweenness.community(its.otu.2_10.Control_negative.1,weights = NULL)
modularity(its.otu.2_10.Control_negative.modularity)

tiff('its.otu.2_10.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.Control_negative.modularity,its.otu.2_10.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.Control_negative.1)$width/5),
                                          vertex.color=V(its.otu.2_10.Control_negative.1)$vertex_color,vertex.frame.color="black",
                                          vertex.label=V(its.otu.2_10.Control_negative.1)$name,vertex.label.color="black",
                                          vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                          edge.color=c("#0b4d2b"),
                                          vertex.size=degree(its.otu.2_10.Control_negative.1),
                                          layout=layout_with_dh(its.otu.2_10.Control_negative.1),main="Fungi(-) 2_10 Control")
legend(x=1.2, y=1, unique(V(its.otu.2_10.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.Control_negative.modularity),2)))

dev.off()

tiff('1its.otu.2_10.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.Control_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.Control_negative.1)$width/5),
     vertex.color=V(its.otu.2_10.Control_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.Control_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.2_10.Control_negative.1),main="Fungi(-) 2_10 Control")
legend(x=1.2, y=1, unique(V(its.otu.2_10.Control_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.Control_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.Control_negative.modularity),2)))

dev.off()

its.otu.2_10.W_positive.modularity <- edge.betweenness.community(its.otu.2_10.W_postive.1,weights = NULL)
modularity(its.otu.2_10.W_positive.modularity)

tiff('its.otu.2_10.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.W_positive.modularity,its.otu.2_10.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.W_postive.1)$width/5),
                                    vertex.color=V(its.otu.2_10.W_postive.1)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(its.otu.2_10.W_postive.1)$name,vertex.label.color="black",
                                    vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#b73616"),
                                    vertex.size=degree(its.otu.2_10.W_postive.1),
                                    layout=layout_nicely(its.otu.2_10.W_postive.1),main="Fungi(+) 2_10 W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.W_positive.modularity),2)))

dev.off()

tiff('1its.otu.2_10.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.W_postive.1)$width/5),
     vertex.color=V(its.otu.2_10.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.2_10.W_postive.1),main="Fungi(+) 2_10 W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.W_positive.modularity),2)))

dev.off()

its.otu.2_10.W_negative.modularity <- edge.betweenness.community(its.otu.2_10.W_negative.1,weights = NULL)
modularity(its.otu.2_10.W_negative.modularity)

tiff('its.otu.2_10.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.W_negative.modularity,its.otu.2_10.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.W_negative.1)$width/5),
                                    vertex.color=V(its.otu.2_10.W_negative.1)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(its.otu.2_10.W_negative.1)$name,vertex.label.color="black",
                                    vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#0b4d2b"),
                                    vertex.size=degree(its.otu.2_10.W_negative.1),
                                    layout=layout_with_dh(its.otu.2_10.W_negative.1),main="Fungi(-) 2_10 W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.W_negative.modularity),2)))

dev.off()

tiff('1its.otu.2_10.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.W_negative.1)$width/5),
     vertex.color=V(its.otu.2_10.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.2_10.W_negative.1),main="Fungi(-) 2_10 W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.W_negative.modularity),2)))

dev.off()

its.otu.2_10.T2W_positive.modularity <- edge.betweenness.community(its.otu.2_10.T2W_postive.2,weights = NULL)
modularity(its.otu.2_10.T2W_positive.modularity)

tiff('its.otu.2_10.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2W_positive.modularity,its.otu.2_10.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2W_postive.2)$width/5),
                                    vertex.color=V(its.otu.2_10.T2W_postive.2)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(its.otu.2_10.T2W_postive.2)$name,vertex.label.color="black",
                                    vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#b73616"),
                                    vertex.size=degree(its.otu.2_10.T2W_postive.2),
                                    layout=layout_nicely(its.otu.2_10.T2W_postive.2),main="Fungi(+) 2_10 T2W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.2_10.T2W_positive.modularity),2)))

dev.off()

tiff('1its.otu.2_10.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2W_postive.2)$width/5),
     vertex.color=V(its.otu.2_10.T2W_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.T2W_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.2_10.T2W_postive.2),main="Fungi(+) 2_10 T2W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.2_10.T2W_positive.modularity),2)))

dev.off()

its.otu.2_10.T2W_negative.modularity <- edge.betweenness.community(its.otu.2_10.T2W_negative.1,weights = NULL)
modularity(its.otu.2_10.T2W_negative.modularity)

tiff('its.otu.2_10.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2W_negative.modularity,its.otu.2_10.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2W_negative.1)$width/5),
                                       vertex.color=V(its.otu.2_10.T2W_negative.1)$vertex_color,vertex.frame.color="black",
                                       vertex.label=V(its.otu.2_10.T2W_negative.1)$name,vertex.label.color="black",
                                       vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                       edge.color=c("#0b4d2b"),
                                       vertex.size=degree(its.otu.2_10.T2W_negative.1),
                                       layout=layout_with_dh(its.otu.2_10.T2W_negative.1),main="Fungi(-) 2_10 T2W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.T2W_negative.modularity),2)))

dev.off()

tiff('1its.otu.2_10.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2W_negative.1)$width/5),
     vertex.color=V(its.otu.2_10.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.2_10.T2W_negative.1),main="Fungi(-) 2_10 T2W")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.T2W_negative.modularity),2)))

dev.off()

its.otu.2_10.T2_positive.modularity <- edge.betweenness.community(its.otu.2_10.T2_postive.2,weights = NULL)
modularity(its.otu.2_10.T2_positive.modularity)

tiff('its.otu.2_10.T2_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2_positive.modularity,its.otu.2_10.T2_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2_postive.2)$width/5),
                                    vertex.color=V(its.otu.2_10.T2_postive.2)$vertex_color,vertex.frame.color="black",
                                    vertex.label=V(its.otu.2_10.T2_postive.2)$name,vertex.label.color="black",
                                    vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                    edge.color=c("#b73616"),
                                    vertex.size=degree(its.otu.2_10.T2_postive.2),
                                    layout=layout_nicely(its.otu.2_10.T2_postive.2),main="Fungi(+) 2_10 T2")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.2_10.T2_positive.modularity),2)))

dev.off()

tiff('1its.otu.2_10.T2_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2_postive.2)$width/5),
     vertex.color=V(its.otu.2_10.T2_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.T2_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.2_10.T2_postive.2),main="Fungi(+) 2_10 T2")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.2_10.T2_positive.modularity),2)))

dev.off()

its.otu.2_10.T2_negative.modularity <- edge.betweenness.community(its.otu.2_10.T2_negative.1,weights = NULL)
modularity(its.otu.2_10.T2_negative.modularity)

tiff('its.otu.2_10.T2_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2_negative.modularity,its.otu.2_10.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2_negative.1)$width/5),
                                     vertex.color=V(its.otu.2_10.T2_negative.1)$vertex_color,vertex.frame.color="black",
                                     vertex.label=V(its.otu.2_10.T2_negative.1)$name,vertex.label.color="black",
                                     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                     edge.color=c("#0b4d2b"),
                                     vertex.size=degree(its.otu.2_10.T2_negative.1),
                                     layout=layout_with_dh(its.otu.2_10.T2_negative.1),main="Fungi(-) 2_10 T2")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.T2_negative.modularity),2)))

dev.off()

tiff('1its.otu.2_10.T2_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.2_10.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.2_10.T2_negative.1)$width/5),
     vertex.color=V(its.otu.2_10.T2_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.2_10.T2_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.2_10.T2_negative.1),main="Fungi(-) 2_10 T2")
legend(x=1.2, y=1, unique(V(its.otu.2_10.T2_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.2_10.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.2_10.T2_negative.modularity),2)))

dev.off()

library(ggpubr)
#put it together
tiff('its.otu.2_10.network.tiff', units="in", width=15, height=7, res=300)
ggarrange(its.otu.2_10.Control_postive.plot,its.otu.2_10.W_postive.plot,its.otu.2_10.T2W_postive.plot,its.otu.2_10.T2_postive.plot,
          its.otu.2_10.Control_negative.plot,its.otu.2_10.W_negative.plot,its.otu.2_10.T2W_negative.plot,its.otu.2_10.T2_negative.plot,
          labels = c("A", "B", "C","D","E","F","G","H"),ncol = 4, nrow = 2,common.legend = TRUE,legend = "bottom",heights = c(1,1))
dev.off()

##################################20-30###############################
its.otu.20_30.Control_postive.modularity <- edge.betweenness.community(its.otu.20_30.Control_postive.1, weights = NULL)
modularity(its.otu.20_30.Control_postive.modularity)

its.otu.20_30.Control_postive.plot <-ggraph(its.otu.20_30.Control_postive.1, layout="linear") +
  geom_edge_fan(color="gray50", width=0.8, alpha=0.5) +
  geom_node_point(color=V(its.otu.20_30.Control_postive.1)$vertex_color, aes(size = V(its.otu.20_30.Control_postive.1)$size*500)) +
  theme_void()

tiff('its.otu.20_30.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.Control_postive.modularity,its.otu.20_30.Control_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.Control_postive.2)$width/5),
                                           vertex.color=V(its.otu.20_30.Control_postive.2)$vertex_color,vertex.frame.color="black",
                                           vertex.label=V(its.otu.20_30.Control_postive.2)$name,vertex.label.color="black",
                                           vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                           edge.color=c("#b73616"),
                                           vertex.size=degree(its.otu.20_30.Control_postive.2),
                                           layout=layout_nicely(its.otu.20_30.Control_postive.2),main="Fungi(+) 20_30 Control")
legend(x=1.2, y=1, unique(V(its.otu.20_30.Control_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.Control_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.Control_postive.modularity),2)))
dev.off()

tiff('1its.otu.20_30.Control_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.Control_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.Control_postive.2)$width/5),
     vertex.color=V(its.otu.20_30.Control_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.Control_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.20_30.Control_postive.2),main="Fungi(+) 20_30 Control")
legend(x=1.2, y=1, unique(V(its.otu.20_30.Control_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.Control_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.Control_postive.modularity),2)))
dev.off()

its.otu.20_30.Control_negative.modularity <- edge.betweenness.community(its.otu.20_30.Control_negative.2,weights = NULL)
modularity(its.otu.20_30.Control_negative.modularity)

tiff('its.otu.20_30.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.Control_negative.modularity,its.otu.20_30.Control_negative.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.Control_negative.2)$width/5),
                                           vertex.color=V(its.otu.20_30.Control_negative.2)$vertex_color,vertex.frame.color="black",
                                           vertex.label=V(its.otu.20_30.Control_negative.2)$name,vertex.label.color="black",
                                           vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                           edge.color=c("#0b4d2b"),
                                           vertex.size=degree(its.otu.20_30.Control_negative.2),
                                           layout=layout_with_dh(its.otu.20_30.Control_negative.2),main="Fungi(-) 20_30 Control")
legend(x=1.2, y=1, unique(V(its.otu.20_30.Control_negative.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.Control_negative.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.Control_negative.modularity),2)))

dev.off()

tiff('1its.otu.20_30.Control_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.Control_negative.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.Control_negative.2)$width/5),
     vertex.color=V(its.otu.20_30.Control_negative.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.Control_negative.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.20_30.Control_negative.2),main="Fungi(-) 20_30 Control")
legend(x=1.2, y=1, unique(V(its.otu.20_30.Control_negative.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.Control_negative.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.Control_negative.modularity),2)))

dev.off()
its.otu.20_30.W_positive.modularity <- edge.betweenness.community(its.otu.20_30.W_postive.1,weights = NULL)
modularity(its.otu.20_30.W_positive.modularity)

tiff('its.otu.20_30.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.W_positive.modularity,its.otu.20_30.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.W_postive.1)$width/5),
                                     vertex.color=V(its.otu.20_30.W_postive.1)$vertex_color,vertex.frame.color="black",
                                     vertex.label=V(its.otu.20_30.W_postive.1)$name,vertex.label.color="black",
                                     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                     edge.color=c("#b73616"),
                                     vertex.size=degree(its.otu.20_30.W_postive.1),
                                     layout=layout_nicely(its.otu.20_30.W_postive.1),main="Fungi(+) 20_30 W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.W_positive.modularity),2)))

dev.off()

tiff('1its.otu.20_30.W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.W_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.W_postive.1)$width/5),
     vertex.color=V(its.otu.20_30.W_postive.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.W_postive.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.20_30.W_postive.1),main="Fungi(+) 20_30 W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.W_postive.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.W_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.W_positive.modularity),2)))

dev.off()

its.otu.20_30.W_negative.modularity <- edge.betweenness.community(its.otu.20_30.W_negative.1,weights = NULL)
modularity(its.otu.20_30.W_negative.modularity)

tiff('its.otu.20_30.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.W_negative.modularity,its.otu.20_30.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.W_negative.1)$width/5),
                                     vertex.color=V(its.otu.20_30.W_negative.1)$vertex_color,vertex.frame.color="black",
                                     vertex.label=V(its.otu.20_30.W_negative.1)$name,vertex.label.color="black",
                                     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                     edge.color=c("#0b4d2b"),
                                     vertex.size=degree(its.otu.20_30.W_negative.1),
                                     layout=layout_with_dh(its.otu.20_30.W_negative.1),main="Fungi(-) 20_30 W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.W_negative.modularity),2)))

dev.off()

tiff('1its.otu.20_30.W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.W_negative.1)$width/5),
     vertex.color=V(its.otu.20_30.W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.20_30.W_negative.1),main="Fungi(-) 20_30 W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.W_negative.modularity),2)))

dev.off()

its.otu.20_30.T2W_positive.modularity <- edge.betweenness.community(its.otu.20_30.T2W_postive.2,weights = NULL)
modularity(its.otu.20_30.T2W_positive.modularity)

tiff('its.otu.20_30.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.T2W_positive.modularity,its.otu.20_30.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2W_postive.2)$width/5),
                                     vertex.color=V(its.otu.20_30.T2W_postive.2)$vertex_color,vertex.frame.color="black",
                                     vertex.label=V(its.otu.20_30.T2W_postive.2)$name,vertex.label.color="black",
                                     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                     edge.color=c("#b73616"),
                                     vertex.size=degree(its.otu.20_30.T2W_postive.2),
                                     layout=layout_nicely(its.otu.20_30.T2W_postive.2),main="Fungi(+) 20_30 T2W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.T2W_positive.modularity),2)))

dev.off()

tiff('1its.otu.20_30.T2W_postive.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.T2W_postive.2,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2W_postive.2)$width/5),
     vertex.color=V(its.otu.20_30.T2W_postive.2)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.T2W_postive.2)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#b73616"),
     vertex.size=5,
     layout=layout_nicely(its.otu.20_30.T2W_postive.2),main="Fungi(+) 20_30 T2W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.T2W_postive.2)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.T2W_postive.2)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.2,paste("Modularity:",round(modularity(its.otu.20_30.T2W_positive.modularity),2)))

dev.off()
its.otu.20_30.T2W_negative.modularity <- edge.betweenness.community(its.otu.20_30.T2W_negative.1,weights = NULL)
modularity(its.otu.20_30.T2W_negative.modularity)

tiff('its.otu.20_30.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.T2W_negative.modularity,its.otu.20_30.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2W_negative.1)$width/5),
                                        vertex.color=V(its.otu.20_30.T2W_negative.1)$vertex_color,vertex.frame.color="black",
                                        vertex.label=V(its.otu.20_30.T2W_negative.1)$name,vertex.label.color="black",
                                        vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
                                        edge.color=c("#0b4d2b"),
                                        vertex.size=degree(its.otu.20_30.T2W_negative.1),
                                        layout=layout_with_dh(its.otu.20_30.T2W_negative.1),main="Fungi(-) 20_30 T2W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.T2W_negative.modularity),2)))

dev.off()

tiff('1its.otu.20_30.T2W_negative.tiff', units="in", width=10.5, height=8, res=300)
plot(its.otu.20_30.T2W_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2W_negative.1)$width/5),
     vertex.color=V(its.otu.20_30.T2W_negative.1)$vertex_color,vertex.frame.color="black",
     vertex.label=V(its.otu.20_30.T2W_negative.1)$name,vertex.label.color="black",
     vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
     edge.color=c("#0b4d2b"),
     vertex.size=5,
     layout=layout_with_dh(its.otu.20_30.T2W_negative.1),main="Fungi(-) 20_30 T2W")
legend(x=1.2, y=1, unique(V(its.otu.20_30.T2W_negative.1)$Phylum), pch=21,
       col="#777777", pt.bg=unique(V(its.otu.20_30.T2W_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.T2W_negative.modularity),2)))

dev.off()
# its.otu.20_30.T2_positive.modularity <- edge.betweenness.community(its.otu.20_30.T2_postive.1,weights = NULL)
# modularity(its.otu.20_30.T2_positive.modularity)
# 
# tiff('its.otu.20_30.T2_postive.tiff', units="in", width=10.5, height=8, res=300)
# plot(its.otu.20_30.T2_positive.modularity,its.otu.20_30.T2_postive.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2_postive.1)$width/5),
#                                      vertex.color=V(its.otu.20_30.T2_postive.1)$vertex_color,vertex.frame.color="black",
#                                      vertex.label=NA,vertex.label.color="black",
#                                      vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
#                                      edge.color=c("#b73616"),
#                                      vertex.size=degree(its.otu.20_30.T2_postive.1),
#                                      layout=layout_nicely(its.otu.20_30.T2_postive.1),main="Fungi(+) 20_30 T2")
# legend(x=1.2, y=1, unique(V(its.otu.20_30.T2_postive.1)$Phylum), pch=21,
#        col="#777777", pt.bg=unique(V(its.otu.20_30.T2_postive.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
# text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.T2_positive.modularity),2)))
# 
# dev.off()

# its.otu.20_30.T2_negative.modularity <- edge.betweenness.community(its.otu.20_30.T2_negative.1,weights = NULL)
# modularity(its.otu.20_30.T2_negative.modularity)
# 
# tiff('its.otu.20_30.T2_negative.tiff', units="in", width=10.5, height=8, res=300)
# its.otu.20_30.T2_negative.plot <-plot(its.otu.20_30.T2_negative.modularity,its.otu.20_30.T2_negative.1,edge.arrow.size=0,edge.width=abs(E(its.otu.20_30.T2_negative.1)$width/5),
#                                       vertex.color=V(its.otu.20_30.T2_negative.1)$vertex_color,vertex.frame.color="black",
#                                       vertex.label=NA,vertex.label.color="black",
#                                       vertex.label.cex=.8,vertex.label.font=2,vertex.label.dist=0,
#                                       edge.color=c("#0b4d2b"),
#                                       vertex.size=V(its.otu.20_30.T2_negative.1)$size*500,
#                                       layout=layout_with_dh(its.otu.20_30.T2_negative.1),main="Fungi(-) 20_30 T2")
# legend(x=1.2, y=1, unique(V(its.otu.20_30.T2_negative.1)$Phylum), pch=21,
#        col="#777777", pt.bg=unique(V(its.otu.20_30.T2_negative.1)$vertex_color), pt.cex=2, cex=.8, bty="n", ncol=1)
# text(x=1.45, y=1.1,paste("Modularity:",round(modularity(its.otu.20_30.T2_negative.modularity),2)))
# 
# dev.off()
# 
# library(ggpubr)
# #put it together
# tiff('its.otu.20_30.network.tiff', units="in", width=15, height=7, res=300)
# ggarrange(its.otu.20_30.Control_postive.plot,its.otu.20_30.W_postive.plot,its.otu.20_30.T2W_postive.plot,its.otu.20_30.T2_postive.plot,
#           its.otu.20_30.Control_negative.plot,its.otu.20_30.W_negative.plot,its.otu.20_30.T2W_negative.plot,its.otu.20_30.T2_negative.plot,
#           labels = c("A", "B", "C","D","E","F","G","H"),ncol = 4, nrow = 2,common.legend = TRUE,legend = "bottom",heights = c(1,1))

summary(its.otu.0_2.Control_postive.1)
summary(its.otu.0_2.Control_negative.1)
summary(its.otu.0_2.T2_postive.1)
summary(its.otu.0_2.T2_negative.1)
summary(its.otu.0_2.T2W_postive.1)
summary(its.otu.0_2.T2W_negative.1)
summary(its.otu.0_2.W_postive.1)
summary(its.otu.0_2.W_negative.1)

summary(its.otu.2_10.Control_postive.1)
summary(its.otu.2_10.Control_negative.1)
summary(its.otu.2_10.T2_postive.1)
summary(its.otu.2_10.T2_negative.1)
summary(its.otu.2_10.T2W_postive.1)
summary(its.otu.2_10.T2W_negative.1)
summary(its.otu.2_10.W_postive.1)
summary(its.otu.2_10.W_negative.1)

summary(its.otu.20_30.Control_postive.1)
summary(its.otu.20_30.Control_negative.1)
summary(its.otu.20_30.T2_postive.1)
summary(its.otu.20_30.T2_negative.1)
summary(its.otu.20_30.T2W_postive.1)
summary(its.otu.20_30.T2W_negative.1)
summary(its.otu.20_30.W_postive.1)
summary(its.otu.20_30.W_negative.1)






# dev.off()