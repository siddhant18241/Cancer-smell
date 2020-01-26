library(circlize)
#install.packages("circlize")
setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/Figures/Fig1_B_chordDiagram")
pdf(file="Fig1_B_TUMOR.pdf")
a2<-read.csv("chord_diagram_tumor.csv")
a2<-as.data.frame(a2)
#a2<-as.factor(a2)
chordDiagram(a2,transparency = 0.7, annotationTrack = "grid",preAllocateTracks = list(track.height = 0.5))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
sector.name = get.cell.meta.data("sector.index")
circos.text(mean(xlim), ylim[1],cex=0.3,sector.name, facing = "clockwise",
niceFacing = TRUE, adj = c(0, 0.6))
}, bg.border = NA)
dev.off()

#chord_diagram_CELLLINE.csv
a2<-read.csv("chord_diagram_CELLLINE.csv")
a2<-as.data.frame(a2)
#a2<-as.factor(a2)
pdf(file="Fig1_B_CELLLINE.pdf")
chordDiagram(a2,transparency = 0.7, annotationTrack = "grid",preAllocateTracks = list(track.height = 0.5))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1],cex=0.3,sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.6))
}, bg.border = NA)
dev.off()

