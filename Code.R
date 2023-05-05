#https://bernatgel.github.io/karyoploter_tutorial//Examples/Rainfall/Rainfall.html

somatic.mutations <- read.table("Pancreas_raw_mutations_data.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
somatic.mutations <- setNames(somatic.mutations, c("sample", "mut.type", "chr", "start", "end", "ref", "alt", "origin"))
head(somatic.mutations)


library(regioneR)
somatic.mutations <- split(somatic.mutations, somatic.mutations$sample)
sm <- somatic.mutations[["APGI_1992"]]
sm.gr <- toGRanges(sm[,c("chr", "start", "end", "mut.type", "ref", "alt")])
seqlevelsStyle(sm.gr) <- "UCSC"
sm.gr


library(karyoploteR)
kp <- plotKaryotype(plot.type=4)
kpPlotRainfall(kp, data = sm.gr)


variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$alt)
kp <- plotKaryotype(plot.type=4)
kpPlotRainfall(kp, data = sm.gr, col=variant.colors)


pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)
kpAddMainTitle(kp, main="Somatic Mutations - APGI_1992", cex=1.2)
kpAxis(kp, ymax = 7, tick.pos = 1:7)
kpPlotRainfall(kp, data = sm.gr, col=variant.colors)
kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04)


kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)
kpAddMainTitle(kp, main="Somatic Mutations - APGI_1992", cex=1.2)
kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.7)
kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.7)
kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.7)
kpPlotDensity(kp, data = sm.gr, r0=0.72, r1=1)
kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)

kata.genes <- toGRanges(data.frame(chr=c("chr11"),
                        start=c(34460472,34642588),
                        end=c(34460561,34684834),
                        labels=c("CAT","EHF"),
                        stringsAsFactors = FALSE))

kp <- plotKaryotype(plot.type=4, plot.params = pp,
                    chromosomes="chr11")
kpAddMainTitle(kp, main="Somatic Mutations - APGI_1992 - Chromosome 11", cex=1.2)
kpAddBaseNumbers(kp)
kpPlotDensity(kp, data = sm.gr, window.size = 10e5, r0=0.62, r1=0.8)
kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.04, r0=0.62, r1=0.8)
kpPlotMarkers(kp, data = kata.genes, text.orientation = "horizontal", r1=1.1, line.color = "#AAAAAA")
kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.6)
kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.6)
kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.6)

