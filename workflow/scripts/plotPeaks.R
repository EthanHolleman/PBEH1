
library(ggplot2)
library(rentrez)
library(ggbio)
library(seqinr)
library(ape)
library(ggplot2)
library(grid)

# function to read GenBank file and create data frame
read_genbank <- function(file_path){
  gb = read.GenBank(file_path, as.character = TRUE)
  features = gb@features
  names(features) = sapply(features, function(x) x$qualifiers$label)
  data.frame(start = unlist(sapply(features, function(x) x$location$range[1])),
             end = unlist(sapply(features, function(x) x$location$range[2])),
             strand = unlist(sapply(features, function(x) x$location$strand)),
             type = unlist(sapply(features, function(x) x$type)),
             gene = names(features))
}

# function to create linear plasmid map
plot_linear_plasmid <- function(file_path, start = NULL, end = NULL){
  
  # read in GenBank file and create data frame
  gb_df = read_genbank(file_path)
  
  # crop data frame if start and end provided
  if (!is.null(start)) {
    gb_df = gb_df[gb_df$end >= start,]
  }
  if (!is.null(end)) {
    gb_df = gb_df[gb_df$start <= end,]
  }
  
  # calculate total length of plasmid
  max_len = max(gb_df$end)
  
  # create plot
  p = ggplot() +
    xlim(0, max_len) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_void() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  # add features
  for (i in 1:nrow(gb_df)){
    row = gb_df[i,]
    if (row$strand == 1) {
      arrow_type = "forward"
    } else {
      arrow_type = "reverse"
    }
    p = p + geom_gene_arrow(aes(xmin = row$start, xmax = row$end, y = 0, 
                                label = row$gene, arrow_type = arrow_type), 
                            height = 0.3, offset = 0.15, 
                            arrow_width = unit(0.15, "cm"))
  }
  
  return(p)
}




# peaks.path.pos <- '/home/ethollem/projects/PBEH1/workflow/output/peakCall/PBEH1_1/1/T7_init_VR_1/PEAKS_GENOME/PCB11_geneT7_INIT_VR_1_Pos_20_0.55_CG.PEAK.genome.bed'
# peaks.path.neg <- '/home/ethollem/projects/PBEH1/workflow/output/peakCall/PBEH1_1/1/T7_init_VR_1/PEAKS_GENOME/PCB11_geneT7_INIT_VR_1_Neg_20_0.55_GC.PEAK.genome.bed'

# peaks.pos.df <- read_tsv_with_headers(peaks.path.pos)
# peaks.pos.df <- read_tsv_with_headers(peaks.path.neg)

genbank.path <- '/home/ethollem/projects/PBEH1/notes/test.gb'

plt <- plot_linear_plasmid(genbank.path, 400, 2000)
ggsave('test.png', plt, dpi=300)