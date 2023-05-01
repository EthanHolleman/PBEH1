library(ggplot2)
library(ggpubr)
library(dplyr)


# helper functions 

remove_header_rows <- function(df) {
  # Get the header row of the data frame
    header <- colnames(df)
    non_header_rows <- apply(df, 1, function(x) !all(x == header))
    df[non_header_rows,]

}

# read input files from snakemake 


barcodes.df.path <- snakemake@input[['barcodes']]
bar.df <- as.data.frame(
    read.table(file = barcodes.df.path, sep = '\t', header = TRUE)
)
bar.df <- remove_header_rows(bar.df)
print(head(bar.df))


read.len.path <- snakemake@input[['read_lengths']]
len.df <- as.data.frame(
    read.table(file = read.len.path, sep = '\t', header = FALSE)
)

colnames(len.df) <- c('ref_name', 'read', 'length')
print(head(len.df))


expected.read.lengths <- snakemake@input[['exp_read_lengths']]
exp.len <- as.data.frame(
    read.table(file = expected.read.lengths, sep = ',', header = TRUE)
)


# merge in read length value
bar.df <- merge(bar.df, len.df, by='read')
print(head(bar.df))


barcode.bar.plt <- function(df){

    # df <- df %>% count(fwd_primer, rev_primer)
    # ggplot(df, aes(x=))

}


# Plot number of reads assigned to each barcode as heatmap type plots and
# include heatmaps that are filtered by good sample identification 

bar.df.sid <- subset(bar.df, sample_id != -1)  # samples that do have id
bar.df.no.sid <- subset(bar.df, sample_id == -1)  # sample without sample id


heatmap.plt <- function(df, title){

    # aggregate data
    df <- df %>% count(fwd_primer, rev_primer)
    print(head(df))

    ggplot(df, aes(x=fwd_primer, y=rev_primer, fill=n), color='black') + 
                geom_tile(color='black') + 
                scale_fill_gradient(low = "white", high = "red") +
                labs(x='Forward primer', y='Rev primer', title=title) +
                theme_pubr() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(legend.key.size = unit(4, 'cm'))


}


barcode.tile.all <- heatmap.plt(bar.df, 'Barcodes all reads')
barcode.tile.sid <- heatmap.plt(bar.df.sid, 'Barcodes reads with SID')
barcode.tile.no.sid <- heatmap.plt(bar.df.no.sid, 'Barcodes reads no SID')


bar.heatmaps <- ggarrange(
    barcode.tile.all, barcode.tile.sid , barcode.tile.no.sid,
    nrow=3, ncol=1
)

ggsave(
    snakemake@output[['heatmaps']], bar.heatmaps, width=10, height=30, dpi=300,
    unit='in'
)

# Plot the read length distribution for reads with fwd and rev barcodes
# as well as 1 or the other

bar.both <- subset(bar.df, fwd_primer!='none' & rev_primer!='none')
bar.fwd.only <- subset(bar.df, fwd_primer!='none' & rev_primer=='none')
bar.rev.only <- subset(bar.df, fwd_primer=='none' & rev_primer!='none')


read.len.dist.plt <- function(df, exp.len, fill, title){

    df$length <- as.numeric(df$length)
    print('===============================================================')
    print(head(df))

    ggplot(df, aes(x=length), fill=fill, color='black') +
        geom_density() + theme_pubr() + labs(x='Read length', title=title) +
        geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
                           linetype="dashed")


}

len.both.plt <- read.len.dist.plt(
    bar.both, exp.len, 'tan', 'Read lengths reads fwd and rev primers')
len.both.fwd <- read.len.dist.plt(
    bar.fwd.only, exp.len, 'tan', 'Read lengths reads fwd primers only')
len.both.rev <- read.len.dist.plt(
    bar.rev.only, exp.len, 'tan', 'Read lengths reads rev primers only')

len.dist.plts <- ggarrange(
    len.both.plt, len.both.fwd, len.both.rev,
    nrow=3, ncol=1
)

ggsave(
    snakemake@output[['length_dist']], len.dist.plts, width=10, height=30, dpi=300,
    unit='in'
)


# make barplot just showing counts of reads with and without barcodes
both <- nrow(bar.both)
fwd.only <- nrow(bar.fwd.only)
rev.only <- nrow(bar.rev.only)

counts.df <- data.frame(
    barcodes <- c('both', 'fwd', 'rev'),
    count <- c(both, fwd.only, rev.only)
)

counts.plt <- ggplot(counts.df, aes(x=barcodes, y=count, fill=barcodes), color='black') +
      geom_bar(stat='identity') + labs(x='Barcode alignment', y='Read count') +
      theme_pubr()
    
ggsave(
    snakemake@output[['counts']], counts.plt, width=5, height=7, dpi=300
)

    








