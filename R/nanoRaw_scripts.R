#!/usr/bin/Rscript

# nanoR version 1.0
# Copyright (C) 2018 Felix Grünberger
#
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# supress warnings
options(warn=-1)

#' Genome coverage of nanopore run
#' Takes in the `sequencing_summary` file of albacore and the genome fasta as input and returns the genome coverage
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param quality Mean quality filter option for every read
#' @param genome_fasta path to genome fasta
#' @return Genome coverage = Nr of bases sequenced/genome length
#' @export
#' @import dplyr data.table seqinr
nano_seq_cov <- function(inputfile, genome_fasta, quality = 0){
  fasta <- seqinr::read.fasta(file = genome_fasta)
  names(fasta) <- "genome"
  return(fread(input = inputfile) %>%
           as.data.table() %>%
           dplyr::filter(mean_qscore_template >= quality) %>%
           summarise(sum(sequence_length_template)/length(fasta$genome))
  )
}


#' Histogram of read lengths of nanopore reads
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the readlength distribution as a histogram using ggplot2
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @param transformation "log10" transformation or leave it empty for regular plot
#' @return Histogram of read lengths
#' @export
#' @import ggplot2 dplyr viridis data.table
nano_histlength <- function(inputfile, viridiscolor = "viridis", transformation = ""){
  if (transformation == ""){
    return(fread(input = inputfile) %>%
      as.data.table() %>%
      ggplot(aes(x=sequence_length_template , fill = ..count..)) +
        geom_histogram(bins = 150, col = "white", lwd = 0.1) +
        theme_light() +
        xlab("Read length") +
        ylab("Number of reads") +
        ggtitle("nano_histlength: Histogram of read lengths") +
        guides(fill = F) +
        scale_fill_viridis(option = viridiscolor))
  }
  if(transformation == "log10"){
    return(fread(input = inputfile) %>%
             as.data.table() %>%
             ggplot(aes(x=sequence_length_template , fill = ..count..)) +
             geom_histogram(bins = 150, col = "white", lwd = 0.1) +
             theme_light() +
             xlab("Read length") +
             ylab("Number of reads") +
             ggtitle("nano_histlength: Histogram of reads lengths after log10 transformation") +
             guides(fill = F) +
             scale_x_log10() +
             scale_fill_viridis(option = viridiscolor))
  }
}




#' Bivariate plot of length against base call quality
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the readquality distribution as a histogram using ggplot2
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @return Bivariate plot of read lengths vs read quality
#' @export
#' @import ggplot2 dplyr viridis data.table ggExtra cowplot
nano_biplot <- function(inputfile, viridiscolor = "viridis"){
    seq_summary <- fread(input = inputfile) %>%
             as.data.table()
    pmain <- ggplot(data = seq_summary, aes(x= sequence_length_template, y = mean_qscore_template)) +
                stat_density2d(aes(alpha=..level.., fill = ..level..), geom="polygon") +
                theme_light() +
                scale_fill_viridis(option = viridiscolor) +
                xlab("Read lengths") +
                ylab("Average read quality") +
                ggtitle("nano_biplot: Read lengths vs Average read quality plot") +
                guides(fill = F, alpha = F)

    xdens <- axis_canvas(pmain, axis = "x") +
      geom_density(data = seq_summary, aes(x = sequence_length_template),
                   alpha = 0.3, size = 0.2, fill = viridis_pal(option = viridiscolor)(10)[1])


    ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
      geom_density(data = seq_summary, aes(x = mean_qscore_template),
                   alpha = 0.3, size = 0.2,fill = viridis_pal(option = viridiscolor)(10)[1]) +
      coord_flip()

    p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")

    p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

    return(ggdraw(p2))
}


#' Heatmap of reads per channel
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the reads generated per channel over time
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @return heatmap of reads per channel
#' @export
#' @import ggplot2 dplyr viridis data.table
nano_channel <- function(inputfile, viridiscolor = "viridis"){
  platelayout <- data.frame (rown = rep (letters[1:16], 32), coln = rep (1:32, each = 16), channel = 1:512)
  return(fread(input = inputfile) %>%
           as.data.table() %>%
    group_by(channel) %>%
    summarise(n = n()) %>%
    left_join(platelayout, by = "channel") %>%
    ggplot(aes(y = factor(rown),x = factor(coln))) +
      geom_tile(color  = "white", aes(fill = n),size = 0.25)  +
      scale_fill_viridis(option = viridiscolor, guide = guide_colorbar())+#,breaks=c(1,10,100,1000,8000)) +
      ggtitle("nano_channel: Number of reads generated by channel") +
      theme_light() +
      theme(legend.position="bottom") +
      theme(legend.title=element_blank()) +
      labs(x=NULL, y = NULL) +
      guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5, ticks = F)) +
      coord_fixed(ratio=1))
}


#' Cumulative yield plot
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the readquality distribution as a histogram using ggplot2
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @return Cumulative yield plot
#' @export
#' @import ggplot2 dplyr viridis data.table
nano_yield <- function(inputfile, viridiscolor = "viridis"){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
    mutate(timeindex = round((start_time + 3600)/3600)) %>%
    group_by(timeindex) %>%
    summarise(yield = sum(sequence_length_template)/1000000000) %>%
    ggplot(aes(x=timeindex, y=cumsum(yield))) +
      geom_area(alpha = 0.6, fill = viridis_pal(option = viridiscolor)(10)[3]) +
      geom_line(alpha = 0.2, col = "white", size = 10 ) +
      theme_light() +
      ggtitle("nano_yield") +
      xlab("Run time (hours)") +
      ylab("Cumulative yield in gigabase"))
}

#' joy distribution of read length
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the readquality distribution as a histogram using ggplot2
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @param transformation "log10" transformation or leave it empty for regular plot
#' @return joy distribution of read length
#' @export
#' @import ggplot2 dplyr viridis data.table ggjoy
nano_joy_length <- function(inputfile,  viridiscolor = "viridis", transformation = ""){
  if(transformation == ""){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
    mutate(timeindex = round((start_time + 3600)/3600)) %>%
    group_by(timeindex) %>%
    ggplot(aes(x = sequence_length_template, y = factor(timeindex), fill = ..x..)) +
      geom_joy_gradient(scale = 5, rel_min_height = 0.01, alpha = 0.5, col = "white") +
      scale_fill_viridis(option = viridiscolor) +
      theme_light() +
      ggtitle("nano_joy_length") +
      ylab("Hours since start") +
      xlab("Read length") +
      guides(fill = F))
  }
  if(transformation == "log10"){
    return(fread(input = inputfile) %>%
             as.data.table() %>%
             mutate(timeindex = round((start_time + 3600)/3600)) %>%
             group_by(timeindex) %>%
             ggplot(aes(x = sequence_length_template, y = factor(timeindex), fill = ..x..)) +
             geom_joy_gradient(scale = 5, rel_min_height = 0.01, alpha = 0.5, col = "white") +
             scale_fill_viridis(option = viridiscolor) +
             theme_light() +
             ylab("Hours since start") +
             ggtitle("nano_joy_length") +
             xlab("Read length") +
             guides(fill = F) +
             scale_x_log10())
  }
}

#' joy distribution of read quality
#' Takes in the `sequencing_summary` file of albacore as input and visualizes the readquality distribution as a histogram using ggplot2
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @return joy distribution of read quality
#' @export
#' @import ggplot2 dplyr viridis data.table ggjoy
nano_joy_quality <- function(inputfile,  viridiscolor = "viridis"){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
    mutate(timeindex = round((start_time + 3600)/3600)) %>%
    group_by(timeindex) %>%
    ggplot(aes(x = mean_qscore_template, y = factor(timeindex), fill = ..x..)) +
      geom_joy_gradient(scale = 5, rel_min_height = 0.01, alpha = 0.5, col = "white") +
      scale_fill_viridis(option = viridiscolor) +
      theme_light() +
      ggtitle("nano_joy_quality") +
      ylab("Hours since start") +
      xlab("Read quality") +
      guides(fill = F))
}

#' Summary of n longest reads and their basecall quality
#' Takes in the `sequencing_summary` file of albacore as input and give you a summary of longest reads and their quality
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param viridiscolor Viridis color palette option
#' @return summary of n longest reads
#' @export
#' @import dplyr data.table
nano_readsum <- function(inputfile, number = 5){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
      dplyr::arrange(desc(mean_qscore_template)) %>%
      head(number) %>%
      mutate(position = 1:number,
             readlength = sequence_length_template,
             basecall_quality = mean_qscore_template) %>%
      dplyr::select(position, basecall_quality,readlength)
  )
}



#' Summary of basic nanopore run statistics
#' Takes in the `sequencing_summary` file of albacore as input and give you a basic summary
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param quality Mean quality filter option for every read
#' @return basic nanopore run statistics
#' @export
#' @import dplyr data.table
nano_stats <- function(inputfile, quality = 0){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
           dplyr::filter(mean_qscore_template >= quality) %>%
           mutate("number of reads"    = length(filename),
                  "total bases"        = sum(sequence_length_template),
                  "median read length" = median(sequence_length_template),
                  "mean read length"   = mean(sequence_length_template)) %>%
           dplyr::select("number of reads", "total bases", "median read length", "mean read length") %>%
           head(1) %>%
           t()
  )
}


#' N50 of nanopore run
#' Takes in the `sequencing_summary` file of albacore as input and returns the N50
#' @param inputfile Filepath to sequencing_summary.txt output of albacore
#' @param quality Mean quality filter option for every read
#' @return N50
#' @export
#' @import dplyr data.table
nano_n50 <- function(inputfile, quality = 0){
  return(fread(input = inputfile) %>%
           as.data.table() %>%
           dplyr::filter(mean_qscore_template >= quality) %>%
           dplyr::arrange(desc(sequence_length_template)) %>%
           mutate(cumulative = cumsum(sequence_length_template)) %>%
           dplyr::select(sequence_length_template, cumulative) %>%
           mutate(N50 = sequence_length_template[which.min(abs(cumulative - max(cumulative)/2))]) %>%
           dplyr::select(N50) %>%
           head(1)
  )
}



