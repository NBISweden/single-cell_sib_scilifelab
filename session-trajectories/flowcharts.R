library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

graph <- grViz("digraph flowchart {
  # node definitions with substituted label text
  node [fontname = Helvetica, shape = rectangle]        
  reads [label = 'Reads \n fastq file']
  genome [label = 'Genome or transcriptome \n fasta file']
  aligner [label = 'Aligner \n STAR', style=filled]
  alignment [label = 'Alignment \n bam file']
  annotation [label = 'Genome annotation \n GTF file']
  quantification [label = 'Quantification & UMI collapse \n velocyto', style=filled]
  spliced [label = 'Spliced counts \n matrix cells vs features']
  unspliced [label = 'Unspliced counts \n matrix cells vs features']
  analysis [label = 'Analysis \n Seurat, velocyto.R, dyno', style = filled]

  # edge definitions with the node IDs
  reads -> aligner;
  genome -> aligner;
  aligner -> alignment;
  alignment -> quantification;
  annotation -> quantification;
  quantification -> spliced;
  quantification -> unspliced;
  spliced -> analysis;
  unspliced -> analysis;
  }
")
graph





graph <- grViz("digraph flowchart {
  # node definitions with substituted label text
  node [fontname = Helvetica, shape = rectangle]        
  reads [label = 'Reads \n fastq file']
  genome [label = 'Genome or transcriptome \n fasta file']
  aligner [label = 'Aligner \n STAR', style=filled]
  alignment [label = 'Alignment \n bam file']
  annotation [label = 'Genome annotation \n GTF file']
  quantification [label = 'Quantification & UMI collapse \n velocyto', style=filled]
  spliced [label = 'Spliced counts \n matrix cells vs features']
  unspliced [label = 'Unspliced counts \n matrix cells vs features']
  analysis [label = 'Analysis \n Seurat, velocyto.R, dyno', style = filled]

  # edge definitions with the node IDs
  reads -> aligner;
  genome -> aligner;
  aligner -> alignment;
  alignment -> quantification;
  annotation -> quantification;
  quantification -> spliced;
  quantification -> unspliced;
  spliced -> analysis;
  unspliced -> analysis;
  }
")
graph
