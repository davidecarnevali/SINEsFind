# SINEsFind

### Abstract

Short Interspersed Element (SINE) retrotransposons are one of the 
most abundant DNA repeat elements in the human genome. They have been found to 
impact the expression of protein-coding genes, but the possible roles in cell 
physiology of their noncoding RNAs, generated by RNA polymerase (Pol) III, are 
just starting to be elucidated. For this reason, Short Interspersed Element 
(SINE) expression profiling is becoming mandatory to obtain a comprehensive 
picture of their regulatory roles. However, their repeated nature and frequent 
location within Pol II-transcribed genes represent a serious obstacle to the 
identification and quantification of genuine, Pol III-derived SINE transcripts 
at single-locus resolution on a genomic scale. Among the recent Next Generation 
Sequencing technologies, only RNA sequencing (RNA-Seq) holds the potential to 
solve these issues, even though both technical and biological matters need to 
be taken into account. A bioinformatic pipeline has been recently set up that, 
by exploiting RNA-seq features and knowledge of SINE transcription mechanisms, 
allows for easy identification and profiling of transcriptionally active genomic

SINEsFind allows to detect free SINE RNAs by using paired-end RNA-Seq data.
It works by comparing the expression coverage upstream and downstream the 
annotad human SINE element (Alus and Mammalian-wide Interpsersed Repeats)

#### Requirements:

 * HTSeq
 * EMBOSS
 * pyBedTools
 * pyBigWig
