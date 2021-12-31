import HTSeq, sys, time, argparse, warnings, csv, re, subprocess
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from pybedtools import BedTool
import pandas as pd
import os

parser = argparse.ArgumentParser(
    description = 'This script takes in an annotation file of human Alu and MIR elements and the corresponding human reference genome sequence\
    to produce the index for SINEs_find.py. The index file is the annotation input file with added columns reporting the lenght, start and end of the alignment\
    of the element sequence to its consensus sequence on the reference genome',
    epilog = 'Written by Davide Carnevali davide.carnevali@unipr.it')
parser.add_argument("annotation", help="Annotation file in GTF format. Should be the same version of the human reference genome file")
parser.add_argument("genome", help="Human reference genome sequence. Should be the same version of the annotation file")
parser.add_argument("output", help="output filename")
args = parser.parse_args()

genome = BedTool(args.genome)
MIR = 'ACAGTATAGCATAGTGGTTAAGAGCACGGACTCTGGAGCCAGACTGCCTGGGTTCGAATCCCGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTTACTTAACCTCTCTGTGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATAGTACCTACCTCATAGGGTTGTTGTGAGGATTAAATGAGTTAATACATGTAAAGCGCTTAGAACAGTGCCTGGCACATAGTAAGCGCTCAATAAATGTTGGTTATTA'
MIRb = 'CAGAGGGGCAGCGTGGTGCAGTGGAAAGAGCACGGGCTTTGGAGTCAGGCAGACCTGGGTTCGAATCCTGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTCACTTAACCTCTCTGAGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATACCTACCTCGCAGGGTTGTTGTGAGGATTAAATGAGATAATGCATGTAAAGCGCTTAGCACAGTGCCTGGCACACAGTAAGCGCTCAATAAATGGTAGCTCTATTATT'
MIRc = 'CGAGGCAGTGTGGTGCAGTGGAAAGAGCACTGGACTTGGAGTCAGGAAGACCTGGGTTCGAGTCCTGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTCACTTAACCTCTCTGAGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATACCTGCCCTGCCTACCTCACAGGGTTGTTGTGAGGATCAAATGAGATAATGTATGTGAAAGCGCTTTGTAAACTGTAAAGTGCTATACAAATGTAAGGGGTTATTATTATT'
MIR3 = 'TTCTGGAAGCAGTATGGTATAGTGGAAAGAACAACTGGACTAGGAGTCAGGAGACCTGGGTTCTAGTCCTAGCTCTGCCACTAACTAGCTGTGTGACCTTGGGCAAGTCACTTCACCTCTCTGGGCCTCAGTTTTCCTCATCTGTAAAATGAGNGGGTTGGACTAGATGATCTCTAAGGTCCCTTCCAGCTCTAACATTCTATGATTCTATGATTCTAAAAAAA'
MIR1_Amn = 'TAGGGAGGCAGTGTGGTCTAGTGGATAGAGCACTGGACTGGGACTCGGGAGACCTGGGTTCNANTCCCGGCTCTGCCACTNGCCNGCTGNGTGACCTTGGGCAAGTCACTTNACCTCTCTGNGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATGATACTGACCTCCTTTGTAAAGTGCTTTGAGATCTACTGATGAAAAGTGCTACATAAGAGCTAGGTATTATTAT'
ALU = 'GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATTGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
FAM = 'GCCGGGCGCGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGTAGTGCGCTATGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAA'
FLAM = 'GCCGGGCGCGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTAAAAAAAA'
FRAM = 'GCCGGGCGCGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCCGTCTCAAAAAAAAAA'

start_time = time.time()
annotation = HTSeq.GFF_Reader(args.annotation)
alu_list = []
char = re.compile('-*')
char2 = re.compile('[-NATGCatgc]*')
    
def needle(chrom, start, end, name, score, strand):
    tmpfile = open("tmpaln.txt", 'w')
    n = 0
    item=BedTool([(chrom, start, end, name, score, strand)])
    item = item.sequence(fi=genome, s=True)
    temp = open(item.seqfn).read().split('\n')[1]
    if name == "MIRb":
        sine_length = 269
        needle_cline = NeedleCommandline(asequence="asis:"+MIRb, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()
                    
    elif name == "MIRc":
        sine_length = 269
        needle_cline = NeedleCommandline(asequence="asis:"+MIRc, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()
                    
    elif name == "MIR3":
        sine_length = 225
        needle_cline = NeedleCommandline(asequence="asis:"+MIR3, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif name == "MIR1_Amn":
        sine_length = 231
        needle_cline = NeedleCommandline(asequence="asis:"+MIR3, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif name == "MIR":
        sine_length = 261
        needle_cline = NeedleCommandline(asequence="asis:"+MIR, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()
        
    elif "Alu" in name:
        sine_length = 313
        needle_cline = NeedleCommandline(asequence="asis:"+ALU, bsequence="asis:"+temp,gapopen=10, gapextend=0.5, outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=tmpfile, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
        child.wait()
        align = AlignIO.read(tmpfile.name, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()
    
    tmpfile.close()
        
    return (aln_start, aln_end, sine_length)

# To mantain GTF notation, which is 1 based, we add +1 to element.iv.start
for element in annotation:
    aln_start, aln_end, sine_length = needle(element.iv.chrom, element.iv.start, element.iv.end, element.attr['gene_id'], int(element.score), element.iv.strand)
    if element.iv.strand == "+":
        alu_list.append([element.iv.chrom, element.source, "exon", element.iv.start + 1, element.iv.end, "0", element.iv.strand, ".", "gene_id " + "\""+element.attr['transcript_id']+"\"; " + "transcript_id " + "\""+element.attr['transcript_id']+"\";",element.iv.start - aln_start, (element.iv.start - aln_start) + aln_end, aln_start, aln_end])
    else:
        alu_list.append([element.iv.chrom, element.source, "exon", element.iv.start + 1, element.iv.end, "0", element.iv.strand, ".", "gene_id " + "\""+element.attr['transcript_id']+"\"; " + "transcript_id " + "\""+element.attr['transcript_id']+"\";",(element.iv.end + aln_start) - aln_end, element.iv.end + aln_start, aln_start, aln_end])

final_list = pd.DataFrame(alu_list)
final_list.to_csv(args.output,sep="\t", header=False,index=False,doublequote=False,quotechar='\'',escapechar='')

print("Finished!")
print("Time elapsed {}".format(time.time() - start_time))