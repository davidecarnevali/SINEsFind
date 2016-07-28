import HTSeq
import pyBigWig
import argparse
import csv
import re
import subprocess
import sys
import time
import warnings
from Bio import AlignIO
from Bio.Emboss.Applications import NeedleCommandline
from pybedtools import BedTool

parser = argparse.ArgumentParser(
    description = 'This script takes a coverage file in BAM or BEDGRAPH format and a SINE annotation file in GTF\
    format to find genuine SINE transcripts. Version 2.1.6 July 2016',
    epilog = 'Written by Davide Carnevali davide.carnevali@nemo.unipr.it')
parser.add_argument("-s", "--stranded",
                    help="Use this option if using a stranded coverage file(s). If using bam file make sure it is\
                    generated with TopHat as this program use the 'XS' tag to identify the strand of the transcript\
                    from which the reads come from",
                    action="store_true")
parser.add_argument("-t", "--filetype", choices=['bam', 'bw'],
                    help="specify coverage file type: default 'bam'.  BigWig stranded files should be comma separated, with plus signal preceding the minus one",
                    default='bam')
parser.add_argument("-bg", "--background", type=int,
                    help="Set how many times the SINE area coverage should be greater than background . Default: 1",
                    default='1')
parser.add_argument("-r", "--ratio", type=float, help="The ratio  left/central or out/central coverage area")
parser.add_argument("-TSS", "--tss", type=int, help="Set the width of Transcription Start Site. Default: 50",
                    default='50')
parser.add_argument("-LR", "--left_region", type=int, help="Set the region size in nt. Default: 100", default='100')
parser.add_argument("-RR", "--right_region", type=int, help="Set the region size in nt. Default: 200", default='200')
parser.add_argument("-OR", "--out_region", type=int, help="Set the region size in nt. Default: 100", default='100')
parser.add_argument("coverage",
                    help="Coverage file to be processed, either in BAM or BigWig format. Using BigWig files the script run much faster (x10). If using BigWig make sure the coverage is made up only of uniquely mapped reads")
parser.add_argument("gtf", help="annotation file in GFF/GTF format")
parser.add_argument("genome", help="reference genome in fasta format")
parser.add_argument("output", help="output filename")
args = parser.parse_args()

genome = BedTool(args.genome)
MIR = 'ACAGTATAGCATAGTGGTTAAGAGCACGGACTCTGGAGCCAGACTGCCTGGGTTCGAATCCCGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTTACTTAACCTCTCTGTGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATAGTACCTACCTCATAGGGTTGTTGTGAGGATTAAATGAGTTAATACATGTAAAGCGCTTAGAACAGTGCCTGGCACATAGTAAGCGCTCAATAAATGTTGGTTATTA'
MIRb = 'CAGAGGGGCAGCGTGGTGCAGTGGAAAGAGCACGGGCTTTGGAGTCAGGCAGACCTGGGTTCGAATCCTGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTCACTTAACCTCTCTGAGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATACCTACCTCGCAGGGTTGTTGTGAGGATTAAATGAGATAATGCATGTAAAGCGCTTAGCACAGTGCCTGGCACACAGTAAGCGCTCAATAAATGGTAGCTCTATTATT'
MIRc = 'CGAGGCAGTGTGGTGCAGTGGAAAGAGCACTGGACTTGGAGTCAGGAAGACCTGGGTTCGAGTCCTGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTCACTTAACCTCTCTGAGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATACCTGCCCTGCCTACCTCACAGGGTTGTTGTGAGGATCAAATGAGATAATGTATGTGAAAGCGCTTTGTAAACTGTAAAGTGCTATACAAATGTAAGGGGTTATTATTATT'
MIR3 = 'TTCTGGAAGCAGTATGGTATAGTGGAAAGAACAACTGGACTAGGAGTCAGGAGACCTGGGTTCTAGTCCTAGCTCTGCCACTAACTAGCTGTGTGACCTTGGGCAAGTCACTTCACCTCTCTGGGCCTCAGTTTTCCTCATCTGTAAAATGAGNGGGTTGGACTAGATGATCTCTAAGGTCCCTTCCAGCTCTAACATTCTATGATTCTATGATTCTAAAAAAA'
ALU = 'GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATTGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTCTCA'

start_time = time.time()
annotation = HTSeq.GFF_Reader(args.gtf)
cvg_plus = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
cvg_minus = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
cvg = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
left = 0
central = 0
right = 0
count_plus = 0
count_minus = 0
count = 0
alu_list = []
char = re.compile('-*')
char2 = re.compile('[-ATGC]*')


# Build the coverage vectors for + and - strand based on XS tag, using uniquely mapped reads

def cvg_bam(file):
    global count_plus
    global count_minus
    for read in file:
        if read.aligned and read.optional_field('NH') == 1:
            if read.optional_field('XS') == '+':
                for cigopt in read.cigar:
                    if cigopt.type == 'M':
                        cvg_plus[cigopt.ref_iv] += 1
                        count_plus += cigopt.ref_iv.length
            elif read.optional_field('XS') == '-':
                for cigopt in read.cigar:
                    if cigopt.type == 'M':
                        cvg_minus[cigopt.ref_iv] += 1
                        count_minus += cigopt.ref_iv.length
    return cvg_plus, cvg_minus, count_plus, count_minus


def cvg_bam_unstranded(file):
    global count
    for read in file:
        if read.aligned and read.optional_field('NH') == 1:
            for cigopt in read.cigar:
                if cigopt.type == 'M':
                    cvg[cigopt.ref_iv] += 1
                    count += cigopt.ref_iv.length
    return cvg, count


# Calculate coverage in the Left/Center/Right arm of the SINEs

def frf_stranded(gtf, peak):
    for element in gtf:
        if element.iv.strand == '+':
            if sum(bw_plus.values(element.iv.chrom, element.iv.start,element.iv.end)) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = int(sum(bw_plus.values(element.iv.chrom, (element.iv.start - aln_start),
                                                                      (element.iv.start - aln_start + aln_end))))
                    right = int(sum(bw_plus.values(element.iv.chrom,(element.iv.start - aln_start + aln_end), (
                                                                        element.iv.start - aln_start + (
                                                                            aln_end + args.right_region)))))
                    left = int(sum(bw_plus.values(element.iv.chrom, (element.iv.start - (aln_start + args.tss + args.left_region)),
                                                                   (element.iv.start - aln_start - args.tss))))
                    out = int(sum(bw_plus.values(element.iv.chrom, (element.iv.start - aln_start + (aln_end + args.right_region)), (element.iv.start - aln_start + (
                        aln_end + args.right_region + args.out_region)))))

                else:
                    continue
            else:
                continue

        elif element.iv.strand == '-':
            if sum(bw_minus.values(element.iv.chrom, element.iv.start, element.iv.end)) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = int(sum(bw_minus.values(element.iv.chrom,(element.iv.end + aln_start - aln_end),
                                                                       (element.iv.end + aln_start))))
                    right = int(sum(bw_minus.values(element.iv.chrom, (element.iv.end + aln_start - (aln_end + args.right_region)),
                                                                     (element.iv.end + aln_start - aln_end))))
                    left = int(sum(bw_minus.values(element.iv.chrom, (element.iv.end + aln_start + args.tss),
                                                              (element.iv.end + aln_start + args.tss + args.left_region))))
                    out = int(sum(bw_minus.values(element.iv.chrom, (element.iv.end + aln_start - (aln_end + args.right_region + args.out_region)), (
                                                                       element.iv.end + aln_start - (aln_end + args.out_region)))))

                else:
                    continue
            else:
                continue
        if args.ratio:
            if left < (central / aln_end) * args.left_region * args.ratio and right < (
                        central / aln_end) * args.right_region and out < (central / aln_end) * args.out_region * args.ratio:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])
        else:
            if left == 0 and right < (central / aln_end) * args.right_region and out == 0:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])


def frf_unstranded(gtf, peak):
    for element in gtf:
        if element.iv.strand == '+':
            if sum(bw.values(element.iv.chrom, element.iv.start,
                                                  element.iv.end)) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = int(sum(bw.values(element.iv.chrom, (element.iv.start - aln_start),
                                                                 (element.iv.start - aln_start + aln_end))))
                    right = int(sum(bw.values(element.iv.chrom,(element.iv.start - aln_start + aln_end),(
                                                                   element.iv.start - aln_start + (
                                                                       aln_end + args.right_region)))))
                    left = int(sum(bw.values(element.iv.chrom, (
                        element.iv.start - (aln_start + args.tss + args.left_region)),
                                                              (element.iv.start - aln_start - args.tss))))
                    out = int(sum(bw.values(element.iv.chrom, (
                        element.iv.start - aln_start + (aln_end + args.right_region)), (element.iv.start - aln_start + (
                        aln_end + args.right_region + args.out_region)))))

                else:
                    continue
            else:
                continue

        elif element.iv.strand == '-':
            if sum(bw.values(element.iv.chrom, element.iv.start,
                                                  element.iv.end)) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = int(sum(bw.values(element.iv.chrom,
                                                                 (element.iv.end + aln_start - aln_end),
                                                                 (element.iv.end + aln_start))))
                    right = int(sum(bw.values(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region)),
                                                               (element.iv.end + aln_start - aln_end))))
                    left = int(sum(bw.values(element.iv.chrom, (element.iv.end + aln_start + args.tss),
                                                              (
                                                                  element.iv.end + aln_start + args.tss + args.left_region))))
                    out = int(sum(bw.values(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region + args.out_region)), (
                                                                 element.iv.end + aln_start - (
                                                                     aln_end + args.out_region)))))

                else:
                    continue
            else:
                continue

        if args.ratio:
            if left < (central / aln_end) * args.left_region * args.ratio and right < (
                        central / aln_end) * args.right_region and out < (central / aln_end) * args.out_region * args.ratio:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])
        else:
            if left == 0 and right < (central / aln_end) * args.right_region and out == 0:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])


def frf_stranded_bam(gtf, peak):
    for element in gtf:
        if element.iv.strand == '+':
            if sum(list(cvg_plus[HTSeq.GenomicInterval(element.iv.chrom, element.iv.start,
                                                       element.iv.end)])) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = sum(list(cvg_plus[HTSeq.GenomicInterval(element.iv.chrom, (element.iv.start - aln_start),
                                                                      (element.iv.start - aln_start + aln_end))]))
                    right = sum(list(cvg_plus[HTSeq.GenomicInterval(element.iv.chrom,
                                                                    (element.iv.start - aln_start + aln_end), (
                                                                        element.iv.start - aln_start + (
                                                                            aln_end + args.right_region)))]))
                    left = sum(list(cvg_plus[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.start - (aln_start + args.tss + args.left_region)),
                                                                   (element.iv.start - aln_start - args.tss))]))
                    out = sum(list(cvg_plus[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.start - aln_start + (aln_end + args.right_region)), (element.iv.start - aln_start + (
                        aln_end + args.right_region + args.out_region)))]))

                else:
                    continue
            else:
                continue

        elif element.iv.strand == '-':
            if sum(list(cvg_minus[HTSeq.GenomicInterval(element.iv.chrom, element.iv.start,
                                                        element.iv.end)])) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = sum(list(cvg_minus[HTSeq.GenomicInterval(element.iv.chrom,
                                                                       (element.iv.end + aln_start - aln_end),
                                                                       (element.iv.end + aln_start))]))
                    right = sum(list(cvg_minus[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region)),
                                                                     (element.iv.end + aln_start - aln_end))]))
                    left = sum(list(cvg_minus[
                                        HTSeq.GenomicInterval(element.iv.chrom, (element.iv.end + aln_start + args.tss),
                                                              (
                                                                  element.iv.end + aln_start + args.tss + args.left_region))]))
                    out = sum(list(cvg_minus[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region + args.out_region)), (
                                                                       element.iv.end + aln_start - (
                                                                           aln_end + args.out_region)))]))

                else:
                    continue
            else:
                continue
        if args.ratio:
            if left < (central / aln_end) * args.left_region * args.ratio and right < (
                        central / aln_end) * args.right_region and out < (central / aln_end) * args.out_region * args.ratio:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])
        else:
            if left == 0 and right < (central / aln_end) * args.right_region and out == 0:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])


def frf_unstranded_bam(gtf, peak):
    for element in gtf:
        if element.iv.strand == '+':
            if sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, element.iv.start,
                                                  element.iv.end)])) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (element.iv.start - aln_start),
                                                                 (element.iv.start - aln_start + aln_end))]))
                    right = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom,
                                                               (element.iv.start - aln_start + aln_end), (
                                                                   element.iv.start - aln_start + (
                                                                       aln_end + args.right_region)))]))
                    left = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.start - (aln_start + args.tss + args.left_region)),
                                                              (element.iv.start - aln_start - args.tss))]))
                    out = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.start - aln_start + (aln_end + args.right_region)), (element.iv.start - aln_start + (
                        aln_end + args.right_region + args.out_region)))]))

                else:
                    continue
            else:
                continue

        elif element.iv.strand == '-':
            if sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, element.iv.start,
                                                  element.iv.end)])) > args.background * (
                        (element.iv.end - element.iv.start) * peak):
                if "MIR" in element.attr['gene_id'] or "Alu" in element.attr['gene_id']:
                    aln_start, aln_end = needle(element.iv.chrom, element.iv.start, element.iv.end,
                                                element.attr['gene_id'], element.score, element.iv.strand)
                    central = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom,
                                                                 (element.iv.end + aln_start - aln_end),
                                                                 (element.iv.end + aln_start))]))
                    right = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region)),
                                                               (element.iv.end + aln_start - aln_end))]))
                    left = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (element.iv.end + aln_start + args.tss),
                                                              (
                                                                  element.iv.end + aln_start + args.tss + args.left_region))]))
                    out = sum(list(cvg[HTSeq.GenomicInterval(element.iv.chrom, (
                        element.iv.end + aln_start - (aln_end + args.right_region + args.out_region)), (
                                                                 element.iv.end + aln_start - (
                                                                     aln_end + args.out_region)))]))

                else:
                    continue
            else:
                continue

        if args.ratio:
            if left < (central / aln_end) * args.left_region * args.ratio and right < (
                        central / aln_end) * args.right_region and out < (central / aln_end) * args.out_region * args.ratio:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])
        else:
            if left == 0 and right < (central / aln_end) * args.right_region and out == 0:
                alu_list.append(
                    [peak, aln_start, aln_end, element.attr['transcript_id'], element.iv.chrom, element.iv.start,
                     element.iv.end, element.iv.strand, left, central, right, out])

# Perform global alignment, with Needle algorithm, of the element to its consensus sequence to define the start/end of the central region

def needle(chrom, start, end, name, score, strand):
    aln_start = 0
    aln_end = 0
    item = BedTool([(chrom, start, end, name, score, strand)])
    item = item.sequence(fi=genome, s=True)
    temp = open(item.seqfn).read().split('\n')[1].upper()
    if name == "MIRb":
        sine_length = 269
        needle_cline = NeedleCommandline(asequence="asis:" + MIRb, bsequence="asis:" + temp, gapopen=10, gapextend=0.5,
                                         outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 shell=(sys.platform != "win32"))
        child.wait()
        align = AlignIO.read(child.stdout, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif name == "MIRc":
        sine_length = 269
        needle_cline = NeedleCommandline(asequence="asis:" + MIRc, bsequence="asis:" + temp, gapopen=10, gapextend=0.5,
                                         outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 shell=(sys.platform != "win32"))
        child.wait()
        align = AlignIO.read(child.stdout, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif name == "MIR3":
        sine_length = 225
        needle_cline = NeedleCommandline(asequence="asis:" + MIR3, bsequence="asis:" + temp, gapopen=10, gapextend=0.5,
                                         outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 shell=(sys.platform != "win32"))
        child.wait()
        align = AlignIO.read(child.stdout, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif name == "MIR":
        sine_length = 261
        needle_cline = NeedleCommandline(asequence="asis:" + MIR, bsequence="asis:" + temp, gapopen=10, gapextend=0.5,
                                         outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 shell=(sys.platform != "win32"))
        child.wait()
        align = AlignIO.read(child.stdout, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    elif "Alu" in name:
        sine_length = 284
        needle_cline = NeedleCommandline(asequence="asis:" + ALU, bsequence="asis:" + temp, gapopen=10, gapextend=0.5,
                                         outfile='stdout')
        child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 shell=(sys.platform != "win32"))
        child.wait()
        align = AlignIO.read(child.stdout, "emboss")
        aln_start = char.search(str(align[1, :].seq)).end()
        aln_end = char2.search(str(align[1, :].seq)).end()

    return aln_start, aln_end


# Check srguments and call specific functions

if args.filetype == 'bam':
    bamfile = HTSeq.BAM_Reader(args.coverage)
    if args.stranded:
        print "Start reading bam file, this will take a while....."
        cvg_bam(bamfile)
        peak = int(round((count_plus + count_minus) / (6 * 10 ** 9), 0))
        cvg_plus.write_bedgraph_file(args.output + "_plus.bg")
        cvg_minus.write_bedgraph_file(args.output + "_minus.bg")
        print "Start applying Flanking Region Filter"
        frf_stranded_bam(annotation, peak)

    else:
        print "Start reading bam file, this will take a while....."
        cvg_bam_unstranded(bamfile)
        peak = int(round(count / (6 * 10 ** 9), 0))
        cvg.write_bedgraph_file(args.output + ".bg")
        print "Start applying Flanking Region Filter"
        frf_unstranded_bam(annotation, peak)

elif args.filetype == 'bw':
    if args.stranded:
        bw_plus = pyBigWig.open(args.coverage.strip().split(",")[0])
        bw_minus = pyBigWig.open(args.coverage.strip().split(",")[1])
        peak = int(round((bw_plus.header()["sumData"] + bw_minus.header()["sumData"])/ (6 * 10 ** 9)))
        print "Start applying Flanking Region Filter"
        frf_stranded(annotation, peak)
    else:
        bw = pyBigWig.open(args.coverage)
        peak = int(round(bw.header()["sumData"] / (6  * 10 ** 9)))
        print "Start applying Flanking Region Filter"
        frf_unstranded(annotation, peak)

# Write to output the Bona Fide SINEs

with open(args.output, 'wb') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(alu_list)

print "Finished!"
print "Time elapsed %s" % (time.time() - start_time)