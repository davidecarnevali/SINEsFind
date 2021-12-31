from __future__ import division
import HTSeq
import pyBigWig
import numpy as np
import argparse
import csv
import time

#from memory_profiler import profile
import math

parser = argparse.ArgumentParser(
    description = 'This script takes in a RNA-Seq coverage file in BAM or \
        BIGWIG format and a SINE annotation file in GTF extended\
            format to find genuine Pol III-derived SINE transcripts.\
                Version 2.4 December 2021',
    epilog = 'Written by Davide Carnevali davide.carnevali@crg.eu')
parser.add_argument("-g", "--ref_genome", choices=['GRCh38', 'GRCh37'],
                    help = "Set the reference genome you are working on.\
                        Default is GRCh38", default='GRCh38')
parser.add_argument("-s", "--stranded", choices=['auto', 'no', 'yes', 'reverse'],
                    help="Use this option if using a stranded coverage file(s).\
                        For first-strand synthesis use 'reverse' while for\
                            second-strand synthesis use 'yes'", default = 'auto')
parser.add_argument("-t", "--filetype", choices=['bam', 'bw'],
                    help="specify coverage file type: default 'bam'. Bamfile\
                        should be sorted by coordinate. BigWig stranded files\
                            should be comma separated, with plus signal\
                                preceding the minus one",
                    default='bam')
parser.add_argument("-bg", "--background", type=int,
                    help="Set how many times the SINE coverage area should be\
                        greater than background . Default: 1",
                    default=1)
parser.add_argument("-lr", "--lratio", type=float, help="The ratio  left/central\
                    coverage area. Default: 0.1", default=0.1)
parser.add_argument("-or", "--oratio", type=float, help="The ratio  out/central\
                    coverage area. Default 0.1", default=0.1)
parser.add_argument("-TSS", "--tss", type=int, help="Set the width of \
                    Transcription Start Site. Default: 15",
                    default=15)
parser.add_argument("-LR", "--left_region", type=int, help="Set the region \
                    size in nt. Default: 100", default=100)
parser.add_argument("-RR", "--right_region", type=int, help="Set the region \
                    size in nt. Default: 300", default=300)
parser.add_argument("-OR", "--out_region", type=int, help="Set the region size\
                    in nt. Default: 100", default=100)
parser.add_argument("-f", "--fraction", type=float, help="Fraction of the \
                    SINE central part that should have an expression coverage\
                        '-bg' times over the background. Default 50 per cent",
                        default=0.5)
parser.add_argument("coverage",
                    help="Coverage file to be processed, either in BAM or\
                        BigWig format. Using BigWig files the script run much\
                            faster (x10). If using BigWig make sure the coverage\
                                is made up only of uniquely mapped reads")
parser.add_argument("gtf", help="The extended annotation file generate using\
                    AnnoGenerate.py")
parser.add_argument("chroms", help="chromosomes length used for mapping\
                    (as in bam header)")
parser.add_argument("output", help="output filename")
args = parser.parse_args()


#memLimit= 8 * 1024 * 1024 * 1024
#resource.setrlimit(resource.RLIMIT_AS, (memLimit, memLimit))
start_time = time.time()
count_sine = 0
alu_list = []
logfile = []
if args.ref_genome == "GRCh38":
    test_region = "chr12:6537505-6538052"
elif args.ref_genome == "GRCh37":
    test_region = "chr12:6646671-6647218"
annotation = [line.strip().split("\t") for line in open(args.gtf, "r")]
for i in annotation:
    i[8] = i[8].strip().split("\"")[3]
chrom_list ={}
with open(args.chroms) as f:
    for line in f:
        (key, val) = line.split()
        if "_" not in key and "chr" in key and "EBV" not in key:
            chrom_list[key] = int(val)
chromosomes = list(chrom_list.items())
chromosomes.sort()
if args.filetype == 'bw':
    if args.stranded == 'reverse' or args.stranded == 'yes':
        bw_plus = pyBigWig.open(args.coverage.strip().split(",")[0])
        bw_minus = pyBigWig.open(args.coverage.strip().split(",")[1])
    elif args.stranded == 'no':
        bw = pyBigWig.open(args.coverage)
# Build the coverage vectors for + and - strand based on XS tag, using uniquely mapped reads
def strdchk(file):
    firststrand = 0
    secondstrand = 0
    total = 0
    for read in file.fetch( region = test_region ):
        if read.proper_pair and read.aligned and not read.failed_platform_qc and read.optional_field('NH') == 1:
            if (read.pe_which == 'second' and read.iv.strand == '+') or (read.pe_which == 'first' and read.iv.strand == '-'):
                firststrand += 1
            elif (read.pe_which == 'second' and read.iv.strand == '-') or (read.pe_which == 'first' and read.iv.strand == '+'):
                secondstrand += 1
            total += 1
    if firststrand/total >= 0.8:
        args.stranded = 'reverse'
    elif secondstrand/total >= 0.8:
        args.stranded = 'yes'
    else:
        args.stranded = 'no'

def cvg_bam(file):
    bedgraphFile_plus = open(args.output+"_plus.bg", 'w')
    bedgraphFile_minus = open(args.output+"_minus.bg", 'w')
    global cvg_plus
    global cvg_minus
    global chromosomes
    run_time = time.time()
    for interval in chromosomes:
        count_plus = 0
        count_minus = 0
        p = HTSeq.GenomicPosition( interval[0], interval[1], "." )
        window = HTSeq.GenomicInterval(p.chrom, 1, p.pos, ".")
        cvg_plus = HTSeq.GenomicArray({p.chrom: p.pos}, stranded=False, typecode="i")
        cvg_minus = HTSeq.GenomicArray({p.chrom: p.pos}, stranded=False, typecode="i")
        for read in file[window]:
            if read.proper_pair and read.aligned and not read.failed_platform_qc and read.optional_field('NH') == 1:
                cig_list = []
                for x in read.cigar:
                    cig_list.append(x.type)
                if 'N' not in cig_list:
                    if args.stranded == 'reverse':
                        if (read.pe_which == 'second' and read.iv.strand == '+') or (read.pe_which == 'first' and read.iv.strand == '-'):
                            for cigopt in read.cigar:
                                if cigopt.type == 'M':
                                    cvg_plus[cigopt.ref_iv] += 1
                                    count_plus += cigopt.ref_iv.length
                        elif (read.pe_which == 'second' and read.iv.strand == '-') or (read.pe_which == 'first' and read.iv.strand == '+'):
                            for cigopt in read.cigar:
                                if cigopt.type == 'M':
                                    cvg_minus[cigopt.ref_iv] += 1
                                    count_minus += cigopt.ref_iv.length
                    elif args.stranded == 'yes':
                        if (read.pe_which == 'first' and read.iv.strand == '+') or (read.pe_which == 'second' and read.iv.strand == '-'):
                            for cigopt in read.cigar:
                                if cigopt.type == 'M':
                                    cvg_plus[cigopt.ref_iv] += 1
                                    count_plus += cigopt.ref_iv.length
                        elif (read.pe_which == 'first' and read.iv.strand == '-') or (read.pe_which == 'second' and read.iv.strand == '+'):
                            for cigopt in read.cigar:
                                if cigopt.type == 'M':
                                    cvg_minus[cigopt.ref_iv] += 1
                                    count_minus += cigopt.ref_iv.length
        cvg_plus.write_bedgraph_file(bedgraphFile_plus)
        cvg_minus.write_bedgraph_file(bedgraphFile_minus)
        peak_plus = int(math.ceil(count_plus / p.pos))
        peak_minus = int(math.ceil(count_minus / p.pos))
        print("Time elapsed for coverage calculation of chromosome {}: {} seconds".format(p.chrom, time.time() - run_time))
        cvg_time = time.time()
        print("Start applying Flanking Region Filter for Alus on chromosome ".format(p.chrom))
        frf_stranded_bam(annotation, peak_plus, peak_minus, p.chrom)
        print("Time elapsed for Pol III Alus discovery on chromosome {}: {} seconds".format(p.chrom, time.time() - cvg_time))
        run_time = time.time()


def cvg_bam_unstranded(file):
    bedgraphFile = open(args.output+".bg", 'w')
    global cvg
    global chromosomes
    run_time = time.time()
    for interval in chromosomes:
        count = 0
        p = HTSeq.GenomicPosition( interval[0], interval[1], "." )
        window = HTSeq.GenomicInterval(p.chrom, 1, p.pos, ".")
        cvg = HTSeq.GenomicArray({p.chrom: p.pos}, stranded=False, typecode="i")
        for read in file[window]:
            if read.proper_pair and read.aligned and not read.failed_platform_qc and read.optional_field('NH') == 1:
                cig_list = []
                for x in read.cigar:
                    cig_list.append(x.type)
                if 'N' not in cig_list:
                    for cigopt in read.cigar:
                        if cigopt.type == 'M':
                            cvg[cigopt.ref_iv] += 1
                            count += cigopt.ref_iv.length
        cvg.write_bedgraph_file(bedgraphFile)
        peak = int(math.ceil(count / p.pos))
        print("Time elapsed for coverage calculation of chromosome {}: {} seconds".format(p.chrom, time.time() - run_time))
        cvg_time = time.time()
        print("Start applying Flanking Region Filter for Alus on chromosome ".format(p.chrom))
        frf_stranded_bam(annotation, peak_plus, peak_minus, p.chrom)
        print("Time elapsed for Pol III Alus discovery on chromosome {}: {} seconds".format(p.chrom, time.time() - cvg_time))
        run_time = time.time()


# Calculate coverage in the Left/Center/Right arm of the SINEs

def frf_stranded(gtf, peak_plus, peak_minus, chr):
    global alu_list
    global count_sine
    for element in gtf:
        if element[0] == chr:
            if element[6] == '+':
                if sum(bw_plus.values(element[0], int(element[3]),int(element[4]))) > args.background * ((int(element[4]) - int(element[3])) * peak_plus):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(np.nan_to_num(bw_plus.values(element[0], int(element[9]),int(element[10])))))
                        right = int(sum(np.nan_to_num(bw_plus.values(element[0],int(element[10]), int(element[10]) + args.right_region))))
                        left = int(sum(np.nan_to_num(bw_plus.values(element[0], int(element[9]) - args.tss - args.left_region, int(element[9]) - args.tss))))
                        out = int(sum(np.nan_to_num(bw_plus.values(element[0], int(element[10]) + args.right_region, int(element[10]) + args.right_region + args.out_region))))
                    else:
                        continue
                else:
                    continue
            elif element[6] == '-':
                if sum(bw_minus.values(element[0], int(element[3]),int(element[4]))) > args.background * ((int(element[4]) - int(element[3])) * peak_minus):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(np.nan_to_num(bw_minus.values(element[0],int(element[9]),int(element[10])))))
                        right = int(sum(np.nan_to_num(bw_minus.values(element[0], int(element[9]) - args.right_region,int(element[9])))))
                        left = int(sum(np.nan_to_num(bw_minus.values(element[0], int(element[10]) + args.tss, int(element[10]) + args.tss + args.left_region))))
                        out = int(sum(np.nan_to_num(bw_minus.values(element[0], int(element[9]) - args.right_region - args.out_region, int(element[9]) - args.right_region))))
                    else:
                        continue
                else:
                    continue
            body_length = int(element[10]) - int(element[9])
            start = int(element[9])
            percentage = 0
            bases_over_bg = 0
            if left <= (central / (int(element[10])-int(element[9]))) * args.left_region * args.lratio and right < (central / (int(element[10])-int(element[9]))) * args.right_region and out <= (central / (int(element[10])-int(element[9]))) * args.out_region * args.oratio:
                while (percentage < args.fraction) and (start < int(element[10])):
                    if element[6] == '+':
                        if (bw_plus.stats(element[0], start, start + 1)[0]>args.background * peak_plus):
                            bases_over_bg += 1
                            percentage = bases_over_bg / body_length
                    elif element[6] == '-':
                        if (bw_minus.stats(element[0], start, start + 1)[0]>args.background * peak_minus):
                            bases_over_bg += 1
                            percentage = bases_over_bg / body_length
                    start += 1
                if percentage >= args.fraction:
                    alu_list.append([element[8], element[0], int(element[3]),int(element[4]), element[6], left, central, right, out,peak_plus, peak_minus, percentage])
        else:
            continue

def frf_unstranded(gtf, peak, chr):
    global count_sine
    for element in gtf:
        if element[0] == chr:
            if sum(bw.values(element[0], int(element[3]),int(element[4]))) > args.background * ((int(element[4]) - int(element[3])) * peak):
                count_sine += 1
                if "MIR" in element[8] or "Alu" in element[8]:
                    if element[6] == '+':
                        central = int(sum(bw.values(element[0], int(element[9]),int(element[10]))))
                        right = int(sum(bw.values(element[0],int(element[10]), int(element[10]) + args.right_region)))
                        left = int(sum(bw.values(element[0], int(element[9]) - args.tss - args.left_region,int(element[9]) - args.tss)))
                        out = int(sum(bw.values(element[0], int(element[10]) + args.right_region, int(element[10]) + args.right_region + args.out_region)))
                    elif element[6] == '-':
                        central = int(sum(bw.values(element[0],int(element[9]),int(element[10]))))
                        right = int(sum(bw.values(element[0], int(element[9]) - args.right_region,int(element[9]))))
                        left = int(sum(bw.values(element[0], int(element[10]) + args.tss, int(element[10]) + args.tss + args.left_region)))
                        out = int(sum(bw.values(element[0], int(element[9]) - args.right_region - args.out_region, int(element[9]) - args.right_region)))
                    body_length = int(element[10]) - int(element[9])
                    start = int(element[9])
                    percentage = 0
                    bases_over_bg = 0
                    if left <= (central / (int(element[10])-int(element[9]))) * args.left_region * args.lratio and right < (central / (int(element[10])-int(element[9]))) * args.right_region and out <= (central / (int(element[10])-int(element[9]))) * args.out_region * args.oratio:
                        while (percentage < args.fraction) and (start < int(element[10])):
                            if (bw.stats(element[0], start, start + 1)[0]>args.background * peak):
                                bases_over_bg += 1
                                percentage = bases_over_bg / body_length
                            start += 1
                        if percentage >= args.fraction:
                            alu_list.append([element[8], element[0], int(element[3]),int(element[4]), element[6], left, central, right, out, peak, percentage])
                else:
                    continue
            else:
                continue
        else:
                continue

def frf_stranded_bam(gtf, peak_plus, peak_minus, chr):
    global count_sine
    for element in gtf:
        if element[0] == chr:
            if element[6] == '+':
                if sum(list(cvg_plus[HTSeq.GenomicInterval(element[0], int(element[3]),int(element[4]))])) > args.background * (
                            (int(element[4]) - int(element[3])) * peak_plus):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(list(cvg_plus[HTSeq.GenomicInterval(element[0],int(element[9]),int(element[10]))])))
                        right = int(sum(list(cvg_plus[HTSeq.GenomicInterval(element[0],int(element[10]), int(element[10]) + args.right_region)])))
                        left = int(sum(list(cvg_plus[HTSeq.GenomicInterval(element[0], int(element[9]) - args.tss - args.left_region,int(element[9]) - args.tss)])))
                        out = int(sum(list(cvg_plus[HTSeq.GenomicInterval(element[0], int(element[10]) + args.right_region, int(element[10]) + args.right_region + args.out_region)])))
                    else:
                        continue
                else:
                    continue
            elif element[6] == '-':
                if sum(list(cvg_minus[HTSeq.GenomicInterval(element[0], int(element[3]),int(element[4]))])) > args.background * (
                            (int(element[4]) - int(element[3])) * peak_minus):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(list(cvg_minus[HTSeq.GenomicInterval(element[0],int(element[9]),int(element[10]))])))
                        right = int(sum(list(cvg_minus[HTSeq.GenomicInterval(element[0], int(element[9]) - args.right_region,int(element[9]))])))
                        left = int(sum(list(cvg_minus[HTSeq.GenomicInterval(element[0], int(element[10]) + args.tss, int(element[10]) + args.tss + args.left_region)])))
                        out = int(sum(list(cvg_minus[HTSeq.GenomicInterval(element[0], int(element[9]) - args.right_region - args.out_region, int(element[9]) - args.right_region)])))
                    else:
                        continue
                else:
                    continue
            body_length = int(element[10]) - int(element[9])
            start = int(element[9])
            percentage = 0
            bases_over_bg = 0
            if left <= (central / (int(element[10])-int(element[9]))) * args.left_region * args.lratio and right < (
                                central / (int(element[10])-int(element[9]))) * args.right_region and out <= (central / (int(element[10])-int(element[9]))) * args.out_region * args.oratio:
                while (percentage < args.fraction) and (start < int(element[10])):
                    if element[6] == "+":
                        if sum(list(cvg_plus[HTSeq.GenomicInterval(element[0], start, start + 1)]))>args.background * peak_plus:
                            bases_over_bg += 1
                            percentage = bases_over_bg / body_length
                    elif element[6] == "-":
                        if sum(list(cvg_minus[HTSeq.GenomicInterval(element[0], start, start + 1)]))>args.background * peak_minus:
                            bases_over_bg += 1
                            percentage = bases_over_bg / body_length
                    start += 1
                if percentage >= args.fraction:
                    alu_list.append([element[8], element[0], int(element[3]),int(element[4]), element[6], left, central, right, out, peak_plus, peak_minus,percentage])

        else:
            continue

def frf_unstranded_bam(gtf, peak, chr):
    global count_sine
    for element in gtf:
        if element[0] == chr:
            if element[6] == '+':
                if sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[3]),int(element[4]))])) > args.background * (
                            (int(element[4]) - int(element[3])) * peak):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(list(cvg[HTSeq.GenomicInterval(element[0],int(element[9]),int(element[10]))])))
                        right = int(sum(list(cvg[HTSeq.GenomicInterval(element[0],int(element[10]), int(element[10]) + args.right_region)])))
                        left = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[9]) - args.tss - args.left_region,int(element[9]) - args.tss)])))
                        out = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[10]) + args.right_region, int(element[10]) + args.right_region + args.out_region)])))
                    else:
                        continue
                else:
                    continue
            elif element[6] == '-':
                if sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[3]),int(element[4]))])) > args.background * ((int(element[4]) - int(element[3])) * peak):
                    count_sine += 1
                    if "MIR" in element[8] or "Alu" in element[8]:
                        central = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[9]), int(element[10]))])))
                        right = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[9]) - args.right_region,int(element[9]))])))
                        left = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[10]) + args.tss, int(element[10]) + args.tss + args.left_region)])))
                        out = int(sum(list(cvg[HTSeq.GenomicInterval(element[0], int(element[9]) - args.right_region - args.out_region, int(element[9]) - args.right_region)])))
                    else:
                        continue
                else:
                    continue
            body_length = int(element[10]) - int(element[9])
            start = int(element[9])
            percentage = 0
            bases_over_bg = 0
            if left <= (central / (int(element[10])-int(element[9]))) * args.left_region * args.lratio and right < (
                                central / (int(element[10])-int(element[9]))) * args.right_region and out <= (central / (int(element[10])-int(element[9]))) * args.out_region * args.oratio:
                while (percentage < args.fraction) and (start < int(element[10])):
                    if sum(list(cvg[HTSeq.GenomicInterval(element[0], start, start + 1)]))>args.background * peak:
                        bases_over_bg += 1
                        percentage = bases_over_bg / body_length
                    start += 1
                if percentage >= args.fraction:
                    alu_list.append([element[8], element[0], int(element[3]),int(element[4]), element[6], left, central, right, out, peak, percentage])

        else:
            continue


#@profile
def main():
    if args.filetype == 'bam':
        bamfile = HTSeq.BAM_Reader(args.coverage)
        if bamfile.get_header_dict()['HD']['SO'] != 'coordinate':
            print("BAM file must be sorted by position. Exiting....")
            raise SystemExit
        else:
            if args.stranded == 'auto':
                strdchk(bamfile)
            if args.stranded == 'reverse' or args.stranded == 'yes':
                print("Strandedness: {}".format(args.stranded))
                cvg_bam(bamfile)
            elif args.stranded == 'no':
                print("Strandedness: {}".format(args.stranded))
                cvg_bam_unstranded(bamfile)

    elif args.filetype == 'bw':
        if args.stranded == 'reverse' or args.stranded == 'yes':
            for chrom in chromosomes:
                peak_plus = int(math.ceil(bw_plus.stats(chrom[0])[0]))
                peak_minus = int(math.ceil(bw_minus.stats(chrom[0])[0]))
                print("Start applying Flanking Region Filter for chromosome {}".format(chrom[0]))
                frf_stranded(annotation, peak_plus, peak_minus, chrom[0])
                print("Time elapsed {}".format((time.time() - start_time)))
        elif args.stranded == 'no':
            for chrom in chromosomes:
                peak = int(math.ceil(bw.stats(chrom[0])[0]))
                print("Start applying Flanking Region Filter for chromosome {}".format(chrom[0]))
                frf_unstranded(annotation, peak, chrom[0])
                print("Time elapsed {}".format((time.time() - start_time)))

    # Write to output the Bona Fide SINEs

    with open(args.output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(alu_list)

    logfile.append(["Running parameters:","-g", args.ref_genome, "-s ", args.stranded, "-t ",args.filetype, "-bg",args.background, "-TSS", args.tss, "-LR",args.left_region, "-RR" ,args.right_region, "-OR ",args.out_region, "-lr", args.lratio, "-or",args.oratio, "-f", args.fraction])
    logfile.append(["Annotation file:", args.gtf])
    with open(args.output + ".log", 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(logfile)

    print("Finished!")
    print("Total Time elapsed {}".format((time.time() - start_time)))
    print("Number of SINEs over background {}".format(count_sine))
# TODO GDC commands

if __name__ == '__main__':
    main()
