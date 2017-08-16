from __future__ import division
import pysam
from Bio import pairwise2
import argparse
from string import maketrans

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('bamfile',help='Consensus BAM file')
parser.add_argument('samplename',help='Sample name')
parser.add_argument('-r',"--reference",default="/gscuser/dspencer/refdata/GRCh37/all_sequences.fa",
                    help='reference genome fasta file')
parser.add_argument('-w',"--window",type=int,default=150,
                    help='window for creating ref and alt sequences')
parser.add_argument('-m',"--minreads",type=int,default=4,
                    help='minimum reads for defining an amplicon')
parser.add_argument('-v',"--minreads4vaf",type=int,default=4,
                    help='minimum reads in an amplicon for amplicon-based vaf calculation')

args = parser.parse_args()

mysample = args.samplename

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)
# open bam file
samfile = pysam.AlignmentFile(args.bamfile,"rb")
# open reffasta
fa = pysam.FastaFile(args.reference)

# window on either side of a variant for indel annotation 
window = args.window
# minimum number of read pairs with a unique insert length to define an amplicon
minampreads = args.minreads
# minimum number of reads in an amplicon needed to be eligible to calculate a VAF
minampreadsforvaf = args.minreads4vaf

vcffile.header.formats.add("TAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position')
vcffile.header.formats.add("SAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position that support the alternate allele')
vcffile.header.formats.add("VAFTYPE", 1, 'String', 'Either MAX (max amplicon VAF) or MEAN (mean of all reads)')
vcffile.header.formats.add("AMPS", 1, 'String', 'Amplicon string, in the format amplicon id 1,num. reference alleles in amplicon 1, num. alternate alleles in amplicon 1;amplicon id 2,num ref, num alt;...')
vcffile.header.formats.add("CVAF", 1, 'Float', 'Calculated VAF, which is either the maximum per-amplicon VAF (among amplicons with >' + str(minampreadsforvaf) + ' reads), or the mean for all reads (if no available amplicons or the max is 1.0 and the mean is <1.0')
vcffile.header.add_line("##ampliconparseroptions={'window':"+str(window)+",'minampreads':"+str(minampreads)+",'minampreadsforvaf':"+str(minampreadsforvaf)+"}")

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print "\n".join(hdr) + "\n" + ("\t".join((hd.split("\t"))[0:9])) + "\t" + mysample

for rec in vcffile.fetch():

    amps = {'ref':[],'alt':[]}
    ampcnt = {}    
    ampstring = list()
    ampnum = 0
    ampvafs = list()

    # if the variant is a substitution
    if len(rec.ref) == len(rec.alts[0]) and len(rec.ref) == 1:
        
        for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos):            
            if pileup.pos == rec.pos-1: # only consider the variant position

                for read in pileup.pileups:

                    # skip positions with indels
                    if not read.is_del and not read.is_refskip:

                        # get the amplicon using the template length and read 1 strand
                        ampname = abs(read.alignment.template_length)
                        if (read.alignment.is_read1 and read.alignment.is_reverse) or (read.alignment.is_read2 and not read.alignment.is_reverse):
                            ampname = -ampname

                        # count reads per amplicon
                        if ampname not in ampcnt:
                            ampcnt[ampname] = 1
                        else:
                            ampcnt[ampname] += 1

                        # count alleles by amplicon
                        if read.alignment.query_sequence[read.query_position] == rec.ref:
                            amps['ref'].append(ampname)
                        elif read.alignment.query_sequence[read.query_position] == rec.alts[0]:
                            amps['alt'].append(ampname)
                        
    else: # if indel

        refseqstart = rec.pos-window-1
        refseqend = rec.pos+window
        refseq = fa.fetch(rec.contig,refseqstart,refseqend)
        altseq = fa.fetch(rec.contig,refseqstart,rec.pos-1) # 0-based

        # if deletion
        if len(rec.ref) > len(rec.alts[0]):
            altseq = altseq + fa.fetch(rec.contig,rec.pos+len(rec.ref)-1,rec.pos+window+len(rec.ref))

        else: # if insertion
            altseq = altseq + rec.alts[0] + fa.fetch(rec.contig,rec.pos,rec.pos+window)

        reads = {}
        for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos):
            if pileup.pos == rec.pos-1:
                for read in pileup.pileups:

                    ampname = abs(read.alignment.template_length)
                    if (read.alignment.is_read1 and read.alignment.is_reverse) or (read.alignment.is_read2 and not read.alignment.is_reverse):
                        ampname = -ampname
            
                    if ampname not in ampcnt:
                        ampcnt[ampname] = 1
                    else:
                        ampcnt[ampname] += 1
                    
                    # if the read has fewer mismatches (NM tag) than the indel length and has no softclipped bases
                    # then assign to reference allele
                    if len(read.alignment.cigartuples) == 1 and read.alignment.cigartuples[0][0] == 0 and read.alignment.cigartuples[0][1] == read.alignment.query_length and read.alignment.get_cigar_stats()[0][10] < abs(len(rec.ref) - len(rec.alts[0])):
                        amps['ref'].append(ampname)

                    # if the record is a deletion and the read has a deletion at this position
                    # with the proper length then assign to the alt allele
                    elif len(rec.ref) > len(rec.alts[0]) and read.indel < 0 and abs(read.indel) == len(rec.ref)-len(rec.alts[0]):
                        amps['alt'].append(ampname)
                        
                    # if the record is an insertion and the read has an insertion
                    # with the proper bases and length then assign to the alt allele
                    elif len(rec.ref) < len(rec.alts[0]) and read.indel > 0 and rec.alts[0] == read.alignment.seq[read.query_position:read.query_position+read.indel+1]:
                        amps['alt'].append(ampname)
                    # if none of the above are satisified, use the alignment method
                    # to assign reads to alleles
                    else:
                        rdseq = read.alignment.query_sequence
                            
                        alnref = pairwise2.align.localms(rdseq,refseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                        alnalt = pairwise2.align.localms(rdseq,altseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                    
                        if alnref >= alnalt:
                            amps['ref'].append(ampname)
                        else:
                            amps['alt'].append(ampname)

    # count number of reads per amplicon and stor per amplicon data if there are enough
    ampwithalt = 0
    for a in ampcnt:
        if ampcnt[a] >= minampreads:
            ampnum += 1
            if amps['alt'].count(a) > 0:
                ampwithalt += 1
                
            ampstring.append(str(a) + "," + str(amps['ref'].count(a)) + "," + str(amps['alt'].count(a)))
            if amps['ref'].count(a) + amps['alt'].count(a) > minampreadsforvaf:
                ampvafs.append(amps['alt'].count(a) / (amps['ref'].count(a) + amps['alt'].count(a)))

    myvaf = 'N/A'
    vaftype = 'MAX'
    if len(ampvafs)>0 and not (max(ampvafs) == 1.0 and len(amps['ref']) > 0):
        myvaf = max(ampvafs)
    else:
        vaftype = 'MEAN'
        myvaf = (len(amps['alt']) / (len(amps['ref']) + len(amps['alt'])))

    mysample=0
    mygt = (0,1)
    for s in rec.samples:
        if rec.samples[s]['GT'] == (1,1) and myvaf > .99:
            mygt = (1,1)
 
    rec.samples[mysample]['GT'] = mygt
    rec.samples[mysample]['NR'] = len(amps['alt']) + len(amps['ref'])
    rec.samples[mysample]['NV'] = len(amps['alt'])
    rec.samples[mysample]['TAMP'] = ampnum
    rec.samples[mysample]['SAMP'] = ampwithalt
    rec.samples[mysample]['VAFTYPE'] = vaftype
    rec.samples[mysample]['AMPS'] = ';'.join(ampstring)
    rec.samples[mysample]['CVAF'] = myvaf

    print "\t".join(str(rec).rstrip().split("\t")[0:10])

# end vcf fetch

fa.close()
vcffile.close()
samfile.close()
