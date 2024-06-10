import pysam
import sys
try:
    import argparse
except:
    sys.path.insert(0,'/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda3/lib/python3.7')


def map_count(bamfiles,seqtype):
    total_reads = 0
    total_basepairs = 0
    mapped = 0
    unmapped=0
    Perfect_match=0
    mismatch = 0
    unique_match = 0
    multi_match = 0
    linenum=0
    for bamfile in bamfiles:
        if bamfile.endswith(".bam"):
            bamreader = pysam.AlignmentFile(bamfile, "rb")
        elif bamfile.endswith(".cram"):
            bamreader =pysam.AlignmentFile(bamfile,'rc')
        else:
            bamreader=pysam.AlignmentFile(bamfile,'r')
        for reads in bamreader:
            linenum += 1
            length = reads.rlen
    if seqtype=='PE':
        total_reads = linenum*2
    else:
        total_reads = linenum
    print(total_reads)
    print(length)



def main():
    parser = argparse.ArgumentParser(description="mapping stat from bam or sam")
    parser.add_argument("-b", "--bam", dest="bam", default="", help="bowtie bam result(or sam), separated by comma ','")
    parser.add_argument("-s","--strand",dest="strand", default="PE",help="stat strand information")
    parser.add_argument("-k", "--key", dest="key",default="", help="prefix of output")
    parser.add_argument("-t", "--seqType", dest="seqType", default="", help="sequencing type,PE or SE")
    args = parser.parse_args()
    bamfiles = [i.strip() for i in args.bam.split(",")]
    map_count(bamfiles,args.seqType)


if __name__ == "__main__":
    main()
