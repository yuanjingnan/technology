from optparse import OptionParser
import os


def genetoid(g2id):
    gene2id={}
    with open(g2id,'r') as fr:
        fr.readline()
        for line in fr.readlines():
            genesymbol=line.strip().split("\t")
            gene2id[genesymbol[0]]=genesymbol[1]
    return gene2id


def trans_old_format(filename,gene2id):
    if filename.endswith(".gene.fpkm.xls"):
        # change the filename to *.new, and write the result to original filename
        newfilename=filename+".new"
        fw = open(newfilename,'w')
        with open(filename,"r") as fr:
            fw.write(fr.readline().strip()+"\t"+"Symbol\n")
            for line in fr.readlines():
                geid = line.strip().split("\t")[0]
                Symbol=''
                if geid in gene2id:
                    Symbol = gene2id[line.strip().split("\t")[0]]
                fw.write(line.strip()+"\t"+Symbol+"\n")
        fw.close()
    elif filename.endswith(".transcript.fpkm.xls"):
        # change the filename to *.new, and write the result to original filename
        newfilename = filename + ".new"
        fw = open(newfilename, 'w')
        with open(filename, "r") as fr:
            fw.write(fr.readline().strip()+"\t"+"Symbol\n")
            for line in fr.readlines():
                nmid = line.strip().split("\t")[1]
                Symbol=''
                if nmid in gene2id:
                    Symbol = gene2id[line.strip().split("\t")[1]]
                fw.write(line.strip()+"\t"+Symbol+"\n")
        fw.close()


def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--infile", dest="infile", help="the input file of geneExp result, *.gene.fpkm.xls or *.transcript.fpkm.xls")
    parser.add_option("-s", "--g2id", dest="g2id",default="/ldfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/mRNAseq_V1/RNA_RNAref_2018a/DataBase/hg19/rsem-build_refMrna/human.gene2symbol.txt",
                      help="the file like /ldfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/mRNAseq_V1/RNA_RNAref_2018a/DataBase/hg19/rsem-build_refMrna/human.gene2symbol.txt")
    (options, args) = parser.parse_args()
    infile = options.infile
    g2id=options.g2id
    gene2id = genetoid(g2id)
    for root,dirs,files in os.walk(infile):
        for filename in files:
            if filename.endswith(".gene.fpkm.xls"):
                genefpkm = os.path.join(root,filename)
                trans_old_format(genefpkm,gene2id)
            elif filename.endswith(".transcript.fpkm.xls"):
                transcriptfpkm = os.path.join(root,filename)
                trans_old_format(transcriptfpkm, gene2id)
            else:
                continue


if __name__ == "__main__":
    main()
