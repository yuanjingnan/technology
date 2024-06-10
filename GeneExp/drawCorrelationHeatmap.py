import os
import sys
import pandas as pd
try:
    import argparse
except:
    sys.path.insert(0,'/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda3/lib/python3.7')
Rscript="/ldfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/mRNAseq_V1/RNA_RNAref_2018a/Bin_CentOS6/R/Rscript"
Convert = "/ldfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/mRNAseq_V1/RNA_RNAref_2018a/Bin_CentOS6/convert"


def get_fpkm(indir,algorithm):
    all_sample_fpkm = pd.DataFrame()
    for fpkm_file in os.listdir(indir):
        if fpkm_file.endswith(".gene.fpkm.xls"):
            sample = fpkm_file.replace(".gene.fpkm.xls","")
            all_path_fpkm = os.path.join(indir,fpkm_file)
            samples_fpkm = pd.read_csv(all_path_fpkm,sep="\t",usecols=['gene_id','FPKM'])
            # gene_id, set it as row_index, change colnames
            samples_fpkm = samples_fpkm.rename(columns={'FPKM':sample}).set_index('gene_id')
            if all_sample_fpkm.empty:
                all_sample_fpkm = samples_fpkm
            else:
                all_sample_fpkm = pd.merge(all_sample_fpkm,samples_fpkm,how='outer',on='gene_id')
        else:
            continue
    # na to 0.01
    all_sample_fpkm = all_sample_fpkm.fillna(0.01)
    all_samples_cor = all_sample_fpkm.corr(method=algorithm).round(decimals=7)
    return all_samples_cor


def plot_heatmap(all_samples_cor,algorithm,outdir):
    if all_samples_cor.shape[1]>6:
        cellsize = all_samples_cor.shape[1]
    else:
        cellsize = 8
    legendName = algorithm+"_value"
    fw = open(os.path.join(outdir,"correlation-heatmap.R"),'w')
    rcode='''
    library(reshape2)
    library(ggplot2)
    library(RColorBrewer)
    x <- read.table("{0}/AllSamples.correlation.xls", sep = "\t", head = T)
    xx = as.matrix(x[,-1])
    rownames(xx) = names(x)[-1]
    xx = melt(xx)
    names(xx)=c("Var1","Var2","{1}");
    pdf("{0}/AllSamples.CorrelationHeatmap.pdf",width={2},height={2})
    ggplot(xx, aes(Var1, Var2, fill={1}))+
    #geom_tile(width=0.8, height=0.8)+
    geom_tile(color='black')+
    geom_text(label=round(xx${1}, 3))+
    scale_fill_gradient(low='#DEEBF7',high='#08519C')+
    theme(axis.text = element_text(angle=30, hjust=1,size=11,vjust=0,color='black'),
    panel.background = element_rect(fill='transparent'),
    panel.grid=element_line(color='grey'),legend.title = element_text(size = 13))+ 
    labs(x="",y="")
    dev.off()
    '''.format(outdir,legendName,cellsize)
    fw.write(rcode)
    fw.close()
    cmd1 = "{0} {1}/correlation-heatmap.R".format(Rscript,outdir)
    cmd2 = "{0} -density 300 -resize 40% {1}/AllSamples.CorrelationHeatmap.pdf {1}/AllSamples.CorrelationHeatmap.png".format(Convert,outdir)
    os.system(cmd1)
    os.system(cmd2)


def main():
    parser = argparse.ArgumentParser(description="Draw Correlation Heatmap between Samples")
    parser.add_argument("-i", "--indir", dest="indir", default="",
                        help="input directory, containing *.Gene.rpkm.xls\(gene expression level files\)")
    parser.add_argument("-R","--Rpath",dest="Rpath",default="",help="path of Rscript, default:")
    parser.add_argument("-c", "--convert", dest="convert", default="", help="path of convert, default:")
    parser.add_argument("-a", "--algorithm", dest="algorithm", default="pearson", help="algorithm for correlation, default: pearson")
    parser.add_argument("-o", "--outdir", dest="outdir", default="", help="output directory, default: current directory")
    opt = parser.parse_args()
    all_samples_cor = get_fpkm(opt.indir,opt.algorithm)
    all_samples_cor.to_csv(os.path.join(opt.outdir,"AllSamples.correlation.xls"),sep="\t")
    plot_heatmap(all_samples_cor,opt.algorithm,opt.outdir)


if __name__ == "__main__":
    main()