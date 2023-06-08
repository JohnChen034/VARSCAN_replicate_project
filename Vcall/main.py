import pandas as pd
import re
import argparse
import math
import numpy as np
from scipy.stats import fisher_exact
import os

def main():
    
    parser = argparse.ArgumentParser(description='variant calling')
    parser.add_argument('mpileupfile', metavar='MPILEUPFILE', type=str,
                        help='path to the input mpileup file')
    parser.add_argument('-mvf', '--min-var-freq', metavar='MINVARFREQ', type=float, default = 0.2,
                        help='minimum frequency threshold for variants (default: 0.2)')
    parser.add_argument('-mc', '--min-coverage', metavar='MINCOVERAGE', type=float, default = 8,
                        help='Minimum read depth at a position to make a call (default: 8)')
    parser.add_argument('-mr', '--min-reads', metavar='MINREADS', type=float, default = 2,
                        help='Minimum supporting reads at a position to call variants (default: 2)')
    parser.add_argument('-maq', '--min-avg-qual', metavar='MINAVGQUAL', type=float, default = 15,
                        help='Minimum base quality at a position to count a read (default: 15)')
    parser.add_argument('-mffh', '--min-freq-for-hom', metavar='MINFREQFORHOM', type=float, default = 0.75,
                        help='Minimum frequency to call homozygote (default: 0.75)')
    parser.add_argument('-p', '--p-value', metavar='PVALUE', type=float, default = 0.99,
                        help='p-value threshold for calling variants (default: 0.99)')
    parser.add_argument('-sf', '--strand-filter', metavar='STRANDFILTER', type=bool, default = True,
                    help='Ignore variants with >90% support on one strand (default: True)')
    parser.add_argument('-out', '--out-path', metavar='OUTPATH', type=str, default = os.getcwd(),
                help='Specify the output path (default: current working directory)')
    parser.add_argument('-f', '--file-name', metavar='FILENAME', type=str, default='output',
                    help='Specify the file name (default: output)')
    
    args = parser.parse_args()
    
    stuff = vars(args)
    
    ###input params
    mpileup_file = stuff["mpileupfile"]
    min_var_frequency = stuff["min_var_freq"] #set minimum var frequency threshold
    min_reads = stuff["min_reads"] #min ALT count filter
    min_coverage = stuff["min_coverage"] #min coverage
    min_avg_qual = stuff["min_avg_qual"] #min quality score
    p_val_threshold = stuff["p_value"] #p-value threshold
    homo_threshold = stuff["min_freq_for_hom"] #Minimum frequency to call homozygote
    if_strand_filt = stuff["strand_filter"] #if perform strand filter
    output_path = stuff["out_path"] #output path
    file_name = stuff["file_name"] #output filename
    
    ###functions 
    def get_INFO(format_row):
        ADP = int(np.mean(format_row.apply(lambda x: int(x[3]))))
        counts = format_row.apply(lambda x: x[0]).values
        WT = np.count_nonzero(counts == '0/0')
        HET = np.count_nonzero(counts == '0/1')
        HOM = np.count_nonzero(counts == '1/1')
        return f"ADP={ADP};WT={WT};HET={HET};HOM={HOM};NC=0"
    
    def filter_strand(row):
        out = []
        for n in range(4, row.shape[0], 3):
            filt_string = row[n]
            ADR = len(re.findall(r'[acgt]', row[n]))
            ADF = len(re.findall(r'[ACGT]', row[n]))
            total = ADR + ADF
            if total == 0:
                out += [True]
            else:
                out += [ADR / total > 0.9 or ADF / total > 0.9]
        return all(out)

    def convert_format(format_lst):
        format_lst[1] = int(format_lst[1])
        format_lst[7] = f"{format_lst[7]:.4e}"
        return ':'.join(str(elem) for elem in format_lst)

    def filter_pval(df, format_df, threshold = p_val_threshold):
        return df[format_df.applymap(lambda x: x[7] <= threshold).any(axis=1)]

    def create_format(row):
        end = row.shape[0]
        count = 1
        for n in range(4, end, 3):
            raw_string = row[n]
            raw_quality = row[n+1]
            if len(raw_string) == 0:
                row[f"Sample{count}"] = [1]*14
            else:
                row[f"Sample{count}"] = FORMAT(raw_string, raw_quality)

            count+=1
        out = row[end:]
        return out

    def FORMAT(raw_string, raw_quality, homo_threshold = homo_threshold, min_avg_qual = min_avg_qual, min_reads = min_reads, min_var_frequency = min_var_frequency):
        filt_string = ''.join([raw_string[i] for i in range(len(raw_quality)) if ord(raw_quality[i]) >= min_avg_qual+33])
        filt_quality = ''.join([i for i in raw_quality if ord(i) >= min_avg_qual + 33])

        SDP = len(raw_string)
        DP = len(filt_string)

        RDF = filt_string.count('.')
        RDR = filt_string.count(',')
        ADR = len(re.findall(r'[acgtn]', filt_string))
        ADF = len(re.findall(r'[ACGTN]', filt_string) )       
        RD = RDF + RDR
        AD = DP - RD
        FREQ = f"{(AD / DP) * 100:.2f}%"
        GT = '0/0' if RD/DP >= homo_threshold else '1/1' if AD/DP >= homo_threshold else '0/1'

        def fischer_test(DP, RD, AD):
            contingency_table = [[0, DP], [AD, RD]]
            p_value = fisher_exact(contingency_table, alternative='less')[1]
            return p_value

        PVAL = fischer_test(DP, RD, AD)

        GQ = -10 * math.log10(PVAL)

        def RABQ(string, quality):
            if RD == 0:
                RBQ = 0
            else:
                indices = [i for i in range(len(string)) if string[i] not in r'[ACGTNacgtn]']
                RBQ = np.mean([ord(quality[i])-33 for i in indices]) if len(indices) > 0 else 0

            if AD == 0:
                ABQ = 0
            else:   
                indices = [i for i in range(len(string)) if string[i] in r'[ACGTNacgtn]']
                ABQ = np.mean([ord(quality[i])-33 for i in indices]) if len(indices) > 0 else 0
            return int(RBQ), int(ABQ)

        RBQ, ABQ = RABQ(filt_string, filt_quality)
        return [GT,GQ,SDP,DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR]

    def filter_base_quality(row, threshold = min_avg_qual):
        for n in range(4, row.shape[0], 3):
            revised_bases = ''
            revised_qualities = ''
            for i in range(row[n-1]):
                if ord(row[n+1][i]) >= threshold + 33:
                    revised_bases += row[n][i]
                    revised_qualities += row[n+1][i]
            row[n+1] = revised_qualities
            row[n] = revised_bases
        return row

    def filter_gap(row):
        for n in range(4, row.shape[0], 3):
            revised_bases = ''
            revised_qualities = ''
            for i in range(row[n-1]):
                if row[n][i] not in "<>":
                    revised_bases += row[n][i]
                    revised_qualities += row[n+1][i]
            row[n+1] = revised_qualities
            row[n] = revised_bases
        return row

    def filter_read_depth(df, threshold = min_coverage):
        filtered_df = df[~(df[range(4, df.shape[1], 3)].apply(lambda x: x.str.len()) < threshold).all(axis=1)]
        return filtered_df



    #helper function for min_var_filt
    def single_var_filt(string, min_var_frequency = 0.2, min_reads = 2):
        lst = []
        varrs = np.unique(re.findall("[ACGTNacgtn]", string))
        for i in varrs:
            count_var = string.count(i)
            ifup = i.isupper()
            if ifup:
                freq = count_var/len(re.findall('[ATCG.]', string))
            else:
                freq = count_var/len(re.findall('[atcg,]', string))
            boo = freq < min_var_frequency and count_var < min_reads
            if boo:
                string = string.replace(i, ".") if ifup else string.replace(i, ",")
            lst += [boo]
        return (string, lst)

    def min_var_filt(df):
        outdf = pd.DataFrame()
        for i in range(4, df.shape[1], 3):
            tup = df[i].apply(single_var_filt)
            new_strs = tup.apply(lambda x: x[0])
            boo_lst = tup.apply(lambda x: x[1])
            df[i] = new_strs
            outdf[i] = boo_lst.apply(lambda x: all(x))
        df = df[~outdf.all(axis=1)]
        return df
    
    #0. Read the mpileup data into a pandas DataFrame
    #print(mpileup_file)
    df = pd.read_csv(mpileup_file, sep='\t', header=None)
    col_len = len(df.columns)


    ### filter process

    ##1. basic filters

    #1.1 remove $^ signs
    for i in range(4, col_len, 3):
        df[i] = df[i].str.replace('\^[\x21-\x7E]', '', regex=True)
        df[i] = df[i].str.replace('$', '', regex=True)

    #1.2 select only positions with 1+ variants read
    df = df[df[range(4, col_len, 3)].apply(lambda x: 
                                               ''.join(x.astype(str)), axis=1).str.contains('[ACGTNacgtn]')]


    #1.3 ignore the positions with insertions and deletions
    for i in ['-[0-9]+[ACGTNacgtn]+', "\+[0-9]+[ACGTNacgtn*#]+", r'[*]']:
        for j in range(4, col_len, 3):
            df = df[~df[j].str.contains(i)]
    
    
    #1.2 select only positions with 1+ variants read
    df = df[df[range(4, col_len, 3)].apply(lambda x: 
                                               ''.join(x.astype(str)), axis=1).str.contains('[ACGTNacgtn]')]
    print('basic filters done')

    #->3.2 prep for pvalue filter and output
    format_df = df.apply(create_format, axis=1)
    print("format_df done")

    ## layer 1
    #2.1 min read depth filter
    df = filter_read_depth(df)
    print("filter_read_depth done")



    #2.2 filter base quality 
    df = df.apply(filter_base_quality, axis =1)
    print("filter_base_quality done")


    #ps: repeat 1.2 select only positions with 1+ variants read
    df = df[df[range(4, col_len, 3)].apply(lambda x: 
                                               ''.join(x.astype(str)), axis=1).str.contains('[ACGTNacgtn]')]
    print("repeat done")
    
    ## Layer 2
    #3.1 filter based on minimum variant frequency
    df = min_var_filt(df)
    print("min_var_filt done")

    format_df = format_df.loc[df.index.tolist()]

    print("filt format_df done")

    #3.2 filter based on min p value
    df = filter_pval(df, format_df)
    print("filter_pval done")

    #3.3 filter based on variants percentage on strand
    if if_strand_filt:
        df = df[~df.apply(filter_strand, axis=1)]
        print('strand filter done')

    #output 
    print("start generating output VCF")
    outdf = pd.DataFrame()
    outdf[['#CHROM','POS']] = df[[0, 1]]
    outdf[['ID']] = '.'
    outdf['REF'] = df[2].apply(lambda x: x.upper())
    outdf['ALT'] = df[range(4, col_len, 3)].apply(lambda x:
                "".join(x.astype(str)), axis=1).str.findall('[ACGTNacgt]').apply(lambda x: x[0].upper())
    outdf['QUAL'] = '.'
    outdf['FILTER'] = 'PASS'

    format_df = format_df.loc[df.index.tolist()]
    # INFO column: 
    outdf['INFO'] = format_df.apply(get_INFO, axis=1)
    
    outdf['FORMAT'] = "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"
    format_df = format_df.applymap(convert_format)
    outdf = pd.concat([outdf, format_df], axis=1)
    
    content = '''##fileformat=VCFv4.1
##source=Vcall
##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 15">
##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">
##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
'''
    file_path = os.path.join(output_path, f"{file_name}.vcf")
    with open(file_path, 'w') as file:
        file.write(content)
        file.write(outdf.to_csv(index=False, sep='\t'))

print(f"successfully output in {file_path}")
if __name__ == '__main__':
    main()