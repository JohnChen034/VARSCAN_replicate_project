import pandas as pd
import re
import argparse

def filter_base_quality(row):
    for n in range(4, row.shape[0], 3):
        revised_bases = ''
        revised_qualities = ''
        threshold = 30
        for i in range(row[n-1]):
            if ord(row[n+1][i]) >= threshold + 33:
                revised_bases += row[n][i]
                revised_qualities += row[n+1][i]
        row[n+1] = revised_qualities
        row[n] = revised_bases
    return row

def main():
    
    parser = argparse.ArgumentParser(description='variant calling')
    parser.add_argument('mpileupfile', metavar='MPILEUPFILE', type=str,
                        help='path to the input mpileup file')
    parser.add_argument('-mvf', '--minvarfreq', metavar='MINVARFREQ', type=float, default = 0.2,
                        help='minimum frequency threshold for variants (default: 0.2)')
    
    args = parser.parse_args()
    
    stuff = vars(args)
    
    #print(stuff)
    
    # Read the mpileup data into a pandas DataFrame
    mpileup_file = stuff["mpileupfile"]
    #print(mpileup_file)
    df = pd.read_csv(mpileup_file, sep='\t', header=None)
    df_original = df

    col_len = len(df.columns)
    
    # filter process
    df = df_original
    #1. remove $^ signs
    for i in range(4, col_len, 3):
        df[i] = df[i].str.replace('\^[\x21-\x7E]', '', regex=True)
        df[i] = df[i].str.replace('$', '', regex=True)
    
    #2. select only positions with 1+ variants read
    df = df[df[range(4, col_len, 3)].apply(lambda x: 
                                               ''.join(x.astype(str)), axis=1).str.contains('[ACGTNacgtn]')]


    #3. ignore the positions with insertions and deletions
    for i in ['-[0-9]+[ACGTNacgtn]+', "\+[0-9]+[ACGTNacgtn*#]+", r'[*]']:
        for j in range(4, col_len, 3):
            df = df[~df[j].str.contains(i)]
    df_pre_filt = df
    
    df = df.apply(filter_base_quality, axis =1)
    df_quality_filt = df
    
    df = df_quality_filt
    #repeat 2. select only positions with 1+ variants read
    df = df[df[range(4, col_len, 3)].apply(lambda x: 
                                               ''.join(x.astype(str)), axis=1).str.contains('[ACGTNacgtn]')]
      
    
    #4. set minimum var frequency and get rid of positoins with variation probabilty less than this threshold
    min_var_frequency = 0.2

    def single_var_filt(string):
        pattern = re.compile(r'[ACGTNacgtn]')
        matches = pattern.findall(string)
        return (len(matches) / len(string)) > min_var_frequency

    def min_var_filt(df):
        df_copy = df.copy()
        for i in range(4, col_len, 3):
            df_copy[i] = df_copy[i].apply(single_var_filt)
        return df[df_copy[range(4, col_len, 3)].apply(lambda x: any(x.values), axis=1)]

    outdf = min_var_filt(df)
    
    outdf[3] = outdf[range(4, col_len, 3)].apply(lambda x:
                "".join(x.astype(str)), axis=1).str.findall('[ACGTNacgt]').apply(lambda x: x[0].upper())
    outdf[2] = outdf[2].apply(lambda x: x.upper())
    outdf = outdf[[0, 1, 2, 3]]
    
    outdf.to_csv("notcompleted.vcf", sep = "\t", index = False, header = False)
    
if __name__ == '__main__':
    main()