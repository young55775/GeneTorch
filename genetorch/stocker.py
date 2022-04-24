import pandas as pd
import io
import os
import shutil


class stockfile:
    def __init__(self, filepath, filenames, outpath):
        self.filepath = filepath
        self.filenames = filenames
        self.temp = filepath + '\\temp'
        self.outpath = outpath
        self.simplist = []
        self.names = []
        self.stock()

    def stock(self):
        os.mkdir(self.temp)
        files = []
        a = os.listdir(self.filepath)
        for i in range(len(a)):
            if os.path.isdir(self.filepath + '/' + a[i]):
                files.append(a[i])

        self.names = []
        for i in range(len(files)):
            self.names.append(files[i].split('_')[0])

        files2 = []
        for i in range(len(files)):
            files2.append(os.listdir(self.filepath + '/' + files[i]))
            for j in range(len(files2[i])):
                if files2[i][j] == self.filenames:
                    shutil.copy(self.filepath + '/' + files[i] + '/' + self.filenames,
                                self.temp + '/' + self.names[i] + '.vcf')

        files3 = os.listdir(self.temp)
        filelist = []
        filename = []

        # read vcf from file using pandas dataframe
        def read_vcf(path_a):
            with open(path_a, 'r') as f:
                lines = [l for l in f if not l.startswith('##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t'
            )

        for file in files3:
            if not os.path.isdir(file):
                f = read_vcf(self.temp + '/' + file)
                filelist.append(f)
                filename.append(os.path.basename(file).split('.')[0])

        # read info
        def simp_file(raw_df):
            b = []
            for i, row in raw_df.iterrows():
                b.append(row['INFO'].split('|'))
            c = []
            for i in range(len(b)):
                c.append([b[i][9], b[i][10], b[i][1], b[i][2], b[i][3], b[i][4], b[i][6], '', ''])
            genelist = pd.DataFrame(c, columns=['nucleotide', 'amino acid', 'type', 'impact', 'gene', 'WBGeneID',
                                                'transcript', 'allele', 'comments'])
            return genelist

        simp_list = []
        for i in range(len(filelist)):
            a = simp_file(filelist[i])
            b = filelist[i].drop(['INFO', 'FILTER', 'ID', 'REF', 'ALT', 'QUAL','FORMAT','unknown'], axis=1)
            c = pd.concat([b, a], axis=1, join='outer')
            simp_list.append(c)
        self.simplist = simp_list

        for i in range(len(self.simplist)):
            self.simplist[i].to_csv(self.outpath + '/' + self.names[i] + '.csv', index=False, sep='\t')
        print('Output csv is in ' + self.outpath)
