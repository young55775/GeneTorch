import pandas as pd
import io
import os
import shutil


def search(df, col, kw):
    return df[col] == kw


class readfile:
    def __init__(self, path=None):
        self.path = path
        self.taglist = []
        self.names = []
        self.readfile()
        self.co_data = pd.DataFrame()
        self.result = pd.DataFrame()
        self.suppressor_group = []
        self.candidate = []

    def readfile(self):
        if self.path != None:
            files = os.listdir(self.path)
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
                ).rename(columns={'#CHROM': 'CHROM'})

            for file in files:
                if not os.path.isdir(file):
                    if os.path.basename(file).split('.')[-1] == 'vcf':
                        f = read_vcf(self.path + '/' + file)
                        filelist.append(f)
                        filename.append(os.path.basename(file).split('.')[0])
            self.names = filename

            def simp_file(raw_df):
                b = []
                for i, row in raw_df.iterrows():
                    b.append(row['INFO'].split('|'))
                c = []
                for i in range(len(b)):
                    c.append([b[i][3], b[i][4], b[i][1], b[i][9], b[i][10]])
                genelist = pd.DataFrame(c, columns=['gene', 'ID', 'type', 'base', 'protein'])
                return genelist

            simp_list = []
            for i in range(len(filelist)):
                simp_list.append(simp_file(filelist[i]))
            self.taglist = simp_list
            for i in range(len(self.taglist)):
                self.taglist[i].insert(loc=len(self.taglist[i].columns), column='tag', value=filename[i])


class getfile:
    def __init__(self, filepath, filenames):
        self.filepath = filepath
        self.filenames = filenames
        self.temp = filepath + '\\temp'
        self.taglist = []
        self.names = []
        self.splitfile()
        self.suppressor_group = []
        self.candidate = []

    def splitfile(self):
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
            ).rename(columns={'#CHROM': 'CHROM'})

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
                c.append([b[i][3], b[i][4], b[i][1], b[i][9], b[i][10]])
            genelist = pd.DataFrame(c, columns=['gene', 'ID', 'type', 'base', 'protein'])
            return genelist

        simp_list = []
        for i in range(len(filelist)):
            simp_list.append(simp_file(filelist[i]))
        self.taglist = simp_list
        for i in range(len(self.taglist)):
            self.taglist[i].insert(loc=len(self.taglist[i].columns), column='tag', value=self.names[i])


def get_impact(a):
    newtaglist = []
    for i in range(len(a.taglist)):
        n = a.taglist[i]
        lista = []
        listb = []
        splice_donor = n[n['type'] == 'splice_donor_variant&intron_variant']['protein'].tolist()
        splice_acceptor = n[n['type'] == 'splice_acceptor_variant&intron_variant']['protein'].tolist()
        for i in range(len(splice_donor)):
            lista.append('X_donor')
        for i in range(len(splice_acceptor)):
            listb.append('X_acceptor')
        c = pd.DataFrame(n[n['type'] == 'splice_donor_variant&intron_variant'])
        d = pd.DataFrame(n[n['type'] == 'splice_acceptor_variant&intron_variant'])
        c['protein'] = lista
        d['protein'] = listb
        n = pd.concat([n, c, d])
        indp = search(n, 'protein', '')
        indq = search(n, 'protein', 'nan')
        ind = search(n, 'type', 'synonymous_variant')
        m = n.loc[ind, :]
        b = n.loc[indp, :]
        e = n.loc[indq, :]
        n = pd.concat([m, n, e]).drop_duplicates(keep=False)
        total = pd.concat([n, b]).drop_duplicates(keep=False)
        newtaglist.append(total)
    a.taglist = newtaglist
    return a.taglist
