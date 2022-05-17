# Gene Torch 1.2.0

[TOC]

## What is new in version 1.2.0
```genetorch.finder.plotpro(readfile,size,c='Reds',sub_c='RdPu',origin = a)```

A more professional version to distinguish gene via binominal distribution instead of overlap. The calculation is a bit slow but more accurate.
The genes you tried will be saved ```in a.candidate``` and all the sample names of the candidate gene will be listed in ```a.suppressor_group```.
You should also notice that in the previous version, ```genetorch.finder.plot(result,size,c,sub_c)``` you shall input a result file. In this function, you may need to input the whole readfile() object or getfile() object instead.

```genetorch.reader.get_impact(readfile)```

This function can now annotate splice_donor and splice_acceptor variation with 'X_donor' and 'X_acceptor' because we noticed that splice variation do have an impact on protein function.

```genetorch.finder.get_p_sup(gene,readfile)```

return a Dataframe including all the p values the gene you chose vs others.

```genetorch.finder.get_p_between(readfile,genea,geneb)```

return the p value between genea and geneb

```genetorch.finder.find(readfile)/filter(readfile,lengthlimit,rid)```

We canceled the threshold in 'filter' because we noticed that it is not very useful and can make our users confused. These two functions will automatically do a ```genetorch.reader.get_impact``` now.
Also, in some cases you may have to reset the  candidate and suppressor_group list. You only need to do a filter or find to do this.

```genetorch.simulator.false_positive(readfile,gene,[0.05,0.04,0.03,0.02,0.01])```

Give the estimated value of the false positive cases using the gene you chose to compare with others. This step will cost about 5 mins for each value.


## What is new in version 1.1.1

### 1. You can get rid of some annoying background variation by just removing them when 'filter'
```genetorch.finder.filter(taglist,lengthlimit,threshold,rid = ['ttn-1','cla-1'])```
#### by setting ```rid = []``` , the genes in the list will no longer exist in your results.

### 2. You can use the intersection function to find multiple candidate suppressors.
```genetorch.finder.plot(result,size,c='Reds',sub_c='RdPu',intersection_factor = 0.03)```
#### c = color map used in the major graph
#### sub_c = color map used in other graphs
#### intersection_num = 0.03
###### For example: If you think gene A is a suppressor, it means that in samples with other mutations,there is less opportunity for them to carry a gene A mutation. intersection_num means the least percent (default 3%) of the sample number you allow to carry more than one suppressor.
###### You can now find the intersection by click on the bubble of the gene on the graph

### Other changes
```genetorch.finder.intersection(result,samplelist)```
###### You can use this function to find genes not in these samples. Return a Dataframe
```genetorch.finder.gsamplename(result,gene)```
###### You can get a list of your sample carry this gene.
```genetorch.reader.get_impact(taglist)```
###### Remove the lines in your data which do not have an amino_acid change

## 1.0.3 update:
##### bug FIXING........
##### filter is more accurate and faster
##### solve the problem that file cannot be decoded by 'utf-8' in MacOS
```pip install genetorch```
### PLEASE INSTALL: pandas, matplotlib, seaborn BEFORE USE
This package can help find candidate genes in high throughput mutagenesis and suppressor screening experiments
without mapping.
Please call variants via freebayes and annotate with snpeff before using this package.
Any advise is welcomed, please contact
 ### [e-mail]:guozhengyang980525@yahoo.co.jp


```genetorch.reader```

```genetorch.reader.readfile(filepath)```

   multiple renamed vcf files must be included in the filepath.
  ```
filepath
|---1.vcf
|---2.vcf
|---3.vcf
|---4.vcf
|---5.vcf
|---6.vcf
```

```genetorch.reader.getfile(filepath,filename)```

   multiple renamed folders must be included in the filepath, and a vcf file with filename must be included
   in the folders, the name of the folder must be splited with '_' to divide the folder name into strain name and
   WGS order name:
examples:
```
filepath
|---cas113_20221011jxskaosdosh---filename.vcf
|---cas114_20221011jxskaosdosh---filename.vcf
|---cas115_20221011jxskaosdosh---filename.vcf
|---cas116_20221011jxskaosdosh---filename.vcf
|---cas117_20221011jxskaosdosh---filename.vcf
|---cas118_20221011jxskaosdosh---filename.vcf
```


a temp folder will be automatically created in the filepath including renamed vcf files:
  ```
filepath
|----temp
|     |---cas113.vcf
|     |---cas114.vcf
|     |---cas115.vcf
|     |---cas116.vcf
|     |---cas117.vcf
|     |---cas118.vcg
|---cas113_20221011jxskaosdosh---filename.vcf
|---cas114_20221011jxskaosdosh---filename.vcf
|---cas115_20221011jxskaosdosh---filename.vcf
|---cas116_20221011jxskaosdosh---filename.vcf
|---cas117_20221011jxskaosdosh---filename.vcf
|---cas118_20221011jxskaosdosh---filename.vcf
```



```a = genetorch.readfile()```
         ```a = genetorch.getfile()```

   a.taglist : a list of Dataframes which included columns: 'gene', 'ID', 'type', 'base', 'protein'，'tag'
   column 'tag' will be filled with strain name

example:
   a.taglist[1]:

|   gene |    ID      |   type   | base | protein |     tag |
| ---- | -------- |  ------ |  -------  |  ----  | -----  |
| ttn-1 | WBGenexxxx  | missense | C<G | Asp666Asn | cas113 |
| cla-1 | WBGenexxxx  | missense | C<G | Asp223Asn | cas114 |


```genetorch.finder```

```genetorch.finder.find(readfile())```

   return a pandas Dataframe item.
example:

| sample | sample size  | gene | variation | variation_number |
|---------- |--------------| ----- | ---- | ----- |
| cas113,cas114 | 2            | ttn-1 | Asp666Asn,1 Glu374Gln,2 | 2 |

   ttn-1 mutation is found in 2 input files, there are three amino acid variation events, 2 * Glu374Gln, 1 * Asp666Asn.
   There are two kinds of variations in this gene.

```genetorch.finder.filter(readfile(),lengthlimit = 0.6,rid = ['ttn-1','cla-1'])```

   return a pandas Dataframe item same as find().

provide filtered data:

lengthlimit:  if  n(amino acid variation) > lengthlimit * n(input files), exclude the variation
 
if the gene has no variation left after filtered, the gene will be deleted from the result Dataframe


```genetorch.finder.plot(result,size,intersection_num)```

   result: result Dataframe from find() or filter()
   size: use the first n lines to show the plot.  

    x: gene name
    y:variation number
    size: size
    color: variation number/size
    automatically show annotation of each bubble


```genetorch.stocker```

   Module for data upload to LIM database

  ``` genetorch.stocker.stock(path,filename,outpath)```

   ## path & filename
   multiple renamed folders must be included in the filepath, and a vcf file with filename must be included
   in the folders, the name of the folder must be splited with '_' to divide the folder name into strain name and
   WGS order name:
examples:
```
filepath
|---cas113_20221011jxskaosdosh---filename.vcf
|---cas114_20221011jxskaosdosh---filename.vcf
|---cas115_20221011jxskaosdosh---filename.vcf
|---cas116_20221011jxskaosdosh---filename.vcf
|---cas117_20221011jxskaosdosh---filename.vcf
|---cas118_20221011jxskaosdosh---filename.vcf
```

a temp folder will be automatically created in the filepath including renamed vcf files:
  ```
filepath
|----temp
|     |---cas113.vcf
|     |---cas114.vcf
|     |---cas115.vcf
|     |---cas116.vcf
|     |---cas117.vcf
|     |---cas118.vcg
|---cas113_20221011jxskaosdosh---filename.vcf
|---cas114_20221011jxskaosdosh---filename.vcf
|---cas115_20221011jxskaosdosh---filename.vcf
|---cas116_20221011jxskaosdosh---filename.vcf
|---cas117_20221011jxskaosdosh---filename.vcf
|---cas118_20221011jxskaosdosh---filename.vcf
```


### outpath:
   outpath is an empty folder.
   output csv files in this folder which can be uploaded to LIM database.

# Examples

```
import genetorch as gt
a = gt.reader.readfile(path) # read vcf files from path
b = gt.reader.getfile(path,filename) # read vcf files with the filename in the individual folders in the path 
                                             # a has the same properties as b
gt.reader.get_impact(a) # get impact results from a
d = gt.finder.find(a) # get the result table
e = gt.finder.filter(a,lengthlimit = 0.6,threshold = 1, rid = ['ttn-1','cla-1']) # get filtered result table
gt.finder.plot(e,1000,1) # draw a picture using the first 1000 lines in the result table,
# allowing 1 sample have more than one suppressor gene
gt.stocker.stockfile(path,filename,outpath) #get stock file

```
