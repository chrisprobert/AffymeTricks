import numpy as np
import os

# generic function to dump dictionary to tsv file
def dumpDict(dictd, filename) :
  f = open(os.path.join("Sample_DE_Data", filename), "w")
  for k in dictd :
    f.write(k + "\t" + "{0:.4f}".format(dictd[k]) + "\n")
  f.close()

geneNames = []
DEGenes = []
inf = open("geneNames.txt", "rU")
for line in inf : geneNames.append(line.strip())
inf = open("diffExpressList.txt", "rU")
for line in inf : DEGenes.append(line.strip())


# create set1 - has genes with FCs all from same distribution
for i in range(0, 20) :
  outname = "set1_sample" + str(i) + ".tsv"
  dictd = {}
  for g in geneNames :
    dictd[g] = np.random.normal(0, .05)
  for g in DEGenes :
    dictd[g] = np.random.normal(0, .05)
  dumpDict(dictd, outname)

# create set2 - has positive distribution on DE genes
for i in range(0, 20) :
  outname = "set2_sample" + str(i) + ".tsv"
  dictd = {}
  for g in geneNames :
    dictd[g] = np.random.normal(0, .05)
  for g in DEGenes :
    dictd[g] = np.abs(np.random.normal(0.05, .05))
  dumpDict(dictd, outname)

# create set3 - has stronger positive distribution
for i in range(0, 20) :
  outname = "set3_sample" + str(i) + ".tsv"
  dictd = {}
  for g in geneNames :
    dictd[g] = np.random.normal(0, .05)
  for g in DEGenes :
    dictd[g] = np.abs(np.random.normal(0.1, .1))
  dumpDict(dictd, outname)
 
