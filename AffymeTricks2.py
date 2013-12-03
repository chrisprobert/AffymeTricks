import parseAffyData
import parseRNASeqData
import operator
import numpy as np

"""

Part of the AffymeTricks package designed for differential
expression analysis on data collected from different
sequencing platforms and/or normalized in different ways.

Takes lists of input files that are DE profiles of a set
of genes between two conditions in samples. The program
requires at least two such lists (e.g. 2 or more sets of
DE data) in order to perform analysis between samples.

"""


# main method, can easily be re-configured
# or invoked via command line statements
def main() :

  SET_SIZE = 20

  # Sample data variables :
  set1f = []
  set2f = []
  set3f = []

  for i in range(0,SET_SIZE) :
    set1f.append("./Sample_DE_Data/set1_sample%s.tsv" % str(i))
    set2f.append("./Sample_DE_Data/set2_sample%s.tsv" % str(i))
    set3f.append("./Sample_DE_Data/set3_sample%s.tsv" % str(i))

  set1 = parseFiles(set1f) 
  set2 = parseFiles(set2f)
  set3 = parseFiles(set3f)

  totalScores = [] # keep the score dictionaries for all samples
  agregScores = {} # gene-by-gene scores
 
  # for testing purposes, assume all sets have the same genes
  # could implement this with upstream filtering / joining if we wanted...
  for s in [set1, set2] :
    magns = {}  # store the magnitudes of FCs
    signs = {}  # store the signs of FCs
    for k in s[0] :
      magns[k] = float(0)
      signs[k] = int(0)
    
    # for each gene, test for correlation accross samples  
    for k in magns :
      for si in s :
        magns[k] += si[k]
        if si[k] < 0 : signs[k] -= 1
        else : signs[k] += 1
      netsign = np.sign(signs[k])
      frac = float(signs[k])
      frac = float(frac/SET_SIZE)
      signs[k] = np.power(frac, 3) # use k = 3, could allow as parameter
      if np.sign(magns[k]) != netsign : signs[k] = 0 # set to 0 if magn and sign disagree
      else : signs[k] = signs[k] * float(magns[k])

    # store these results in totalScores list
    totalScores.append(signs)
    
  # now, go through and find genes with highest scores accross all samples
  for k in totalScores[0] :
    agregScores[k] = float(0)
    for s in totalScores :
      agregScores[k] += float(s[k])
      
  sortedAggreg = sorted(agregScores.iteritems(), key=operator.itemgetter(1))[::-1]
  
  header = "Gene\tScore"
  print "Correlation scores of top 10 aggregate genes"
  print header
  for g in range(0, 10) :
    line = sortedAggreg[g][0] + "\t" + "{0:.5f}".format(sortedAggreg[g][1])
    print line

def parseFiles( files ) :
  ret = []
  for filename in files :
    d = {}
    f = open(filename, "rU")
    for line in f :
      l = line.split()
      d[l[0]] = float(l[1].strip())
    ret.append(d)
  return ret
      

if __name__ == "__main__" : main()
