import parseAffyData
import parseRNASeqData
import operator

"""

Main program of the modular AffymeTricks package.
Designed to allow modular expansion by providing object
based access and integration of different data sources

"""


# main method, can easily be re-configured
# or invoked via command line statements
def main() :

  # Sample data variables :
  gsm1 = ["GSM1.tsv"]
  RNA1 = ["RNA1.fpkm"]
  
  gsmresults = parseAffyData.parseFiles(gsm1, "map.txt")
  rnaseqrest = parseRNASeqData.parseFiles(RNA1)

  affyCounts = []
  rnaCounts  = []
  for k in gsmresults : affyCounts.append(float(0))
  for k in rnaseqrest : rnaCounts.append(float(0))

  for i in range(0, len(gsmresults)) :
    todel = []
    for k in gsmresults[i] :
      try:
        affyCounts[i] += float(gsmresults[i][k])
      except ValueError :
        todel.append(k)
        continue
    for k in todel :
      del gsmresults[i][k]
    for k in gsmresults[i] :
      gsmresults[i][k] = float(gsmresults[i][k]) * 1000000 / affyCounts[i]


  for i in range(0, len(rnaseqrest)) :
    todel = []
    for k in rnaseqrest[i] :
      try:  
        rnaCounts[i] += float(rnaseqrest[i][k])
      except ValueError :
        todel.append(k)
        continue
    for k in todel :
      del rnaseqrest[i][k]
    for k in rnaseqrest[i] :
      rnaseqrest[i][k] = float(rnaseqrest[i][k]) * 1000000 / rnaCounts[i]

  aggreg = {}
  for res in gsmresults :
    for k in res :
      if k not in aggreg : aggreg[k] = float(res[k])
      else: aggreg[k] += float(res[k])
  
  for res in rnaseqrest :
    for k in res :
      if k not in aggreg : aggreg[k] = float(res[k])
      else : aggreg[k] += float(res[k])

  sortedAggreg = sorted(aggreg.iteritems(), key=operator.itemgetter(1))[::-1]
  sampleNames = gsm1 + RNA1
  
  for g in range(100, 110) :
    gene = sortedAggreg[g][0]
    print gene
  
  header = "Gene"
  for n in sampleNames : header += "\t" + n
  print "Expected values of Fragments Per Million (FPM) for most expressed genes accross all samples:"
  print header
  for g in range(0, 100) :
    gene = sortedAggreg[g][0]
    line = gene
    for n in range(0,len(gsm1)) :
      count = float(0)
      if gene in gsmresults[n] : count = gsmresults[n][gene]
      line += "\t" + "{0:.2f}".format(count)
    for n in range(0,len(RNA1)) :
      count = float(0)
      if gene in rnaseqrest[n] : count = rnaseqrest[n][gene]
      line += "\t" + "{0:.2f}".format(count)
    print line
      

if __name__ == "__main__" : main()   
