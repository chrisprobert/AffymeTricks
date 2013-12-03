"""

functions to parse data from a GEO RNA-seq FPKM file

input should be a GEO RNA-seq FPKM file

"""

# Driver method to take a list of affy files and output them
# with standard gene names and normalized count values

def parseFiles(fileNames) :
  fileValues = []

  for fileName in fileNames :
    thisFileValues = parseRNASeq(open(fileName, "rU"))
    fileValues.append(thisFileValues)
  
  return fileValues



# function to parse an affy file
# arguments: handle for the affy file, gene ID dictionary
# returns: dicitonary of geneIDs and counts

def parseRNASeq(fileHandle) :
  
  ret = {}
  
  for line in fileHandle:

    l = line.split()
    if len(l) < 10 : continue    
    
    geneID = l[3]
    FPKM = l[9]
    ret[geneID] = FPKM
  return ret

