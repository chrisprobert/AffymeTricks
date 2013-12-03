"""

functions to parse data from an affay array

input should be of the form "Affy_ID intensity"

"""

# Driver method to take a list of affy files and output them
# with standard gene names and normalized count values

def parseFiles(fileNames, geneDictFile) :
  geneDict = getGeneDict(geneDictFile)
  fileValues = []

  for fileName in fileNames :
    thisFileValues = parseAffy(open(fileName, "rU"), geneDict)
    fileValues.append(thisFileValues)
  
  return fileValues


# function to parse an affy file
# arguments: handle for the affy file, gene ID dictionary
# returns: dicitonary of geneIDs and counts

def parseAffy(fileHandle, geneDict) :
  
  ret = {}
  
  for line in fileHandle:

    l = line.split()
    if len(l) != 2 : continue    
    if l[0] not in geneDict : continue
   
    geneName = geneDict[l[0]] 
    if geneName not in ret: ret[geneName] = 0
    ret[geneName] += float(l[1].strip())

  return ret



# helper method to get a dict of entries in the
# gene name convsion table
def getGeneDict(geneDictFile) :
  retDict = {}
  f = open(geneDictFile, "rU")
  for line in f :
    l = line.split()
    if len(l) != 3 : continue
    retDict[l[0]] = l[2]
  
  return retDict
