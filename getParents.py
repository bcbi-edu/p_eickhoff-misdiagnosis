'''
Created on Jun 22, 2019

@author: cindyli
'''

###Takes a list of CUIs that result from /search?string=PT,IN,MIN,PIN&inputType=tty&sabs=RXNORM,USPMG.
###For each of theses CUIs, output the atoms
###If there is a PIN from RxNorm in the same CUI as a USPMG atom, find the RxNorm associated single ingredient form and
###remove it from the final list
###Then go through again and remove USP atoms that are in the same CUI as an RxNorm atom

from basicversion.authentication import Authentication

import requests
import simplejson
import argparse
import collections
from collections import OrderedDict

parser = argparse.ArgumentParser(description='process user given parameters')
parser.add_argument("-k", "--apikey", required = True, dest = "apikey", help = "f0ea53a8-d79a-463f-af1b-24fd82e5b9b3")
parser.add_argument("-v", "--version", required =  False, dest="version", default = "current", help = "enter version example-2015AA")
parser.add_argument("-o", "--outputfile", required = False, dest = "outputfile", help = "enter a name for your output file")
parser.add_argument("-s", "--sabs", required = False, dest="sabs",help = "enter a comma-separated list of vocabularies, like MSH, SNOMEDCT_US, or RXNORM")
parser.add_argument("-t", "--ttys", required = False, dest="ttys",help = "enter a comma-separated list of term types, like PT,SY,IN")
parser.add_argument("-i", "--inputfile", required = False, dest = "inputfile", help = "enter a name for your input file")


args = parser.parse_args(["--apikey", "f0ea53a8-d79a-463f-af1b-24fd82e5b9b3"])
apikey = args.apikey
#version = args.version
#outputfile = args.outputfile
#inputfile = args.inputfile
sabs = args.sabs
ttys = args.ttys
AuthClient = Authentication(apikey) 

###################################
#get TGT for our session
###################################

tgt = AuthClient.gettgt()
base_uri = "https://uts-ws.nlm.nih.gov"
rxnorm = "https://rxnav.nlm.nih.gov"
pageNumber=1
pageCount=1
auis = {}

 
def rxnorm_get(path,query):
    r = requests.get(rxnorm+path, params=query)
    r.encoding = 'utf-8'
    #print(r.url)
    return simplejson.loads(r.text)
 
def uts_get(path,query, full):
    if full:
        r = requests.get(path, params=query)
        r.encoding = 'utf-8'
        #print(r.url)
        return simplejson.loads(r.text)
    else:
        r = requests.get(base_uri + path, params=query)
        r.encoding = 'utf-8'
        #print(r.url)
        try:
            return simplejson.loads(r.text)
        except:
            print("page not found error")
  
def getSingleIngredientForm(rxcui):
    path = "/REST/rxcui/"+rxcui+"/related.json?tty=IN"
    query = {}
    #print(r.url)
    results=rxnorm_get(path,query)
    return results
  
def retrieveConceptAtoms(cui):
   
    path = "/rest/content/current/CUI/"+cui+"/atoms/preferred"
    query = {"ticket":AuthClient.getst(tgt)}
    if sabs:
        query["sabs"] = sabs
    if ttys:
        query["ttys"] = ttys
     

     
    results = uts_get(path,query, False)
    return results
    
#with open(inputfile, 'r') as f:
def getParents(line):
    cui = line.strip()
    json = retrieveConceptAtoms(cui)
    #print(json["result"]["relations"])
    #print(json["result"]["ui"])
    return json["result"]["relations"]
        
    
def getParList(parURI):
    query = {"ticket":AuthClient.getst(tgt)}
    if sabs:
        query["sabs"] = sabs
    if ttys:
        query["ttys"] = ttys
    results = uts_get(parURI, query, True)
    
    #print(results["result"])
    return results["result"]
    
    
    

def parCUIList(cui):
    cuiListPar = [cui]
    cuiListChd = [cui]
    cuiListSim = [cui]
    try:
        #print("kmn")
        par = getParents(cui)
        #print("hello")
        #print(par)
        #print("bye")
        
        if par != "NONE":   
            atomList = getParList(par)
            #print(atomList)
            for atom in atomList:
                rel = atom.get("relationLabel")
                if rel == "CHD" or rel == "RN":
                    auiLink = atom["relatedId"]
                    query = {"ticket":AuthClient.getst(tgt)}
                    results = uts_get(auiLink, query, True)
                    #print(results["result"])
                    try:             
                        concept = results["result"]["concept"]
                    except:
                        prefAtom = results["result"]["defaultPreferredAtom"]
                        query = {"ticket":AuthClient.getst(tgt)}
                        finalResults = uts_get(prefAtom, query, True)
                        concept = finalResults["result"]["concept"]
                    conLen = len(concept)
                    simCUI = concept[conLen - 8:]
                    if simCUI not in cuiListPar:
                        cuiListPar.append(simCUI)
                if rel == "RB" or rel == "PAR":
                    auiLink = atom["relatedId"]
                    query = {"ticket":AuthClient.getst(tgt)}
                    results = uts_get(auiLink, query, True)
                    #print(results["result"])
                    try:             
                        concept = results["result"]["concept"]
                    except:
                        prefAtom = results["result"]["defaultPreferredAtom"]
                        query = {"ticket":AuthClient.getst(tgt)}
                        finalResults = uts_get(prefAtom, query, True)
                        concept = finalResults["result"]["concept"]
                    conLen = len(concept)
                    simCUI = concept[conLen - 8:]
                    if simCUI not in cuiListChd:
                        cuiListChd.append(simCUI)
                if rel == "RL" or rel == "SY":
                    auiLink = atom["relatedId"]
                    query = {"ticket":AuthClient.getst(tgt)}
                    results = uts_get(auiLink, query, True)
                    #print(results["result"])
                    try:             
                        concept = results["result"]["concept"]
                    except:
                        prefAtom = results["result"]["defaultPreferredAtom"]
                        query = {"ticket":AuthClient.getst(tgt)}
                        finalResults = uts_get(prefAtom, query, True)
                        concept = finalResults["result"]["concept"]
                    conLen = len(concept)
                    simCUI = concept[conLen - 8:]
                    if simCUI not in cuiListSim:
                        cuiListSim.append(simCUI)
        return cuiListPar, cuiListChd, cuiListSim  
    except:
        return cuiListPar, cuiListChd, cuiListSim 

        
#C0013153
#C1700185
        
def main():
    #print(simCUIList("C0030581"))
    #x, y, z = parCUIList("C0030581")
    #x, y, z = parCUIList("C0001432")
    x, y, z = parCUIList("C0205642")
    print(x)
    print(y)
    print(z)
    
if __name__ == '__main__':
    main()
