'''
Created on Jun 22, 2019

'''

###Takes a list of CUIs that result from /search?string=PT,IN,MIN,PIN&inputType=tty&sabs=RXNORM,USPMG.
###For each of theses CUIs, output the atoms
###If there is a PIN from RxNorm in the same CUI as a USPMG atom, find the RxNorm associated single ingredient form and
###remove it from the final list
###Then go through again and remove USP atoms that are in the same CUI as an RxNorm atom

# adapted from uts-rest-api samples
# Authentication.py can be found uts-rest-api

from Authentication import Authentication

import requests
import simplejson
import argparse

parser = argparse.ArgumentParser(description='process user given parameters')
parser.add_argument("-k", "--apikey", required = True, dest = "apikey")
parser.add_argument("-v", "--version", required =  False, dest="version", default = "current", help = "enter version example-2015AA")
parser.add_argument("-o", "--outputfile", required = False, dest = "outputfile", help = "enter a name for your output file")
parser.add_argument("-s", "--sabs", required = False, dest="sabs",help = "enter a comma-separated list of vocabularies, like MSH, SNOMEDCT_US, or RXNORM")
parser.add_argument("-t", "--ttys", required = False, dest="ttys",help = "enter a comma-separated list of term types, like PT,SY,IN")
parser.add_argument("-i", "--inputfile", required = False, dest = "inputfile", help = "enter a name for your input file")


args = parser.parse_args(["--apikey", "<key>"])
apikey = args.apikey
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

 
def uts_get(path,query, full):
    if full:
        r = requests.get(path, params=query)
        r.encoding = 'utf-8'
        return simplejson.loads(r.text)
    else:
        r = requests.get(base_uri + path, params=query)
        r.encoding = 'utf-8'
        try:
            return simplejson.loads(r.text)
        except:
            print("page not found error")
  
  
def retrieveConceptAtoms(cui):
  
    path = "/rest/content/current/CUI/"+cui+"/atoms/preferred"
    query = {"ticket":AuthClient.getst(tgt)}
    if sabs:
        query["sabs"] = sabs
    if ttys:
        query["ttys"] = ttys
     
    results = uts_get(path,query, False)
    return results
    

def get_link(line):
    cui = line.strip()
    json = retrieveConceptAtoms(cui)
    return json["result"]["relations"]
        
    
def get_rel_list(par_uri):
    query = {"ticket":AuthClient.getst(tgt)}
    if sabs:
        query["sabs"] = sabs
    if ttys:
        query["ttys"] = ttys
        
    results = uts_get(par_uri, query, True)

    return results["result"]


def rel_cui(label, aui_link):
    
    query = {"ticket":AuthClient.getst(tgt)}
    results = uts_get(aui_link, query, True)
    try:             
        concept = results["result"]["concept"]
    except:
        prefAtom = results["result"]["defaultPreferredAtom"]
        query = {"ticket":AuthClient.getst(tgt)}
        finalResults = uts_get(prefAtom, query, True)
        concept = finalResults["result"]["concept"]
    
    con_len = len(concept)
    simCUI = concept[con_len - 8:]
    
    return simCUI    
    

def rel_cui_list(cui):
    par_list = set([cui])
    chd_list = set([cui])
    sim_list = set([cui])
    
    try:
        par = get_link(cui)
        if par != "NONE":   
            atomList = get_rel_list(par)
            for atom in atomList:
                label = atom.get("relationLabel")
                aui_link = atom["relatedId"]
                relative = rel_cui(label, aui_link)
                if label == "CHD" or label == "RN":
                    par_list.add(relative)
                elif label == "PAR" or label == "RB":
                    chd_list.add(relative)
                elif label == "SY" or label == "RL":
                    sim_list.add(relative)
    except:
        pass
    
    return par_list, chd_list, sim_list

        
def main():
    return 1;
    
if __name__ == '__main__':
    main()
