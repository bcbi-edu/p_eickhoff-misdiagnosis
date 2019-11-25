'''
Created on Jul 9, 2019

'''

from basicversion import getParents
from quickumls import QuickUMLS
import re
import networkx as nx
import matplotlib as mpl
mpl.use('tkagg')
import matplotlib.pyplot as plt
#from matplotlib import cm

matcher = QuickUMLS(quickumls_fp='umls', overlapping_criteria="length", accepted_semtypes={"T047", "T191"})
generalFile = 'diseases+count.txt'
zeroFile = 'zero-zero.txt'
oneFile = 'one-one-test1.txt'
titleFile = 'titles.txt'
graphStatsFile = 'graph-degrees-weighted.txt'

file = open(titleFile, 'r')

disFile = open(generalFile, 'a')

zeroTozero = open(zeroFile, 'a')

oneToone = open(oneFile, 'a')

graphStats = open(graphStatsFile, 'a')

oneDict = {}
oneNormDict = {}
nonOneDis = {}
countsDict = {}
cuiToDis = {}
parCUIs = []
chdToPar = {}

G = nx.DiGraph()

counter = 0
onesCounter = 0

"""
given a title, gets the matches and stores them as two lists, before and after, 
and returns these two lists as a tuple
return (before, after), where before is the list of matches before the pattern 
and after is the list of matches after the pattern
"""
def handleMatches(title):
    m = matcher.match(title, best_match=True, ignore_syntax=False)
    terms = []
    mStart = re.search('(?:misdiagnosed as | masquerading as)', title).start()
    for elem in m:
        start = 0
        end = 0
        sim = 0
        cui = 0
        term = ""
    
        for posMatch in elem:
            if posMatch.get("end") - posMatch.get("start") > end - start:
                start = posMatch.get("start")
                end = posMatch.get("end")
                sim = posMatch.get("similarity")
                cui = posMatch.get("cui")
                term = posMatch.get("term").lower()
            elif posMatch.get("end") - posMatch.get("start") == end - start:
                if posMatch.get("similarity") > sim:
                    start = posMatch.get("start")
                    end = posMatch.get("end")
                    sim = posMatch.get("similarity")
                    cui = posMatch.get("cui")
                    term = posMatch.get("term").lower()
                elif posMatch.get("similarity") == sim:
                    if posMatch.get("cui") in cuiToDis:
                        start = posMatch.get("start")
                        end = posMatch.get("end")
                        sim = posMatch.get("similarity")
                        cui = posMatch.get("cui") 
                        term = posMatch.get("term").lower()
                        break
                    elif posMatch.get("cui") < cui:
                        start = posMatch.get("start")
                        end = posMatch.get("end")
                        sim = posMatch.get("similarity")
                        cui = posMatch.get("cui") 
                        term = posMatch.get("term").lower()
          
        newDict = {"start": start, "end": end,  "disease": term, "cui": cui}
        terms.append(newDict)
    before = list(filter(lambda x: x.get("start") < mStart, terms))
    after = list(filter(lambda x: x.get("start") > mStart, terms))
    if len(before) == 1 and len(after) == 1:
        cuiToDis[before[0].get("cui")] = before[0].get("disease")
        cuiToDis[after[0].get("cui")] = after[0].get("disease")
    return (before, after)

"""
updates the counts of pairs of diseases
"""
def count_one_one(befDis, aftDis):
    if befDis in oneDict:
        if aftDis in oneDict[befDis]:
            oneDict[befDis][aftDis] += 1
        else:
            oneDict[befDis][aftDis] = 1
    else:
        oneDict[befDis] = {aftDis: 1}

"""
updates the counts of a disease being misdiagnosed    
"""        
def getNormal(befDis):
    if befDis in oneNormDict:
        oneNormDict[befDis] = oneNormDict[befDis] + 1
    else:
        oneNormDict[befDis] = 1
        

"""
    first checks if synonymous relationship with existing CUI
    if only one disease has SY, sends it to check parents of other CUI
    if both have SY, uses SY CUI
    if neither, check parents of both
"""           
def handle_relations(befDis, aftDis):
    simDis1Par, simDis1Chd, simDis1Syn = getParents.rel_cui_list(befDis)
    simDis2Par, simDis2Chd, simDis2Syn = getParents.rel_cui_list(aftDis)
    
    dis1 = befDis
    dis2 = aftDis
    
    for simBef in simDis1Syn:
        if simBef != befDis and cuiToDis.get(simBef) != None:
            dis1 = simBef
    for simAft in simDis2Syn:
        if simAft != aftDis and cuiToDis.get(simAft) != None:
            dis2 = simAft
    if dis1 != befDis and dis2 != aftDis:
        getNormal(dis1)
        count_one_one(dis1, dis2)
    elif dis1 != befDis:
        check_parents([dis1], simDis2Par, dis1, dis2, [dis1], simDis2Chd)
    elif dis2 != aftDis:
        check_parents(simDis1Par, [dis2], dis1, dis2, simDis1Chd, [dis2])
    else:
        check_parents(simDis1Par, simDis2Par, befDis, aftDis, simDis1Chd, simDis2Chd)
        
    
"""
    checks if parent CUIs are in dictionary and adds counts to parents
    if children CUIs are in dictionary, adds children counts to it and 
    removes children
    
    parList1 and chdList1 correspond to befDis
    parList2 and chdList2 correspond to aftDis
"""
def check_parents(parList1, parList2, befDis, aftDis, chdList1, chdList2):
    dis1 = befDis
    dis2 = aftDis
    if befDis in chdToPar:
        dis1 = chdToPar[befDis]
        cuiToDis.pop(befDis, None)
    else:
        for simBef in parList1:
            if simBef != befDis and cuiToDis.get(simBef) != None:
                dis1 = simBef
                cuiToDis.pop(befDis, None)
    if aftDis in chdToPar:
        dis2 = chdToPar[aftDis]
        cuiToDis.pop(aftDis, None)
    else:
        for simAft in parList2:
            if simAft != aftDis and cuiToDis.get(simAft) != None:
                dis2 = simAft
                cuiToDis.pop(aftDis, None)
    handleChildren(dis1, chdList1)
    handleChildren(dis2, chdList2)
    getNormal(dis1)
    count_one_one(dis1, dis2)
    
    
"""
    checks the children in chdList1 of befDis to see if it already exists and 
    adds to parent if yes
    
    for when befDis has a child but that child does not have befDis as a parent
"""                    
def handleChildren(befDis, chdList1):
    for chd in chdList1:
        if chd != befDis and chd in cuiToDis:
            chdToPar[chd] = befDis
            del cuiToDis[chd]
            if chd in oneNormDict:
                if oneNormDict.get(befDis):
                    oneNormDict[befDis] += oneNormDict.pop(chd)
                else:
                    oneNormDict[befDis] = oneNormDict.pop(chd)
                x = oneDict.get(befDis, {})
                x.update(oneDict.pop(chd))
                oneDict[befDis] =  x
            for v in oneDict.values():
                if chd in v:
                    if befDis in v:
                        v[befDis] += v.pop(chd)
                    else:
                        v[befDis] = v.pop(chd)
        
                 
            
 
"""
    writes to the zero-zero file if QuickUMLS extracts zero diseases on one 
    side of the pattern
    befList - list of disease names before the misdiagnosed/masquerading 
    aftList - list of diases names after the misdiagnosed/masquerading
    title - the title of the article
"""       
def handleZeros(befList, aftList, title):
    zeroTozero.write(title + "\n     Before:\n")
    for dis1 in befList:
        zeroTozero.write("        -" + dis1.get("disease") + "\n")
    zeroTozero.write("     After:\n")
    for dis2 in aftList:
        zeroTozero.write("        -" + dis2.get("disease") + "\n")
    zeroTozero.write("----------------------------------\n")
    
"""
    adds counts for titles with irregular match patterns
"""
def addingCounts(dis1, dis2, disDict):
    pair = (dis1, dis2)
    if pair in disDict:
        disDict[pair] = disDict[pair] + 1
    else:
        disDict[pair] = 1

"""
    writes the title and the pairs to the general file while adding counts
"""    
def handle_all(bef, aft, title, disDict):
    disFile.write(title + "\n")
    for dis1 in bef:
        for dis2 in aft:
            disFile.write("     " + dis1.get("disease") + " --> " + dis2.get("disease") + "\n")
            addingCounts(dis1.get("disease"), dis2.get("disease"), disDict)
    disFile.write("-------------------------------\n")
    
"""
    counts the number of each type of pattern 
    (aka 0 diseases misdiagnosed as 1 disease, 1 disease misdiagnosed as 1 disease, etc)
"""
def getCounts(bef, aft, disTodisCounts):
    pair = (len(bef), len(aft))
    if pair in disTodisCounts:
        disTodisCounts[pair] = disTodisCounts[pair] + 1
     
"""
    normalizes the one-one pair counts and writes it to the file
    adds nodes to graph
"""       
def normalize():
    global counter
    global onesCounter
    for befDis in sorted(oneDict):
        for aftDis in oneDict[befDis]:
            if befDis in cuiToDis and aftDis in cuiToDis:
                pairTerms = (cuiToDis[befDis], cuiToDis[aftDis])
                ratio = (oneDict[befDis][aftDis]/oneNormDict[befDis])
                G.add_node(befDis)
                G.add_node(aftDis)
                G.add_edge(befDis, aftDis, weight = ratio)
                counter = counter + 1
                onesCounter = onesCounter + 1
                oneToone.write(" --> ".join((befDis, aftDis)) + "\n" + 
                               " --> ".join(pairTerms) + 
                               "\n    normalized count: " + ratio + 
                               "\n      total of first: " + 
                               str(oneNormDict[befDis]) +  "\n\n")

"""
    writes the counts of each disease pair from all titles in the general file
"""     
def writeCounts():
    global counter
    for pair in sorted(nonOneDis):
        disFile.write(pair[0] + " --> " + pair[1] + " " + str(nonOneDis[pair]) + "\n")
        counter = counter + nonOneDis[pair]

"""
   writes the counts of the pattern types (from getCounts) to the general file 
"""
def writeTotNum():
    global counter
    disFile.write(str(counter) + "\n")
    for pair in countsDict:
        disFile.write(str(pair[0][0]) + "diseases --> " + str(pair[1][0]) + "\n")

"""
    goes through each title, getting the matches, filters out the ones without 
    exactly 1 cui on either side to get counts for
"""
def runMatches():
    for line in file:
        matches = handleMatches(line.lower())
        getCounts(matches[0], matches[1], countsDict)
        if len(matches[0]) == 1 and len(matches[1]) == 1:
            count_one_one(matches[0][0].get("cui"), matches[1][0].get("cui"))
            handle_relations(matches[0][0].get("cui"), matches[1][0].get("cui"))
        #elif len(matches[0]) == 0 or len(matches[1]) == 0:
         #   handleZeros(matches[0], matches[1], line)
        #handle_all(matches[0], matches[1], line.lower(), nonOneDis)

"""
    to get a graph of related nodes 
    (aka all the nodes connected to one disease and all the nodes connected to 
    its connections, and so on)
    
    for visualization
"""     
def getGroups(pair, graph):
    #group = {pair: oneDict.get(pair)}
    firstDis = filter(lambda x: x[0] == pair[0] or x[1] == pair[0] and x != pair, oneDict)
    secDis = filter(lambda x: x[0] == pair[1] or x[1] == pair[1] and x != pair, oneDict)
    if len(firstDis) > 0:
        for item in firstDis:
            if item not in firstDis:
                graph.add_node(item[0])
                graph.add_node(item[1])
                graph.add_edge(item[0], item[1], oneDict[item]/oneNormDict[item[0]])
                getGroups(item)
    if len(secDis) > 0:
        for item in secDis:
            if item not in secDis:
                graph.add_node(item[0])
                graph.add_node(item[1])
                graph.add_edge(item[0], item[1], oneDict[item]/oneNormDict[item[0]])
                getGroups(item)

                


runMatches()
normalize()
#writeCounts()
print("total pairs: " + str(counter) + ", ones pairs: " + str(onesCounter))


""" graph visualization and degree analysis """

pos=nx.spring_layout(G)
fig = plt.figure(3,figsize=(12,12))
edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
d = nx.out_degree_centrality(G)
inDeg = nx.in_degree_centrality(G)

numNodes = 0
totOutDeg = 0
totInDeg = 0

for v in d.values():
    numNodes += 1
    totOutDeg += v
    
for v in inDeg.values():
    totInDeg += v

graphStats.write("avg out-degree-centrality: " + str(totOutDeg/numNodes) + 
                 "\n\navg in-degree-centrality: " + str(totInDeg/numNodes) + "\n\n")

totBet = nx.betweenness_centrality(G, k=numNodes)
bet = 0
for n in totBet.values():
    bet += n
#graphStats.write("avg betweeness-centrality: " + str(bet / numNodes) + "\n\n")

totIn = 0
totOut = 0

inDeg = {}
outDeg = {}

for n in list(G.nodes):
    #nIn = G.in_degree(n, weight = "weight")
    inDeg[n] = G.in_degree(n)
    totIn += G.in_degree(n)
    
    #nOut = G.out_degree(n, weight = "weight")
    outDeg[n] = G.out_degree(n)
    totOut += G.out_degree(n)
    #graphStats.write(n +  "\nout-degree: " + str(nOut) + "\nin-degree: " + str(nIn) + "\nbetweeness centrality: " + str(totBet.get(n)) + "\n\n")

graphStats.write("avg in-degree: " + str(totIn/numNodes) + 
                 "\navg out-degree: " + str(totOut/numNodes) + "\n\n")

nx.draw(G,pos, edge_color=weights, width=1, node_color = 'b', with_labels=False, 
        node_size=[v * 10000 + 10 for v in d.values()], font_size=8, 
        edge_cmap=plt.cm.get_cmap('Purples'))

fig.set_facecolor("skyblue")
plt.show()

file.close()
disFile.close()
zeroTozero.close()
oneToone.close()
graphStats.close()
