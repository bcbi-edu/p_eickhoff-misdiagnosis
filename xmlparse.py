import os
import xml.etree.cElementTree as ET
import re
import gzip
#import zipfile
#from io import BytesIO
#import xml.dom.minidom as minidom
#from xml.dom.minidom import parse, parseString

#directory = '/Users/cindyli/Documents/CS/Research/workspace-biomed/pubmeddata/basicversion/data'


titles = []


def getTitles(directory):  
    for member in os.listdir(directory):
        if member.endswith(".gz"):
            with gzip.open(directory + "/" + member, 'rb') as f:
                art = f.read()
                data = ET.fromstring(art)
                for child in data.iter('PubmedArticle'):
                    title = child.find('./MedlineCitation/Article/ArticleTitle')
                    if not title is None:
                        if not title.text is None:
                            amatch = re.findall('.+(?:misdiagnosed as | masquerading as).+', title.text)
                            if amatch:                         
                                titles.extend(amatch)
        

         

getTitles('/Users/cindyli/Documents/CS/Research/workspace-biomed/pubmeddata/basicversion/data/pubmed/')

with open('/Users/cindyli/Documents/CS/Research/workspace-biomed/pubmeddata/basicversion/titles.txt', 'a') as file:
    for elt in titles: 
        file.write(elt + "\n")

            
       



