import requests, sys
import json
import wget
import urllib.request
 
server = "https://parasite.wormbase.org"
ext = "/rest-15/info/genomes/taxonomy/Platyhelminthes"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json", "Accept" : ""})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
#decoded = r.json()
#print(repr(decoded))

data = json.loads(r.text)
for i, line in enumerate(data):
    resultObject = data[i]['name']
    x = resultObject.split("_")
    url = f'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/{x[0]}_{x[1]}/{x[2].upper()}/{x[0]}_{x[1]}.{x[2].upper()}.WBPS15.protein.fa.gz'
    #urllib.request.urlretrieve(url, '/lustre/scratch118/infgen/team133/as57/annotation/')
    filename = wget.download(url)
    print(url)
