import re
import argparse

ap = argparse.ArgumentParser()
ap.add_argument ("gff_pc", type=str )
ap.add_argument("gff_lncRNA", type=str)
args = ap.parse_args()

fn1 = args.gff_pc
fn2 = args.gff_lncRNA
pc_information = list()

def gffparse(line):
  if re.match("#", line):
    return
  line = line.rstrip()
  fields = line.split("\t")
  #print (fields[2])
  if len(fields) != 9:
    return
  the_feature = {
     'scaffold' : fields[0],
     'source'   : fields[1],
     'type'     : fields[2],
     'start'    : fields[3],
     'end'      : fields[4], 
     'score'    : fields[5],
     'strand'   : fields[6],
     'phase'    : fields[7],
     #'id'       : re.match('ID=([^;]*)', fields[8]).group(),
     
  }
  the_feature['attrib'] = dict(item.split("=") for item in fields[8].split(";"))  
  return the_feature

blacklist_transcript = list()
for line1 in (open(fn1)):
    if re.match("#", line1):
        continue
    field = line1.strip().split("\t")
    if field[2] == "gene":
        pc_information.append((field[0],int(field[3]) - 1000 ,int (field[4]) + 1000 ,field[6]))
    #print (pc_information[0][0],pc_information[0][1],pc_information[0][2])

for line in (open(fn2)):
  #print (line)
    my_gff = gffparse(line)
    if re.match("#", line):
        continue
    
    if (my_gff["type"] == 'transcript'):
        for pc in pc_information:
            if (my_gff["scaffold"] == pc[0] and int (my_gff["start"]) >=  int(pc[1]) and int (my_gff["end"]) <= int (pc[2]) and my_gff["strand"] == pc[3]) or (my_gff["scaffold"] == pc[0] and  int (my_gff["start"]) <=  int(pc[1]) and int (my_gff["end"]) >= int (pc[1]) and my_gff["strand"] == pc[3]) or (my_gff["scaffold"] == pc[0] and  int (my_gff["start"]) <=  int(pc[2]) and int (my_gff["end"]) >= int (pc[2]) and my_gff["strand"] == pc[3]):
                id = my_gff["attrib"].get("ID")
                blacklist_transcript.append(id)
                #print (blacklist_transcript)
               #print (my_gff["scaffold"], pc[0], my_gff["start"],my_gff["end"], pc[1],pc[2],my_gff["strand"],pc[3])
#print (len(blacklist_transcript))
for line in (open(fn2)):
    my_gff = gffparse(line)
    if re.match("#", line):
        continue
    if (my_gff["type"] == 'transcript'):
        if (my_gff["attrib"].get("ID") not in blacklist_transcript):
            print (line, end="")
            pass
    elif (my_gff["attrib"].get("Parent") not in blacklist_transcript ):
        print (line, end="")
        pass
    

