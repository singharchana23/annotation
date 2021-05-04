import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument ("gff_pc", type=str) # protein_coding gff
ap.add_argument ("gff_lncRNA", type=str) # lncRNA 
args = ap.parse_args()

fn1 = args.gff_pc
fn2 = args.gff_lncRNA

def gffparse(line):
    if re.match ("#", line):
        return
    fields = line.rstrip().split("\t")
    the_feature = {
        'scaffold' : fields[0],
        'source'   : fields[1],
        'type'     : fields[2],
        'start'    : fields[3],
        'end'      : fields[4], 
        'score'    : fields[5],
        'strand'   : fields[6],
        'phase'    : fields[7],
        'attribute' : fields[8],
    }
    the_feature['attri'] = dict(item.split("=") for item in fields[8].split(";"))  
    return the_feature
def overlap (coordinate_pc, coordinate_lc):
    sense = False
    antisense = False
    intergenic = False 
    at = ""
    for pc in coordinate_pc:
            if (coordinate_lc[0] == pc[0]):
                if (int (coordinate_lc[1]) >=  int(pc[1]) and int (coordinate_lc[2]) <= int (pc[2]) and coordinate_lc[3] == pc[3]) or (int (coordinate_lc[1]) <=  int(pc[1]) and int (coordinate_lc[2]) >= int (pc[1]) and coordinate_lc[3] == pc[3]) or (int (coordinate_lc[1]) <=  int(pc[2]) and int (coordinate_lc[2]) >= int (pc[2]) and coordinate_lc[3] == pc[3]):
                    sense = True
                    break
                elif (int (coordinate_lc[1]) >=  int(pc[1]) and int (coordinate_lc[2]) <= int (pc[2]) and coordinate_lc[3] is not pc[3]) or (int (coordinate_lc[1]) <=  int(pc[1]) and int (coordinate_lc[2]) >= int (pc[1]) and coordinate_lc[3] is not pc[3]) or (int (coordinate_lc[1]) <=  int(pc[2]) and int (coordinate_lc[2]) >= int (pc[2]) and coordinate_lc[3] is not pc[3]):
                    antisense = True
                    break
                else:
                    intergenic = True
            else:
                continue
    if sense:
        at = "ID="+str(coordinate_lc[4].get("ID"))+";"+str("Parent="+coordinate_lc[4].get("parent"))+";type=sence_lncRNA"
        return at
    if antisense:
        at = "ID="+str(coordinate_lc[4].get("ID"))+";"+"Parent="+str(coordinate_lc[4].get("parent"))+";type=antisense_lncRNA"
        return at
    if intergenic:
        at = at = "ID="+str(coordinate_lc[4].get("ID"))+";"+"Parent="+str(coordinate_lc[4].get("parent"))+";type=lincRNA"
        return at

pc_coordinates =list()
lc_coordinates =list()
for line1 in open((fn1)):
    if re.match("#", line1):
        continue
    field = line1.rstrip().split("\t")
    if (field[2] == 'mRNA'):
        pc_coordinates.append((field[0],int(field[3]),int(field[4]),field[6]))

for line2 in open((fn2)):
    if re.match("#", line2):
        continue
    gff = gffparse(line2)
    if (gff['type'] == 'gene' or gff['type'] == 'exon'):
        print (gff['scaffold'], gff['source'], gff['type'],gff['start'],gff['end'],gff['score'],gff['strand'],gff['phase'],gff['attribute'],sep='\t' )
    if gff['type'] == 'lncRNA':
        lc_coordinates = ((gff['scaffold'],int(gff['start']),int(gff['end']),gff['strand'],gff['attri'] ))
        attrib = overlap (pc_coordinates , lc_coordinates)
        print (gff['scaffold'], gff['source'], gff['type'],gff['start'],gff['end'],gff['score'],gff['strand'],gff['phase'],attrib,sep='\t' )

    



