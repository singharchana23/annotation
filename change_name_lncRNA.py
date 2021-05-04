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
    }
    the_feature['attri'] = dict(item.split("=") for item in fields[8].split(";"))  
    return the_feature
pc_id =list()
for line1 in open((fn1)):
    if re.match("#", line1):
        continue
    field = line1.rstrip().split("\t")
    
    if field[2] == "gene":
        my_gff = gffparse(line1)
        id_num = my_gff["attri"].get("ID").split("_")
        pc_id.append(int(id_num[1]))
pc_id.sort()
pc_id_min = pc_id[0]
pc_id_max = int(pc_id[-1])
start_id = pc_id_max + 1000
id_pc = str(start_id)
names =list()

for line2 in open((fn2)):
    if re.match("#", line2):
        continue
    field2 = line2.rstrip().split("\t")
    if field2[2] == 'gene':
        name = re.search("Name=([^;]+)",field2[8])
        #parent = re.search("Parent=([^;]+)",field2[8])
        id_lc = re.search("ID=([^;]+)",field2[8]).group(1)
        
        gene_dict = {}
        gene_dict ={
            id_lc : {
                'new_name' : start_id + len(gene_dict)*10,

            }
        }
        attributes = "ID=Smp_"+str(gene_dict[id_lc]['new_name'])+";Name=Smp_"+str(gene_dict[id_lc]['new_name'])
        print (field2[0], 'WBPS', field2[2], field2[3], field2[4], field2[5], field2[6], field2[7], attributes, sep='\t')
        start_id +=10
        count_mrna = 1
        start_mrna = ""
        #start_exon = ""
    
    
    if field2[2] == 'mRNA':
        #count_exon =1
        parent = re.search("Parent=([^;]+)",field2[8]).group(1)
        id_lc = re.search("ID=([^;]+)",field2[8]).group(1)
        if (parent == start_mrna):
            count_mrna+=1
        mrna_dict ={
            id_lc : {
                'new_name': gene_dict[parent]['new_name'],
                'exon' : 0
            }
        }
        attributes = "ID=Smp_"+str(mrna_dict[id_lc]['new_name'])+"."+str(count_mrna)+";Parent=Smp_"+str(mrna_dict[id_lc]['new_name'])
        for_exon = ";Parent=Smp_"+str(mrna_dict[id_lc]['new_name'])+"."+str(count_mrna)
        print (field2[0], 'WBPS', 'lncRNA', field2[3], field2[4], field2[5], field2[6], field2[7], attributes, sep='\t')
        start_mrna = parent

    if field2[2] == 'exon':
        parent = re.search("Parent=([^;]+)",field2[8]).group(1)
        #if (parent == start_exon):
            #count_exon+= 1
        attributes = "ID=Smp_"+str(mrna_dict[parent]['new_name'])+"."+str(count_mrna)+for_exon
        print (field2[0], 'WBPS', field2[2], field2[3], field2[4], field2[5], field2[6], field2[7], attributes, sep='\t')
        #start_exon = parent
    




# mrna_dict ={
#     'Sarah_35207_2_15.3.1' : {
#         'new_name': gene_dict['RLOC_00000001']['new_name']+str(gene_dict['RLOC_00000001']['mrna']+1),
#     'exon' : 0,
#     }
# }

