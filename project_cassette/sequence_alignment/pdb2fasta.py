import sys

three_to_one = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
                'PHE':'F','GLY':'G','HIS':'H','ILE':'I',
                'LYS':'K','LEU':'L','MET':'M','ASN':'N',
                'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
                'THR':'T','VAL':'V','TRP':'W','TYR':'Y',
                'DA':'','DC':'','DT':'','DG':'','DU':''}
                

pdb_name = sys.argv[1]
last_res_num=-999
seq = ''
with open(pdb_name) as of:
    for line in of:
        if line[:5] != 'ATOM ':
            continue
        #print line
        this_res_num = int(line[23:28])
        if this_res_num != last_res_num:
            #print this_res_num
            last_res_num = this_res_num
            resname = line[17:20].strip()
            seq += three_to_one[resname]
print '>' + pdb_name.split('/')[-1].replace('.pdb','')
print seq
