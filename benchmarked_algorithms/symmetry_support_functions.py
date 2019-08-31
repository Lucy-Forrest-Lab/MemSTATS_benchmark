# Name: symmetry_support_functions.py
# Author: Antoniya Aleksandrova
# Date: 22 May 2019
# Language: Python 3.5
# All the functions used to read the raw symmetry information

import os
import numpy as np
import pickle as pkl
import json
import subprocess
import operator
from shutil import copyfile, copy2
import re
import time
start_time = time.time()

import argparse
from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, Selection, Select
from Bio.PDB.Polypeptide import is_aa

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.PDB.vectors import *
from itertools import combinations


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if all([x==0 for x in vector])==True:
      print("Warning: the zero vector is being processed by the unit_vector function.\n")
      return vector
    else:
      return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle



def from3to1_general(resname):
    f3t1 = {'ALA' : 'A',
            'ARG' : 'R',
            'ASN' : 'N',
            'ASP' : 'D',
            'CYS' : 'C',
            'GLN' : 'Q',
            'GLU' : 'E',
            'GLY' : 'G',
            'HIS' : 'H',
            'ILE' : 'I',
            'LEU' : 'L',
            'LYS' : 'K',
            'MET' : 'M',
            'PHE' : 'F',
            'PRO' : 'P',
            'SER' : 'S',
            'THR' : 'T',
            'TRP' : 'W',
            'TYR' : 'Y',
            'VAL' : 'V',
            'MSE' : 'M'}

    if resname in list(f3t1.keys()):
        return f3t1[resname]
    else:
        return '0'


######### Alignment functions 

def get_pdb_sequence_with_chains(structure):
    """
    Retrieves the AA sequence from a PDB structure. It's a list that looks like [(5, 'R', 'A'), (6, 'E', 'A'), (7, 'H', 'A'), (8, 'W', 'A'),...]
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'),r.get_parent().get_id(),r.id[0],r.id[2])
    seq = [_aainfo(r) for r in structure.get_residues() if (is_aa(r) and r.has_id('CA'))]
    return seq


def parse_structure(spath):
    """Parses a PDB/cif structure"""

    if not os.path.isfile(spath):
        return IOError('File not found: {0}'.format(spath))

    if spath.endswith(('pdb', 'ent')):
        parser = PDBParser(QUIET=True)
    elif spath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(spath))

    sname = os.path.basename(spath.split('.')[0])
    return parser.get_structure(sname, spath)


def strip_tm_chains(wkdir, inputf, oriented_enc, chains):
    """
    Extract only the transmembrane chains atoms that are properly specified and order them as they appear in the PDB file.
    Note that chains can be either an array/list or a string.
    """
    f = open(oriented_enc, 'r')
    altloc = ' '
    flag = 0
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip() == 'CA' and line[16:17] != ' ' and (
                float(line[54:60]) > 0.5 or flag == 0):
            altloc = line[16:17]
            flag = 1
    f.seek(0)
    o = open(wkdir + inputf + "_tmp.pdb", "w")
    num = 0
    LINELEM = "{:76s}{:>2s}\n"
    atomnames = []
    old_resid = '-88'
    for line in f:
        if (line.startswith("ATOM") or line.startswith("TER ")) and (
                line[21:22] in chains or chains == "") and from3to1_general(line[17:20].strip()) != 0:
            if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                    line[16:17] != altloc and line[16:17] != ' '):  # sort out disordered atoms
                continue
            elif (old_resid == '-88' or line[22:27] != old_resid):
                old_resid = line[22:27]
                atomnames = []
            atomnames.append(line[12:16])

            if line[76:78].strip() != '':  # ensure that the element symbol is included
                o.write(line)
                num += 1
            else:
                if line.startswith("TER "):
                    o.write(line)
                else:
                    atom = line[12:16].strip()
                    elem = atom[0]
                    o.write(LINELEM.format(line[0:76], elem))
                    num += 1
        if line.startswith("HETATM") and line[17:20] == 'MSE' and (line[21:22] in chains or chains == ""):
            if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                    line[16:17] != altloc and line[16:17] != ' '):
                continue
            elif (old_resid == '-88' or line[22:27] != old_resid):
                old_resid = line[22:27]
                atomnames = []
            atomnames.append(line[12:16])

            if line[76:78].strip() != '':
                o.write("ATOM  " + line[6:])
                num += 1
            else:
                atom = line[12:16].strip()
                elem = atom[0]
                o.write("ATOM  " + line[6:76] + " " + elem + "\n")
                num += 1

    o.write("END\n")
    f.close()
    o.close()
    return num

def strip_tm_chains_in_order(wkdir,inputf,oriented_enc,chains):
    """
    Extract only the transmembrane chains atoms that are properly specified and order the chains as requested in the variable chains.
    Note that chains is an array here, eg. ['A','B','C'].
    """
    f=open(oriented_enc,'r')
    altloc=' '
    flag=0
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip()=='CA' and line[16:17]!=' ' and (float(line[54:60])>0.5 or flag==0):
            altloc=line[16:17]
    f.seek(0)
    o=open(wkdir+inputf+"_tmp.pdb","w")
    num=0
    LINELEM="{:76s}{:>2s}\n"
    for chain in chains:
        f.seek(0)
        atomnames=[]
        old_resid='-88'
        for line in f:
            if (line.startswith("ATOM") or line.startswith("TER ")) and line[21:22]==chain and from3to1_general(line[17:20].strip())!=0:
                if old_resid!='-88' and line[22:27]==old_resid and line[12:16] in atomnames or (line[16:17]!=altloc and line[16:17]!=' '): #sort out disordered atoms
                    continue
                elif (old_resid=='-88' or line[22:27]!=old_resid):
                    old_resid=line[22:27]
                    atomnames=[]
                atomnames.append(line[12:16])

                if line[76:78].strip()!='': #ensure that the element symbol is included
                    o.write(line)
                    num+=1
                else:
                    if line.startswith("TER "):
                        o.write(line)
                    else:
                        atom=line[12:16].strip()
                        elem=atom[0]
                        o.write(LINELEM.format(line[0:76],elem))
                        num+=1
            if line.startswith("HETATM") and line[17:20]=='MSE' and line[21:22]==chain:
                if old_resid!='-88' and line[22:27]==old_resid and line[12:16] in atomnames or (line[16:17]!=altloc and line[16:17]!=' '):
                    continue
                elif (old_resid=='-88' or line[22:27]!=old_resid):
                    old_resid=line[22:27]
                    atomnames=[]
                atomnames.append(line[12:16])

                if line[76:78].strip()!='':
                    o.write("ATOM  "+line[6:])
                    num+=1
                else:
                    atom=line[12:16].strip()
                    elem=atom[0]
                    o.write("ATOM  "+line[6:76]+" "+elem+"\n")
                    num+=1
	   
    o.write("END\n")
    f.close()
    o.close()
    return num

#### Reading CE-Symm and SymD information ####
def cesymm_data(cesymm_out_dir, inputf):
    """
    Crawl through the CE-Symm result files and collect all the necessary information for a structure
    """
    ce_order='na'
    sym_level='na'
    ce_rmsd='na'
    tmScore='na'
    cce_type='na'
    num_repeats='na'
    internal="na"
    quaternary="na"
    angle="na"
    axis_angle="na"
    detected_symmetry='na'
    coverage='na'
    length='na'
    core_length='na'
    translation='na'
    unrefined_tmscore='na'
    unrefined_rmsd='na'
    seed='na'
    symm_type='na'
    repeats='na'
    repeat_length='na'
    closed_opened='na'
    repeat_level='na'
    if os.path.isfile(cesymm_out_dir+inputf+".axes") and os.path.isfile(cesymm_out_dir+inputf+"_stdout.out") and os.path.getsize(cesymm_out_dir+inputf+"_stdout.out") > 5:
        stdout=open(cesymm_out_dir+inputf+"_stdout.out","r")
        for ln in stdout:
            if inputf[0:4] in ln:
                fields=ln.split()
                num_repeats=fields[1] #17
                ce_order=fields[2] #18
                ind=0
                try:
                    sym_level=int(fields[4])
                except ValueError:
                    ind=-1
                    sym_level=int(fields[ind+4]) #20
                    ce_order='na'
                closed_opened=fields[ind+5]
                tmScore=float(fields[ind+10]) #24
                ce_rmsd=float(fields[ind+11]) #25
                unrefined_tmscore=float(fields[ind+8])
                unrefined_rmsd=float(fields[ind+9])
                repeat_length=int(fields[ind+12])
                core_length=int(fields[ind+13]) #29
                length=int(fields[ind+14]) #30
                coverage=fields[ind+15] #31
                if length > 0 and float(coverage)-core_length/length > 0.01: # correct errors in CE-Symm coverage calculation
                    coverage=str(round(core_length/length,2))
                cce_type="na"
                angle=fields[ind+6]
                internal="No"
                quaternary="No"
                axis_angle="na"
                detected_symmetry='na' 
                translation=fields[ind+7]
                with open(cesymm_out_dir+inputf+".seed","r") as f:
                        seed=f.readline().strip() 
                symm_type="na"
                repeats="na" 
                repeat_level="na"
                if sym_level>0:
                      cce_type=""
                      angle=""
                      axis_angle=""
                      translation=""
                      symm_type=""
                      repeats=""
                      closed_opened=""
                      repeat_level=""
                      f=open(cesymm_out_dir+inputf+".axes","r")
                      for line in f:
                        if inputf[0:4] in line:
                          flds=line.split()
                          closed_opened=closed_opened+flds[2]+';'
                          angle=angle+flds[4]+";"
                          translation=translation+flds[5]+";"
                          repeats=repeats+flds[8].replace(";",",")+";"
                          repeat_level=repeat_level+flds[1]+";"
                          pt1=[float(x) for x in flds[6].split(",")]
                          pt1=np.array(pt1)
                          pt2=[float(x) for x in flds[7].split(",")]
                          pt2=np.array(pt2)
                          u_tmp=pt1-pt2
                          u=unit_vector(u_tmp)
                          axis_angle_tmp=angle_between(u,(0,0,1))*180/np.pi
                          if axis_angle_tmp>90:
                            axis_angle_tmp=axis_angle_tmp-180
                          axis_angle=axis_angle+str(float('%.2f'% round(axis_angle_tmp,2)))+';'
                          if np.abs(np.dot(u, (0,0,1)))>0.5:
                            cce_type=cce_type+"Parallel;"
                          else:
                            cce_type=cce_type+"Antiparallel;"
                          rep=flds[8].split(")")
                          chain=[rep[0][pos+1] for pos, char in enumerate(rep[0]) if char == '.']	#finds all chain names in first symmetry bracket (they are preceded by .)
                          if chain.count(chain[0])!=len(chain): #if all repeats are from the same chain, no quaternary structure symmetry was detected
                            detected_symmetry="Quaternary"
                            quaternary='Yes'
                          else:
                            detected_symmetry="Internal"
                            internal='Yes'
                          symm_type=symm_type+detected_symmetry+";"	
                      f.close()
                      repeats=repeats.replace(inputf[0:4].upper()+".","")
                      
	
    else:
      print("CE-Symm files missing for %s" % inputf)


    cesymm_dic={}
    cesymm_dic['size']=length
    cesymm_dic['coverage']=coverage
    cesymm_dic['unit_angle']=angle
    cesymm_dic['unit_translation']=translation
    cesymm_dic['refined_tmscore']=str(tmScore)
    cesymm_dic['refined_rmsd']=str(ce_rmsd)
    cesymm_dic['unrefined_tmscore']=str(unrefined_tmscore)
    cesymm_dic['unrefined_rmsd']=str(unrefined_rmsd)
    cesymm_dic['repeats_number']=num_repeats
    cesymm_dic['topology']=cce_type
    cesymm_dic['symmetry_order']=ce_order
    cesymm_dic['symmetry_levels']=str(sym_level)
    cesymm_dic['aligned_length']=core_length
    cesymm_dic['seed']=seed
    cesymm_dic['internal_symmetry']=internal
    cesymm_dic['quaternary_symmetry']=quaternary
    cesymm_dic['axis_angle_with_membrane_normal']=axis_angle
    cesymm_dic['symmetry_type']=symm_type
    cesymm_dic['repeats']=repeats
    cesymm_dic['repeat_length']=repeat_length 
    cesymm_dic['closed_opened']=closed_opened
    cesymm_dic['repeat_level']=repeat_level
    return cesymm_dic
    

def symd_data(symd_out_dir,inputf):
    """
    Crawl through the SymD result files and collect all the necessary information for a structure
    """
    if os.path.isfile(symd_out_dir+inputf+"-info.txt") and os.path.isfile(symd_out_dir+inputf+"-trfm.pdb"):
        f=open(symd_out_dir+inputf+'-info.txt','r')
        for line in f:
            if line.startswith("! Protein Size:"):
                flds=line.split()
                size=int(flds[3])
            if line.startswith("! Best IS:"):
                flds=line.split(':')
                best_is=flds[1].strip()
            if line.startswith("! N-aligned at best IS:"):
                flds=line.split(':')
                aligned=int(flds[1])
                align_per=str(round(1.0*aligned/size,2))
            if line.startswith("! TM/Nr at best IS:"):
                flds=line.split(':')
                symd_tm=float(flds[1])
                symd_tm=str(round(symd_tm,2))
            if line.startswith("! Z(TsC), Z(Ts), Z(TM) at best IS:"):
                flds=line.split(',')
                z_tm=flds[-1].strip()
            if line.startswith("! RMSD at best IS:"):
                flds=line.split(':')
                symd_rmsd=float(flds[1])  
                symd_rmsd=str(round(symd_rmsd,3))  
            if line.startswith("! Highest-scoring angle:"):
                flds=line.split()
                symd_angle=float(flds[3])
            if line.startswith("! Highest-scoring p-transl:"):
                flds=line.split()
                symd_translation=flds[3].strip()
            if line.startswith("! Derived unit angle:"):
                flds=line.split()
                symd_unit_angle=float(flds[4])
                if symd_unit_angle!=0:
                    symd_order=int(round(360.0/symd_unit_angle))
                else:
                    symd_order='na'
            if line.startswith("! Derived unit p-transl:"):
                flds=line.split()
                symd_unit_translation=flds[4].strip() 
                break
        f.close()
        f=open(symd_out_dir+inputf+'-trfm.pdb','r')
        counter=0
        flag=0
        pts=[]
        for line in f:
            if line.startswith("MODEL        3"):
                flag=1
            if flag==1 and 0<=counter<2 and line.startswith("ATOM") and line[13:16].strip()=='CA':
                pts.append([float(line[30:38]), float(line[38:46]),float(line[46:54])])
                counter+=1
        f.close()
        if len(pts)==2:                 
            u=np.array(pts[0])-np.array(pts[1])
            u=unit_vector(u)
            axis_angle=angle_between(u,(0,0,1))*180/np.pi
            if axis_angle>90:
                axis_angle=axis_angle-180
            axis_angle=str(float('%.2f' % axis_angle))
            if np.abs(np.dot(u, (0,0,1)))>0.5:
                symd_type="orthogonal"
            else:
                symd_type="parallel"
        else:
            symd_order,symd_rmsd,symd_angle,symd_type,align_per,aligned,symd_tm,z_tm,symd_unit_translation,symd_translation,symd_unit_angle,size,best_is,axis_angle='na','na','na','na','na','na','na','na','na','na','na','na','na','na','na'

    else:
        print("Missing SymD files for %s" % inputf)
        symd_order='na'
        symd_rmsd='na'
        symd_angle='na'
        symd_type='na'
        align_per='na'
        aligned='na'
        symd_tm='na'
        z_tm='na'
        symd_unit_translation='na'
        symd_translation='na'
        symd_unit_angle='na'
        size='na'
        best_is='na'
        axis_angle='na'
    
    symd_dic={}
    symd_dic['size']=str(size)
    symd_dic['coverage']=str(align_per)
    symd_dic['unit_angle']=str(symd_unit_angle)
    symd_dic['unit_translation']=symd_unit_translation
    symd_dic['tmscore']=symd_tm
    symd_dic['rmsd']=symd_rmsd
    symd_dic['initial_shift']=best_is
    symd_dic['z_tmscore']=z_tm
    symd_dic['is_angle']=str(symd_angle)
    symd_dic['is_translation']=symd_translation
    symd_dic['aligned_length']=str(aligned)
    symd_dic['symmetry_order']=str(symd_order)
    symd_dic['topology']=symd_type
    symd_dic['axis_angle_with_membrane_normal']=axis_angle
    return symd_dic


def ananas_data(out_dir, inputf, oriented_struct):
    """ Read AnAnaS raw data and select highest-order symmetry that includes the chains present in the structure"""
    fname = out_dir + inputf + "_ananas"
    ananas_dic = {}
    if os.path.isfile(fname + ".json") and len(open(fname + ".json").readlines( )) > 2 and os.path.isfile(fname + ".out"):
        # Figure out the names and residue ranges of the chains and connect these to the AnAnaS numbering
        resseq_A = get_pdb_sequence_with_chains(oriented_struct)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
        chain_resids = {} # eg. {0: ('A', 2, 334), 1:('B', 5, 175)}
        chains = []
        [chains.append(i[2]) for i in resseq_A if i[2] not in chains]
        for i, chain in enumerate(chains):
            residues = list(filter(lambda x: x[2] == chain, resseq_A))
            chain_resids[i] = (chain, residues[0], residues[-1])

        # Read info from AnAnaS
        info = open(fname + ".out", 'r')
        for line in info:
            if line.startswith("Number of chains read"):
                chains_num = int(line.strip().split()[-1])
        info.close()
        with open(fname + ".json") as jf:
            data = json.load(jf)
        for i,symm in enumerate(data):
            if 'c' in symm['group'] and (int(symm['group'][1:]) > chains_num or (int(symm['group'][1:]) < chains_num and chains_num % int(symm['group'][1:]) !=0)): #e.g. c9 to indicate symmetry group C9
                print('WARNING: [AnAnaS] Inconsistent chain number or cyclic order')
                continue
            elif 'd' in symm['group'] and (int(symm['group'][1:])*2 > chains_num or (int(symm['group'][1:])*2 < chains_num and chains_num % int(symm['group'][1:])*2 !=0)): #e.g. d6 to indicate symmetry group D6
                print('WARNING: [AnAnaS] Inconsistent chain number or dihedral order')
                continue
            ananas_dic['symmetry_order'] = symm['group'].upper() # the global order (e.g. D6)
            ananas_dic['avg_rmsd'] = symm['Average_RMSD']
            ananas_dic['symmetries'] = []

            # Find one (the first) highest order permutation
            for j, t in enumerate(symm['transforms']):
                level={}
                if t['ORDER'] == int(ananas_dic['symmetry_order'][1:]):
                    level['order'] = t['ORDER']

                    # estimate topology and angle with membrane
                    u = np.array(t['P0']) - np.array(t['P1'])
                    level['topology'], axis_angle = determine_repeat_topology(u)

                    # Trace which chain corresponds to which and form the symmetry-related tuples
                    repeats_tuples = [] # eg. [[0,2],[1,3]]
                    for m in range(len(t['permutation'])):
                        tup = []
                        if len([True for k in repeats_tuples if (m in k or t['permutation'][m] in k)]) == 0:
                            tup.append(m)
                            ind = t['permutation'][m]
                            while ind not in tup:
                                tup.append(ind)
                                ind = t['permutation'][ind]
                            repeats_tuples.append(tup)
                    level['repeats_chain_indexes'] = np.transpose(repeats_tuples)
                    level['repeats'] = [[(chain_resids[m][0], chain_resids[m][1][0], chain_resids[m][2][0]) for m in ind] for ind in level['repeats_chain_indexes']] # eg. [[('A', 1, 453), ('B', 44, 295)], [('C', 1, 453), ('D', 44, 295)]], i.e. [Rep1, Rep2]
                    repeats_text = [[chain_resids[m][0] + '_' + str(chain_resids[m][1][0])
                                     + '-' + str(chain_resids[m][2][0]) for m in ind] for ind in repeats_tuples]
                    level['repeats_text'] = ''.join(['(' + ','.join(m) + ')' for m in repeats_text]) # CE-Symm compatible format, i.e. (A_1-453,C_1-453)(B_44-295,D_44-295)

                    level['symm_units'] = 'Quaternary'
                    level['group'] = 'Cyclic'
                    level['permutation'] = t['permutation']
                    ananas_dic['symmetries'].append(level)
                    break
            # for dihedral symmetries, extract the two-fold symmetries
            if ananas_dic['symmetry_order'][0] == 'D':
                chains_per_repeat = chains_num // (int(ananas_dic['symmetry_order'][1:]) * 2)
                two_fold_symms = np.array(ananas_dic['symmetries'][0]['repeats_chain_indexes'])
                two_fold_symms = two_fold_symms.reshape(two_fold_symms.shape[0], two_fold_symms.shape[1]//chains_per_repeat, chains_per_repeat)
                for correspondence in two_fold_symms:
                    for t in symm['transforms']:
                        if [t['permutation'][m] for m in correspondence[0]] == list(correspondence[1]):
                            level = {}
                            level['order'] = t['ORDER']
                            assert level['order'] == len(correspondence)

                            # estimate topology and angle with membrane
                            u = np.array(t['P0']) - np.array(t['P1'])
                            level['topology'], axis_angle = determine_repeat_topology(u)
                            level['repeats'] = [[(chain_resids[m][0], chain_resids[m][1][0], chain_resids[m][2][0]) for m in r] for r in correspondence]
                            repeats_text = [[chain_resids[m][0] + '_' + str(chain_resids[m][1][0])
                                        + '-' + str(chain_resids[m][2][0]) for m in r] for r in correspondence]
                            repeats_text = np.transpose(repeats_text)
                            level['repeats_text'] = ''.join(['(' + ','.join(m) + ')' for m in repeats_text]) # CE-Symm compatible format
                            level['symm_units'] = 'Quaternary'
                            level['group'] = 'Cyclic'
                            level['permutation'] = t['permutation']
                            ananas_dic['symmetries'].append(level)
                            break
            break
    else:
        if os.path.isfile(fname + ".json") and len(open(fname + ".json").readlines( )) <= 2:
            print("AnAnaS found no symmetry %s" % inputf)
        else:
            print("Missing AnAnaS files for %s" % inputf)
    if len(ananas_dic) == 0:
        ananas_dic = {'symmetry_order': None, 'avg_rmsd': None, 'symmetries': [{'order': None, 'group': None, 'topology': None, 'repeats': None, 'symm_units': None}]}

    return ananas_dic

def determine_repeat_topology(u):
    u = unit_vector(u)
    axis_angle = angle_between(u, (0, 0, 1)) * 180 / np.pi
    if axis_angle > 90:
        axis_angle = axis_angle - 180
    if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
        return "Parallel", str(float('%.2f' % axis_angle))
    else:
        return "Antiparallel", str(float('%.2f' % axis_angle))

def get_repeat_resid(cesymm_out_dir,inputf,reference):
    """
    CE-Symm lists repeats in a way that makes it unclear where the repeat ends in multi-chain structures. 
    For example, 2A65.A_9-235 indicates that the repeat starts with residue 9 in chain A but the end might 
    be residue 235 of chain C, encompassing chain B on the way. Therefore, we need to double-check the 
    repeat listed in the fasta alignments with the sequence from the pdb file and determine what ranges of 
    residues constitute each repeat. Here, algn_map is an array of the type 
    [[(9,'A'),(235,'A'),(1,'B'),(235,'B'),(1,'C'),(235,'C')],[...]], so that within each [] we have the ranges 
    for one particular repeat. 
     
    """
    resseq_A = get_pdb_sequence_with_chains(reference) # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm=list(SeqIO.parse(cesymm_out_dir+inputf+".fasta", "fasta"))
    repeat_start=[]
    algn_map=[[] for i in range(len(cesm))]
    for num,r in enumerate(cesm):
        r_id=r.id[len(inputf[0:4]):]
        r_id=r_id.split('_')
        ch=r_id[0][-1]
        repeat=r_id[1].strip()
        first_digit=repeat[0] # handling negative residues
        repeat=repeat[1:].split('-')
        repeat[0]=first_digit+repeat[0]

        repeat_start=[i for i, sublist in enumerate(resseq_A) if (int(repeat[0])==sublist[0] and ch==sublist[2])] #finds the index of the residue that starts the repeat
        count = 0
        small_lett = 0
        for i,c in enumerate(r.seq):
            if c!='-' and c!='/':
              if c.isupper() and resseq_A[repeat_start[0]+count][1] == c and count==0:
                algn_map[num].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
              elif resseq_A[repeat_start[0]+count][2]!=resseq_A[repeat_start[0]+count-small_lett-1][2] and c.isupper() and resseq_A[repeat_start[0]+count][1] == c:
                algn_map[num].append((resseq_A[repeat_start[0]+count-small_lett-1][0],resseq_A[repeat_start[0]+count-small_lett-1][2]))
                algn_map[num].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
              elif len(r.seq)==i+1 or r.seq[i+1]=='/': # note that small letters at the endo of repeats are not counted by CE-Symm in the repeat id (it's a bug!) - they should be
                algn_map[num].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
                break
              count+=1
              if c.isupper():
                  small_lett = 0
              else:
                  small_lett+=1
            if c=='-' and (len(r.seq)==i+1 or r.seq[i+1]=='/') and len(algn_map[num])%2==1: # handles cases like ...KAF--
                algn_map[num].append((resseq_A[repeat_start[0]+count-1][0],resseq_A[repeat_start[0]+count-1][2]))
        if len(algn_map[num])%2==1:
            algn_map[num].append(algn_map[num][-1])
    return algn_map  

def get_symd_aligned_resid(symd_out_dir,inputf,reference):
    resseq_A = get_pdb_sequence_with_chains(reference) # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    symd_fasta=SeqIO.to_dict(SeqIO.parse(symd_out_dir+inputf+"-best.fasta", "fasta"))
    i=0
    aligned=[]
    aligned_chains=[]
    end='0'  
    flag=0
    for c in symd_fasta[inputf+'-original'].seq:
        if c!='-':
            if c.isupper():
                if resseq_A[i][2] not in aligned_chains:
                    if end!='0' and last_char.isupper():
                        aligned[ind].append(str(resseq_A[i-1][0]))
                    aligned_chains.append(resseq_A[i][2])
                    if i==(len(resseq_A)-1):  #this provisions for a single capital letter at the very end of the sequence
                        aligned.append([str(resseq_A[i][0]),str(resseq_A[i][0])])
                    else:
                        aligned.append([str(resseq_A[i][0])])
                        ind=aligned_chains.index(resseq_A[i][2])
                        end=str(resseq_A[i][0])
                    flag=1
                else:
                    ind=aligned_chains.index(resseq_A[i][2])
                    end=str(resseq_A[i][0])
                    flag=1
                    if i==(len(resseq_A)-1):
                        aligned[ind].append(end)
                    if last_char.islower():  #accounts for unaligned residues (breaks) in a single chain
                        aligned[ind].append(end)
                        flag=1
            i+=1
            last_char=c
        if c.islower() and end!='0' and flag==1:
            aligned[ind].append(end)
            flag=0   
    for i,x in enumerate(aligned): # handles a single capital letter in the beginning of a chain sequence, eg. 3tdo
        if len(x)%2==1:
            aligned[i].append(x[-1])  
    return aligned, aligned_chains # a list of the residue ids and a list of the chains of all the capital letter residues in the alignment
                                   # eg: aligned=[[2,10,20,33],[1,1,5,18]] aligned_chains=['A','L']

