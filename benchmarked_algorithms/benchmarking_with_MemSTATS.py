# Name: benchmarking_with_MemSTATS.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5
# Date: 22 May 2019
# Description: Takes each representative symmetry from MemSTATS_benchmark.xlsx
#              and assesses whether CE-Symm 2.2, SymD 1.6, AnAnaS and EncoMPASS identify it correctly.
#              Make sure to define the paths of key folders where indicated below
# Prerequisites: All output files from processing the MemSTATS pdbs with the algorithms
#                Information about the transmembrane regions of the proteins obtained with PPM (Lomize, AL et. al. 2011, Lomize, MA et al. 2012)
#                Files can be found at
# Usage: python benchmark_review_symmetry.py


import pandas as pd
import datetime
import copy
import pickle as pkl
from symmetry_support_functions import *
import glob

def read_analysis_files(cesymm_out_dir, inputf):
    analysis_dic = {}
    cce_type = ""
    angle = ""
    axis_angle = ""
    translation = ""
    symm_type = ""
    repeats = ""
    closed_opened = ""
    repeat_level = ""
    if os.path.isfile(cesymm_out_dir + inputf + ".axes"):
        f = open(cesymm_out_dir + inputf + ".axes", "r")
        for line in f:
            if inputf[0:4] in line:
                flds = line.split()
                closed_opened = closed_opened + flds[2] + ';'
                angle = angle + flds[4] + ";"
                translation = translation + flds[5] + ";"
                repeats = repeats + flds[8].replace(";", ",") + ";"
                repeat_level = repeat_level + flds[1] + ";"
                pt1 = [float(x) for x in flds[6].split(",")]
                pt1 = np.array(pt1)
                pt2 = [float(x) for x in flds[7].split(",")]
                pt2 = np.array(pt2)
                u_tmp = pt1 - pt2
                u = unit_vector(u_tmp)
                axis_angle_tmp = angle_between(u, (0, 0, 1)) * 180 / np.pi
                if axis_angle_tmp > 90:
                    axis_angle_tmp = axis_angle_tmp - 180
                axis_angle = axis_angle + str(float('%.2f' % round(axis_angle_tmp, 2))) + ';'
                if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
                    cce_type = cce_type + "Parallel;"
                else:
                    cce_type = cce_type + "Antiparallel;"
                rep = flds[8].split(")")
                chain = [rep[0][pos + 1] for pos, char in enumerate(rep[0]) if
                         char == '.']  # finds all chain names in first symmetry bracket (they are preceded by .)
                if chain.count(chain[0]) != len(
                        chain):  # if all repeats are from the same chain, no quaternary structure symmetry was detected
                    detected_symmetry = "Quaternary"
                else:
                    detected_symmetry = "Internal"
                symm_type = symm_type + detected_symmetry + ";"
        f.close()
        repeats = repeats.replace(inputf[0:4].upper() + ".", "")
        repeats = repeats.replace(inputf[0:4] + ".", "")
    else:
        print("Analysis files missing (looking for %s in %s)" % (inputf, cesymm_out_dir))
        cce_type = 'na'
        angle = "na"
        axis_angle = "na"
        translation = 'na'
        symm_type = 'na'
        repeats = 'na'
        closed_opened = 'na'
        repeat_level = 'na'

    analysis_dic['unit_angle']=angle
    analysis_dic['unit_translation']=translation
    analysis_dic['topology']=cce_type
    analysis_dic['axis_angle_with_membrane_normal']=axis_angle
    analysis_dic['symmetry_type']=symm_type
    analysis_dic['repeats']=repeats
    analysis_dic['closed_opened']=closed_opened
    analysis_dic['repeat_level']=repeat_level
    return analysis_dic

def benchmark_set_symmetry(fields, col_names):
    # Takes a line from the benchmark set table and returns a list of dictionaries for each symmetry in the entry
    bs_symm = []
    order = fields[col_names['Order']].split(';')
    topology = fields[col_names['Repeat-Topology']].split(';')
    group = fields[col_names['Closed/Open-Symmetry-Group']].split(';')
    symm_units = fields[col_names['Internal/Quaternary']].split(';')
    described = fields[col_names['Described']].split(';')
    repeats = fields[col_names['Repeats']].split(';') # eg: repeats = ['(C_2-280,M_2-280,Q_2-280)(D_-1-190,N_-1-190,R_-1-190)','(C_1-90,C_102-191)']
    interdig = fields[col_names['Interdigitating']].split(';')
    assert len(order)==len(topology) and len(order)==len(group) and len(order)==len(symm_units) and len(order)==len(symm_units) and len(order)==len(described), "Inconsistencies in the benchmark set."
    repeats_list = [] # we want it in the form eg: [[[('C',2,280),('D',-1,190)],[('M',2,280),('N',-1,190)],[('Q',2,280),('R',-1,190)]], [[('C',1,90)],[('C',102,191)]]], i.e. All <Symmetry < Repeat1, Repeat2, Repeat3
    if order[0] == 'None':
        bs_symm.append({'order': None, 'topology': None, 'group': None, 'symm_units': None,
                        'described': described[0], 'repeats': None, 'interdigitating': None, 'repeats_text': repeats})
    else:
        for i in range(0,len(order)):
            reps = repeats[i].strip(')/(/\n')
            reps = reps.split(')(')  # eg: reps = ['C_2-280,M_2-280,Q_2-280' , 'D_-1-190,N_-1-190,R_-1-190']
            repeats_list.append([[] for z in range(len(reps[0].split(',')))]) # repeats_list = [[[],[],[]]]
            for pair in reps:
                entry = pair.split(',')  # eg: entry = ['C_2-280','M_2-280','Q_2-280']
                for ind, k in enumerate(entry):
                    ch = k.split('_')[0][0]
                    beg = k.split('_')[1]
                    first_digit = beg[0]  # in case the residue has a negative id (i.e. -1)
                    end = beg[1:].split('-')[1] # eg: end = 280
                    beg = beg[1:].split('-')[0]
                    beg = first_digit + beg  # eg: beg = 2
                    repeats_list[i][ind].append((ch,int(beg),int(end)))

            bs_symm.append({'order':int(order[i]), 'topology':topology[i], 'group':group[i],'symm_units':symm_units[i],
                            'described':described[i], 'repeats':repeats_list[i], 'interdigitating': interdig[i],
                            'repeats_text': repeats[i]})
    return bs_symm

def cesymm_symmetry(pdb, tm_chains, order_matters, locations, pdb_suff):
    cesymm_out_dir = locations['FSYSPATH']['cesymm']
    if os.path.isfile(cesymm_out_dir+pdb+"_stdout.out"):
        cesymm_dic = cesymm_data(cesymm_out_dir, pdb)
    else:
        cesymm_dic = read_analysis_files(cesymm_out_dir,pdb)
    cs_symm = []
    if cesymm_dic['topology'] == 'na' or ('symmetry_order' in cesymm_dic and (cesymm_dic['symmetry_order']=='C1' or (cesymm_dic['symmetry_levels']!='na' and int(cesymm_dic['symmetry_levels'])<1))):
        cs_symm.append({'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None})
    else:
        topology = cesymm_dic['topology'][:-1].split(';')
        repeats = cesymm_dic['repeats'][:-1].split(';')
        repeats = cesymm_dic['repeats'][:-1].split(';')
        closed_opened = cesymm_dic['closed_opened'][:-1].split(';')
        symmetry_type = cesymm_dic['symmetry_type'][:-1].split(';')
        repeats_list = []

        # Figure out the actual ranges of each repeat
        wkdir = locations['OPT']['temppath']
        if not os.path.isdir(wkdir):
            os.mkdir(wkdir)
        oriented_opm = locations['FSYSPATH']['whole'] + pdb[0:4] + pdb_suff + '.pdb'
        if order_matters=='yes':
            num = strip_tm_chains_in_order(wkdir, pdb[0:4], oriented_opm, tm_chains)
        else:
            num = strip_tm_chains(wkdir, pdb[0:4], oriented_opm, tm_chains)  # should not strip them in order because they were not run that way
        oriented = wkdir + pdb[0:4] + "_tmp.pdb"
        reference = parse_structure(oriented)[0]
        repeat_selection = get_repeat_resid(cesymm_out_dir,pdb,reference) # eg: [[(9,'A'),(235,'A'),(1,'B'),(235,'B'),(1,'C'),(235,'C')],[...]], so that within each [] we have the ranges
        os.remove(oriented)

        for i in range(0,len(topology)):
            reps = repeats[i].strip(')/(/\n')
            reps = reps.split(')(') # eg: reps = ['C_2-280,M_2-280,Q_2-280' , 'D_-1-190,N_-1-190,R_-1-190']
            repeats_list.append([[] for z in range(len(reps[0].split(',')))]) # repeats_list = [[[],[],[]]]
            order = len(repeats_list[i])
            if closed_opened[i]=='CLOSED':
                group = 'Cyclic'
            else:
                group = 'Repeated'
            for pair in reps:
                entry = pair.split(',') # eg: entry = ['C_2-280','M_2-280','Q_2-280']
                for ind, k in enumerate(entry):
                    ch = k.split('_')[0][0]
                    beg = k.split('_')[1]
                    first_digit = beg[0] # in case the residue has a negative id (i.e. -1)
                    beg = beg[1:].split('-')[0]
                    beg = first_digit + beg # eg: beg = 2
                    locator = [i for i, sublist in enumerate(repeat_selection) if ((int(beg), ch) == sublist[0])]
                    for z in range(0,int(len(repeat_selection[locator[0]]) / 2)):
                        repeats_list[i][ind].append((repeat_selection[locator[0]][z * 2][1],repeat_selection[locator[0]][z * 2][0],repeat_selection[locator[0]][z * 2 + 1][0]))

            cs_symm.append({'order':order, 'topology':topology[i], 'group':group, 'symm_units':symmetry_type[i], 'repeats':repeats_list[i]})


    return cs_symm

def symd_symmetry(pdb, tm_chains, order_matters, locations, pdb_suff, zScoreThreshold):
    symd_out_dir = locations['FSYSPATH']['symd']
    symd_dic = symd_data(symd_out_dir, pdb)
    symd_symm = []
    if symd_dic['topology'] == 'na' or symd_dic['unit_angle'] =='na' or float(symd_dic['z_tmscore']) < zScoreThreshold:
        if symd_dic['unit_angle'] == 'na':
            print('Warnging: SymD results might be missing.')
        symd_symm.append({'order': None, 'topology': None, 'group': None, 'symm_units': None, 'aligned': None})
    else:
        if float(symd_dic['unit_translation']) > 3:
            group = 'Repeated'
        else:
            group = 'Cyclic'

        if symd_dic['topology'] == 'orthogonal':
            topology = 'Parallel'
        else:
            topology = 'Antiparallel'

        # Figure out the actual ranges of the aligned pieces
        wkdir = locations['OPT']['temppath']
        if not os.path.isdir(wkdir):
            os.mkdir(wkdir)
        oriented_opm = locations['FSYSPATH']['whole'] + pdb[0:4] + pdb_suff + '.pdb'
        if order_matters=='yes':
            strip_tm_chains_in_order(wkdir, pdb[0:4], oriented_opm, tm_chains)
        else:
            strip_tm_chains(wkdir, pdb[0:4], oriented_opm, tm_chains)  # should not strip them in order because they were not run that way
        oriented = wkdir + pdb[0:4] + "_tmp.pdb"
        reference = parse_structure(oriented)[0]
        aligned_ids, aligned_chains = get_symd_aligned_resid(symd_out_dir, pdb, reference) #eg: aligned_ids=[[2,10,20,33],[1,1,5,18]] aligned_chains=['A','L']
        os.remove(oriented)
        aligned = []  # eg: aligned = [('A', {2,3,4,5,...,10,20,21,22,...,33}),('L',{1,5,6,7,8,...,18})]
        for i,c in enumerate(aligned_chains):
            ids_set = set()
            for k in range(0,len(aligned_ids[i])//2):
                s = set(range(int(aligned_ids[i][k*2]), int(aligned_ids[i][k*2 + 1])))
                ids_set = ids_set.union(s)
            aligned.append((c,ids_set))

        chain_coverage = [1 for i in aligned if len(i[1]) > 5]
        if len(chain_coverage) > 1:
            symm_units = 'Quaternary'
        else:
            symm_units = 'Internal'

        symd_symm.append({'order': int(symd_dic['symmetry_order']), 'topology': topology, 'group': group, 'symm_units': symm_units, 'aligned': aligned})

    return symd_symm

def ananas_symmetry(pdb, tm_chains, locations, pdb_suff, RmsdThreshold=10):
    ananas_out_dir = locations['FSYSPATH']['ananas']
    wkdir = locations['OPT']['temppath']
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    oriented_opm = locations['FSYSPATH']['whole'] + pdb[0:4] + pdb_suff + '.pdb'
    strip_tm_chains(wkdir, pdb[0:4], oriented_opm, tm_chains)
    oriented = wkdir + pdb[0:4] + "_tmp.pdb"
    oriented_struct = parse_structure(oriented)[0]
    ananas_dic = ananas_data(ananas_out_dir, pdb, oriented_struct)
    if ananas_dic['avg_rmsd'] == None or ananas_dic['avg_rmsd'] < RmsdThreshold:
        return ananas_dic['symmetries']
    else:
        return [{'order': None, 'group': None, 'topology': None, 'repeats': None, 'symm_units': None}]
        
    

def read_transfer_symmetry(cesymm_out_dir,pdb):
    transfer_symm = []
    cesymm_dic = read_analysis_files(cesymm_out_dir, pdb)
    if cesymm_dic['topology'] == 'na':
        transfer_symm.append({'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None})
    else:
        topology = cesymm_dic['topology'][:-1].split(';')
        repeats = cesymm_dic['repeats'][:-1].split(';')
        repeats = cesymm_dic['repeats'][:-1].split(';')
        closed_opened = cesymm_dic['closed_opened'][:-1].split(';')
        symmetry_type = cesymm_dic['symmetry_type'][:-1].split(';')
        repeats_list = []


        for i in range(0,len(topology)):
            reps = repeats[i].strip(')/(/\n')
            reps = reps.split(')(') # eg: reps = ['C_2-280,M_2-280,Q_2-280' , 'D_-1-190,N_-1-190,R_-1-190']
            repeats_list.append([[] for z in range(len(reps[0].split(',')))]) # repeats_list = [[[],[],[]]]
            order = len(repeats_list[0])
            if closed_opened[i]=='CLOSED':
                group = 'Cyclic'
            else:
                group = 'Repeated'
            for pair in reps:
                entry = pair.split(',') # eg: entry = ['C_2-280','M_2-280','Q_2-280']
                for ind, k in enumerate(entry):
                    ch = k.split('_')[0][0]
                    beg = k.split('_')[1]
                    first_digit = beg[0] # in case the residue has a negative id (i.e. -1)
                    end = beg[1:].split('-')[1] # eg: end = 280
                    beg = beg[1:].split('-')[0]
                    beg = first_digit + beg # eg: beg = 2
                    repeats_list[i][ind].append((ch,int(beg),int(end)))

            transfer_symm.append({'order':order, 'topology':topology[i], 'group':group, 'symm_units':symmetry_type[i], 'repeats':repeats_list[i]})

    return transfer_symm

def writeout_benchmark_results(benchmark_dic, out_name, bs, zScoreThreshold, quatsymm_file=''):
    out = open(out_name, 'w')
    out.write('# TN = True negative, FN = false negative, complete = method detected the entire range of the symmetric repeats, partial = method detected only part of the repeats\n')
    bs.seek(0)
    if quatsymm_file != '':
        quatsymm_data = pkl.load(open(quatsymm_file, 'rb'))
    for line in bs:
        fields = line.split('\t')
        inputf = fields[2]
        pdb = fields[2][0:4]
        if pdb == 'PDB':
            out.write(line.strip() + '\t' + 'CE-Symm' + '\t' + 'EncoMPASS'+ '\t' + 'SymD' + '-' + str(zScoreThreshold) + '\t' + 'AnAnaS' + '\t' + 'QuatSymm' + '\n')
            col_names = {}
            for i in range(len(fields)):
                col_names[fields[i]] = i
        elif pdb != 'PDB' and pdb != '' and inputf in benchmark_dic:
            line2 = line.strip()
            repeats = fields[col_names['Repeats']].split(';')  # eg: repeats = ['(C_2-280,M_2-280,Q_2-280)(D_-1-190,N_-1-190,R_-1-190)','(C_1-90,C_102-191)']
            cesymm = []
            encompass = []
            symd = []
            ananas = []
            quatsymm = []
            for rep in repeats:
                for b in benchmark_dic[inputf]:
                    if b['repeats_text'] == rep:
                        cesymm.append(b['repeat_coverage']['CE-Symm'])
                        encompass.append(b['repeat_coverage']['EncoMPASS'])
                        symd.append(b['repeat_coverage']['SymD'])
                        if 'AnAnaS' in b['repeat_coverage']:
                            ananas.append(b['repeat_coverage']['AnAnaS'])
                if quatsymm_file != '':
                    for b in quatsymm_data[inputf]:
                        if b['repeats_text'] == rep or (rep == 'None' and b['repeats_text'][0] == rep):
                            quatsymm.append(b['quatsymm'])

            if len(cesymm)==0:
                cesymm.append('TN')
            if len(symd)==0:
                symd.append('TN')
            if len(encompass)==0:
                encompass.append('TN')
            if len(ananas)==0 and quatsymm_file != '':
                ananas.append('TN')
            elif len(ananas)==0 and quatsymm_file == '':
                ananas.append('-')
                quatsymm.append('-')
            line2 = line2 + '\t' + ';'.join(cesymm) + '\t' + ';'.join(encompass)+ '\t' + ';'.join(symd) + '\t' + ';'.join(ananas) + '\t' + ';'.join(quatsymm) + '\n'
            out.write(line2)
    out.close()
    return

def in_benchmark_counter(symm_dic, inputf, b, method, submethod):
    """
    Deduce whether a given benchmark symmetry has been recognized by a symmetry detection method
    :param symm_dic: a dictionary that contains all the symmetry-related information (benchmark or method-derived) for a given structure
    :param inputf: the name of the structure (whole pdb or chain)
    :param b: a dictionary with the properties of the benchmarked symmetry
    :param method: a symmetry detection method, eg. CE-Symm/EncoMPASS/SymD/etc.
    :param submethod: 'analysis' or 'transfer' depending on which EncoMPASS submethod was used; '' if FN
    :return: whether the benchmarked symmetry was identified by the method (positive) and
    whether the repeats were fully described (coverage)
    """
    def overlap_gap_calc(b_rep, cs_rep, overlap):
        # Determine what the difference is in coverage between two repeat descriptions
        diff = b_rep.difference(cs_rep)
        overlap = overlap + len(b_rep.intersection(cs_rep))
        gap = 0
        if len(diff) > 0:
            prev = list(diff)[0]
        for d in diff:
            if d - prev > 5:
                gap = 0
            else:
                gap += 1
            if gap > 20:
                break
            prev = d
        return overlap, gap

    if submethod == '':
        method_dic = symm_dic[method]
    else:
        method_dic = symm_dic[method][submethod]

    positive = 'no'
    coverage = 'FN' # FN = false negative
    order_exception = 0

    for k, c in enumerate(method_dic):
        cs_positive = 'no'
        cs_coverage = 'FN'

        if b['symm_units'] == c['symm_units'] and b['topology'] == c['topology'] and b['order'] == c[
            'order']:
            if b['order'] == None:
                positive = 'yes'
                coverage = 'complete'
                cs_positive = 'yes'
                cs_coverage = 'complete'
            else:
                for ind, bracket in enumerate(b['repeats']):  # bracket=[('D', 43, 140), ('D', 259, 362)]
                    for entry in bracket:  # entry=('D', 43, 140)
                        gap = 0
                        b_rep = set(range(entry[1], entry[2] + 1))
                        overlap = 0
                        size = len(b_rep)
                        chain_accounted_for = 0
                        if method == 'SymD': # SymD does not report repeat ranges; instead it says which residues were aligned and which were not
                            cs_entry = c['aligned'][0]
                            if entry[0] == cs_entry[0]:
                                chain_accounted_for = 1
                                cs_rep = cs_entry[1]
                                overlap, gap = overlap_gap_calc(b_rep, cs_rep, overlap)
                            if chain_accounted_for == 0:
                                gap += len(b_rep)
                        else:
                            for rep in range(len(c['repeats'])): # other than SymD, other method report repeat ranges
                                for cs_entry in c['repeats'][rep]:
                                    if entry[0] == cs_entry[0]:
                                        chain_accounted_for = 1
                                        cs_rep = set(range(cs_entry[1], cs_entry[2] + 1))
                                        ### subfunction
                                        overlap, gap = overlap_gap_calc(b_rep, cs_rep, overlap)
                            if chain_accounted_for == 0:
                                gap += len(b_rep)

                        if overlap > 0:
                            positive = 'yes'
                            cs_positive = 'yes'
                            if gap > 20 and size - overlap > 20:
                                coverage = 'partial'
                                cs_coverage = 'partial'
                            else:
                                if coverage != 'partial':
                                    coverage = 'complete'
                                    cs_coverage = 'complete'

        elif b['symm_units'] == 'Quaternary' and b['symm_units'] == c['symm_units'] and b['topology'] == c['topology'] \
                and int(b['order']) % int(c['order']) == 0:
            order_exception += 1
        if cs_positive == 'yes' or (order_exception > 0 and c['symm_units'] == 'Quaternary' and 'in_benchmark' not in c):
            c['in_benchmark'] = cs_coverage

    if order_exception > 2:
        print("Warning! Order exception used for this structure and repeat: ", inputf, b['repeats']) # This allows cases in which CE-Symm has broken down a quaternary symmetry into multiple symmetries that in effect describe the same thing to be counted as true; for example, breaking down a C4 into multiple C2s
        positive = 'yes'
        coverage = 'complete'

    return positive, coverage

def count_internal_symmetry(locations, bs, tm_archive, analysis_chain_file, out, omit_type, transfer, pdb_suff, tab):
    print("####### Counting internal symmetry excluding " + omit_type + " proteins. Inferred symmetries considered: " + transfer + " ########")
    out.write("####### Counting internal symmetry excluding " + omit_type + " proteins. Inferred symmetries considered: " + transfer + " ########\n")

    # Define counters
    counters = {}
    template_counters = {'Correct': [], 'Completely Correct': [], 'True Positives': [], 'False Positives': [],
                         'True Negatives': [], 'False Negatives': []}
    methods = ['CE-Symm', 'EncoMPASS', 'SymD']
    for method in methods:
        counters[method] = copy.deepcopy(template_counters)
    counters['Benchmark'] = {'Unique Symmetries': []}

    symm={}

    bs.seek(0)
    benchmark_dic = {}
    for line in bs:
        fields = line.split('\t')
        inputf = fields[2]
        pdb = fields[2][0:4]
        if pdb == 'PDB':
            col_names = {}
            for i in range(len(fields)):
                col_names[fields[i]]=i
        elif pdb != 'PDB' and pdb != '' and fields[col_names['TM-Chains']] != '':
            if len(tm_archive[pdb]['tmchains']) > 0:
                tm_chains_list = ';'.join(sorted(tm_archive[pdb]['tmchains']))
                prot_type = tm_archive[pdb]['class']
            else:
                print("Warning! No TM chains for %s" % pdb)
                tm_chains_list = ''
            symm[inputf] = {}

            if prot_type == omit_type:
                continue
            if tm_chains_list != fields[col_names['TM-Chains']]:
                print("Warning! PDB " + pdb + ": possible change in the biological unit. " +  tm_chains_list + ' != ' + fields[col_names['TM-Chains']])
            if all([True for x in fields[col_names['TM-Chains']] if x in tm_chains_list]) == False:
                print("Warning! PDB "+pdb+" will be skipped: difference between TM chains in database and in set. " + tm_chains_list + ' != ' + fields[col_names['TM-Chains']])
                continue

            # Read the symmetry data in the benchmark set
            print('pdb: ',inputf)
            symm[inputf]['benchmark_set'] = benchmark_set_symmetry(fields, col_names)
            benchmark_dic[inputf] = symm[inputf]['benchmark_set']
            if symm[inputf]['benchmark_set'][0]['order']==None:
               symm[inputf]['benchmark_set'][0]['chains']='-'
            symm[inputf]['general'] = {'str_unique_chains': fields[col_names['Structurally-Unique-Chains']], 'xray_title': fields[col_names['Structure-Title']], 'resolution': fields[col_names['Resolution']],
                                    'review_abbr': fields[col_names['Fold-Abbreviation']], 'review_name': fields[col_names['Fold-Name']]}

            # Read the symmetry data from CE-Symm
            order_matters='no' # Does the order of the chains in the inputf matter?
            symm[inputf]['CE-Symm']=[]
            if len(tm_chains_list)==1:
                symm[inputf]['CE-Symm'] = cesymm_symmetry(pdb, tm_chains_list, order_matters, locations, pdb_suff)
            else:
                ch = inputf[5:]
                chains = cesymm_symmetry(inputf, ch, order_matters, locations, pdb_suff)
                # Add the chain to the dictionary eliminating any redundant "None" entries
                if len(symm[inputf]['CE-Symm'])>0 and symm[inputf]['CE-Symm'][0]['order']==None:
                    symm[inputf]['CE-Symm']=[]
                if chains[0]['order']!=None:
                    for c in symm[inputf]['CE-Symm']:  # Eliminate repeated chain symmetries due to whole complex runs
                        if c['symm_units']=='Internal' and c['repeats'][0][0][0]==ch:
                            symm[inputf]['CE-Symm'].remove(c)
                    symm[inputf]['CE-Symm'] = symm[inputf]['CE-Symm'] + chains
            if len(symm[inputf]['CE-Symm']) == 0:
                symm[inputf]['CE-Symm'] = [{'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None}]

            # Read the symmetry data from SymD
            symm[inputf]['SymD'] = []
            if len(tm_chains_list) == 1:
                symm[inputf]['SymD'] = symd_symmetry(pdb, tm_chains_list, order_matters, locations, pdb_suff,
                                                     zScoreThreshold)
            else:
                ch = inputf[5:]
                chains = symd_symmetry(inputf, ch, order_matters, locations, pdb_suff, zScoreThreshold)
                # Add the chain to the dictionary eliminating any redundant "None" entries
                if len(symm[inputf]['SymD']) > 0 and symm[inputf]['SymD'][0]['order'] == None:
                    symm[inputf]['SymD'] = []
                if chains[0]['order'] != None:
                    for c in symm[inputf]['SymD']:  # Eliminate repeated chain symmetries due to whole complex runs
                        if c['symm_units'] == 'Internal' and c['repeats'][0][0][0] == ch:
                            symm[inputf]['SymD'].remove(c)
                    symm[inputf]['SymD'] = symm[inputf]['SymD'] + chains
            if len(symm[inputf]['SymD']) == 0:
                symm[inputf]['SymD'] = [
                    {'order': None, 'topology': None, 'group': None, 'symm_units': None, 'aligned': None}]

            # Read the symmetry data from EncoMPASS
            ### Symmetry in complex
            whole_inputf = [{'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None}]
            trans_chains = []
            order_matters='yes'
            if transfer == 'yes':
                for f in glob.glob(locations['FSYSPATH']['transfer'] + pdb + '_transfer.axes'):
                    fn = f.split('/')[-1][:-5]
                    trans_chains = read_transfer_symmetry(locations['FSYSPATH']['transfer'], fn)
            ### Symmetry in chain
            chains = []
            if len(fields[col_names['TM-Chains']].strip(';'))>0:
                ch = inputf[5:]
                print(ch)
                ### Symmetry from transfer (lower priority than analysis symmetry)
                trans_chain = []
                if transfer == 'yes':
                    for f in glob.glob(locations['FSYSPATH']['transfer']+inputf+'*transfer.axes'):
                        fn = f.split('/')[-1][:-5]
                        trans_chain = trans_chain + read_transfer_symmetry(locations['FSYSPATH']['transfer'],fn)
                trans_chains = trans_chains + trans_chain
                ### Symmetry from analysis
                chain=[]
                analysis_chain_file.seek(0)
                for a in analysis_chain_file:
                    if a.startswith(inputf):
                        chain = []
                        entry = a.split()
                        filename = entry[16].split('&amp;')
                        filepath = entry[17].split('&amp;')
                        for f in range(len(filename)):
                            chain = chain + cesymm_symmetry(filename[f], ch, order_matters, locations, pdb_suff)
                        break
                chains = chains + chain
            symm[inputf]['EncoMPASS']={}
            if (whole_inputf[0]['symm_units'] == 'Internal' and len(chains)>0 and chains[0]['symm_units'] == 'Internal') or (len(chains)>0 and whole_inputf[0]['order']==None):
                symm[inputf]['EncoMPASS']['analysis'] = chains
            else:
                if whole_inputf[0]['order']==None and len(trans_chains)>0: # used to penalize FP transfer symmetries when no analysis symmetries were detected
                    symm[inputf]['EncoMPASS']['analysis'] = trans_chains
                else:
                    symm[inputf]['EncoMPASS']['analysis']= whole_inputf + chains
            symm[inputf]['EncoMPASS']['transfer'] = trans_chains



            # Record which of the benchmark symmetries were found by each method
            for j, b in enumerate(symm[inputf]['benchmark_set']):
                counters['Benchmark']['Unique Symmetries'].append(inputf)
                b['detected'], b['repeat_coverage'] = {}, {}
                for k, method in enumerate(methods):
                    if method == 'EncoMPASS':
                        positive, coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = 'analysis')
                        tmp_positive, tmp_coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = 'transfer')
                        if tmp_coverage == 'complete' or positive == 'no':
                            positive = tmp_positive
                            coverage = tmp_coverage
                    else:
                        positive, coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = '')


                    b['detected'][method] = positive
                    b['repeat_coverage'][method] = coverage
                    if b['symm_units']!='Quaternary' and (b['symm_units']==None or b['repeats'][0][0][0] in symm[inputf]['general']['str_unique_chains']):
                        if positive == 'yes':
                            counters[method]['Correct'].append(inputf)
                            if coverage == 'complete':
                                counters[method]['Completely Correct'].append(inputf)
                        elif b['symm_units'] != None:
                            counters[method]['False Negatives'].append(inputf)


            for k, method in enumerate(methods):
                order_exception = 0
                chains_with_symmetries = []
                if method == 'EncoMPASS':
                    method_symm_dic = symm[inputf][method]['analysis']
                else:
                    method_symm_dic = symm[inputf][method]
                for k, c in enumerate(method_symm_dic):
                    if method == 'SymD':
                        repeats = c['aligned']
                    else:
                        repeats = c['repeats']
                    if 'in_benchmark' in c and c['order'] != None:
                        if c['in_benchmark'] == 'FN' and c['symm_units'] == 'Quaternary':
                            if 3 > order_exception > 1:
                                counters[method]['True Positives'].append(inputf)
                            order_exception += 1
                        elif c['symm_units'] == 'Quaternary' or repeats[0][0][0] in symm[inputf]['general']['str_unique_chains']:
                            counters[method]['True Positives'].append(inputf)
                            if c['symm_units'] == 'Internal':
                                chains_with_symmetries.append(repeats[0][0][0])
                    elif 'in_benchmark' not in c and c['symm_units']=='Internal' and repeats[0][0][0] in symm[inputf]['general']['str_unique_chains']:
                        counters[method]['False Positives'].append(inputf)
                    elif 'in_benchmark' in c and c['order'] == None:
                        if 'transfer' in symm[inputf][method] and len(symm[inputf][method]['transfer']) == 0:
                            counters[method]['True Negatives'].append(inputf)
                        elif not 'transfer' in symm[inputf][method]:
                            counters[method]['True Negatives'].append(inputf)

                if method == 'EncoMPASS':
                    for k, c in enumerate(symm[inputf][method]['transfer']):
                        if c['repeats'][0][0][0] not in chains_with_symmetries and 'in_benchmark' in c and c['repeats'][0][0][
                            0] in symm[inputf]['general']['str_unique_chains']:
                            counters[method]['True Positives'].append(inputf)
                # examples to monitor: 4or2_A


            print('benchmark_set: ', symm[inputf]['benchmark_set'])
            print('cesymm: ', symm[inputf]['CE-Symm'])
            print('analysis: ', symm[inputf]['EncoMPASS']['analysis'])
            print('transfer: ', symm[inputf]['EncoMPASS']['transfer'], '\n\n')
            print('Number of benchmark symmetries: ', str(len(counters['Benchmark']['Unique Symmetries'])), '\n')
            for method in methods:
                if method == 'SymD':
                    print('SymD z-score Threshold: ', str(zScoreThreshold), '\n')
                print(method, 'Correct:', str(len(counters[method]['Correct'])), '\n')
                print(method, 'Completely Correct:', str(len(counters[method]['Completely Correct'])), '\n')
                print(method, 'True Positives:', str(len(counters[method]['True Positives'])), '\n')
                print(method, 'False Positives:', str(len(counters[method]['False Positives'])), '\n')
                print(method, 'True Negatives:', str(len(counters[method]['True Negatives'])), '\n')
                print(method, 'False Negatives:', str(len(counters[method]['False Negatives'])),','.join(counters[method]['False Negatives']), '\n')

    out.write('Number of benchmark symmetries: ' + str(len(counters['Benchmark']['Unique Symmetries'])) + '\n')
    for method in methods:
        if method == 'SymD':
            out.write('SymD z-score Threshold: ' + str(zScoreThreshold) + '\n')
        out.write(method + ' Correct: ' + str(len(counters[method]['Correct'])) + '\n')
        out.write(method + ' Completely Correct: ' + str(len(counters[method]['Completely Correct'])) + '\n')
        out.write(method + ' True Positives: ' + str(len(counters[method]['True Positives'])) + '\n')
        out.write(method + ' False Positives: ' + str(len(counters[method]['False Positives'])) + '\n')
        out.write(method + ' True Negatives: ' + str(len(counters[method]['True Negatives'])) + '\n')
        out.write(method + ' False Negatives: ' + str(len(counters[method]['False Negatives'])) + '\n')
    out.write('\n')

    print('Total number of benchmark symmetries', len(counters['Benchmark']['Unique Symmetries']))
    for method in methods:
        for p in sorted(counters[method]):
            print(method, p, len(counters[method][p]))


    if omit_type == 'beta':
        type_name = 'alpha'
    else:
        type_name = 'beta'
    if transfer == 'yes':
        bs_out_name = out_dir + 'MemSTATS-' + tab + '_' + type_name + '_internal-inferred_zscore' + str(zScoreThreshold) + '_' + str(datetime.date.today()) + '_results.tsv'
    else:
        bs_out_name = out_dir + 'MemSTATS-' + tab + '_' + type_name + '_internal_zscore' + str(zScoreThreshold) + '_' + str(datetime.date.today()) + '_results.tsv'
    writeout_benchmark_results(benchmark_dic, bs_out_name, bs, zScoreThreshold)
    return counters


def count_quat_symmetry(locations, bs, tm_archive, analysis_whole_file, out, omit_type,pdb_suff, tab):
    print("####### Counting quaternary symmetry excluding " + omit_type + " proteins")
    out.write("####### Counting quaternary symmetry excluding " + omit_type + " proteins" + "\n")

    # Define counters
    counters = {}
    template_counters = {'Correct': [], 'Completely Correct': [], 'True Positives': [], 'False Positives': [],
                         'True Negatives': [], 'False Negatives': []}
    methods = ['CE-Symm', 'EncoMPASS', 'SymD', 'AnAnaS']
    for method in methods:
        counters[method] = copy.deepcopy(template_counters)
    counters['Benchmark'] = {'Unique Symmetries': []}

    symm = {}

    bs.seek(0)
    benchmark_dic ={}
    for line in bs:
        fields = line.split('\t')
        inputf = fields[2]
        pdb = fields[2][0:4]
        if pdb == 'PDB':
            col_names = {}
            for i in range(len(fields)):
                col_names[fields[i]] = i
        elif pdb != 'PDB' and pdb != '' and fields[col_names['TM-Chains']] != '':
            if len(tm_archive[pdb]['tmchains']) > 0:
                tm_chains_list = ';'.join(sorted(tm_archive[inputf]['tmchains']))
                prot_type = tm_archive[inputf]['class']
            else:
                print("Warning! No TM chains for %s" % pdb)
                tm_chains_list = ''
            symm[inputf] = {}

            if prot_type == omit_type:
                continue
            if tm_chains_list != fields[col_names['TM-Chains']]:
                print(
                    "Warning! PDB " + pdb + ": difference between TM chains in database and in set - this pdb will be skipped. ", tm_chains_list, fields[col_names['TM-Chains']])
                continue

            # Read the symmetry data in the benchmark set
            print('pdb: ', inputf)
            symm[inputf]['benchmark_set'] = benchmark_set_symmetry(fields, col_names)
            benchmark_dic[inputf] = symm[inputf]['benchmark_set']
            if symm[inputf]['benchmark_set'][0]['order'] == None:
                symm[inputf]['benchmark_set'][0]['chains'] = '-'

            symm[inputf]['general'] = {'str_unique_chains': fields[col_names['Structurally-Unique-Chains']],
                                       'xray_title': fields[col_names['Structure-Title']],
                                       'resolution': fields[col_names['Resolution']],
                                       'review_abbr': fields[col_names['Fold-Abbreviation']],
                                       'review_name': fields[col_names['Fold-Name']]}

            # Read the symmetry data from CE-Symm
            order_matters = 'no'  # Does the order of the chains in the pdb file matter?
            symm[inputf]['CE-Symm'] = cesymm_symmetry(pdb, tm_chains_list, order_matters, locations, pdb_suff)
            for c in symm[inputf]['CE-Symm']:  # Eliminate internal symmetries
                if c['symm_units'] == 'Internal':
                    symm[inputf]['CE-Symm'].remove(c)
            if len(symm[inputf]['CE-Symm']) == 0:
                symm[inputf]['CE-Symm'] = [
                    {'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None}]
                

            # Read the symmetry data from SymD
            symm[inputf]['SymD'] = symd_symmetry(pdb, tm_chains_list, order_matters, locations, pdb_suff,
                                                 zScoreThreshold)
            for c in symm[inputf]['SymD']:  # Eliminate internal symmetries
                if c['symm_units'] == 'Internal':
                    symm[inputf]['SymD'].remove(c)
            if len(symm[inputf]['SymD']) == 0:
                symm[inputf]['SymD'] = [
                    {'order': None, 'topology': None, 'group': None, 'symm_units': None, 'aligned': None}]
                
            # Read the symmetry data from AnAnaS
            symm[inputf]['AnAnaS'] = ananas_symmetry(pdb, tm_chains_list, locations, pdb_suff)


            # Read the symmetry data from EncoMPASS
            ### Symmetry in complex
            whole_pdb = []
            trans_chains = []
            order_matters = 'yes'
            analysis_whole_file.seek(0)

            for a in analysis_whole_file:
                if a.startswith(pdb):
                    entry = a.split()
                    filename = entry[16]
                    chains_in_order = entry[1]
                    filepath = entry[17]
                    whole_pdb = cesymm_symmetry(filepath, chains_in_order, order_matters, locations, pdb_suff)
                    break
            for c in whole_pdb:  # Eliminate internal symmetries
                if c['symm_units'] == 'Internal':
                    whole_pdb.remove(c)
            if len(whole_pdb) == 0:
                whole_pdb = [{'order': None, 'topology': None, 'group': None, 'symm_units': None, 'repeats': None}]
            symm[pdb]['EncoMPASS'] = {}
            symm[pdb]['EncoMPASS']['analysis'] = whole_pdb
            symm[inputf]['EncoMPASS']['transfer'] = trans_chains

            # Record which of the benchmark symmetries were found by each method
            for j, b in enumerate(symm[inputf]['benchmark_set']):
                b['detected'], b['repeat_coverage'] = {}, {}
                for k, method in enumerate(methods):
                    if method == 'EncoMPASS':
                        positive, coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = 'analysis')
                        tmp_positive, tmp_coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = 'transfer')
                        if tmp_coverage == 'complete' or positive == 'no':
                            positive = tmp_positive
                            coverage = tmp_coverage
                    else:
                        positive, coverage = in_benchmark_counter(symm[inputf], inputf, b, method, submethod = '')


                    b['detected'][method] = positive
                    b['repeat_coverage'][method] = coverage
                    if (b['symm_units'] == None or b['symm_units'] == 'Quaternary'):
                        if k == 0:
                            counters['Benchmark']['Unique Symmetries'].append(inputf)
                        if positive == 'yes':
                            counters[method]['Correct'].append(inputf)
                            if coverage == 'complete':
                                counters[method]['Completely Correct'].append(inputf)
                        elif b['symm_units'] != None:
                            counters[method]['False Negatives'].append(inputf)

            for k, method in enumerate(methods):
                # While detection by transfer does not exist yet for quaternary symmetries in EncoMPASS, the logic has
                # been kept here so that it can be easily added at a later stage.
                order_exception = 0
                chains_with_symmetries = []
                if method == 'EncoMPASS':
                    method_symm_dic = symm[inputf][method]['analysis']
                else:
                    method_symm_dic = symm[inputf][method]
                for k, c in enumerate(method_symm_dic):
                    if 'in_benchmark' in c and c['order'] != None:
                        if c['in_benchmark'] == 'FN' and c['symm_units'] == 'Quaternary':
                            if 3 > order_exception > 1:
                                counters[method]['True Positives'].append(inputf)
                            order_exception += 1
                        elif c['symm_units'] == 'Quaternary' or c['repeats'][0][0][0] in symm[inputf]['general']['str_unique_chains']:
                            counters[method]['True Positives'].append(inputf)
                            if c['symm_units'] == 'Internal':
                                chains_with_symmetries.append(c['repeats'][0][0][0])
                    elif c['order'] != None and 'in_benchmark' not in c and c['symm_units'] != 'Internal':
                        counters[method]['False Positives'].append(inputf)
                    elif 'in_benchmark' in c and c['order'] == None:
                        if 'transfer' in symm[inputf][method] and len(symm[inputf][method]['transfer']) == 0:
                            counters[method]['True Negatives'].append(inputf)
                        elif not 'transfer' in symm[inputf][method]:
                            counters[method]['True Negatives'].append(inputf)

                if method == 'EncoMPASS':
                    for k, c in enumerate(symm[inputf][method]['transfer']):
                        if c['repeats'][0][0][0] not in chains_with_symmetries and 'in_benchmark' in c and c['repeats'][0][0][
                            0] in symm[inputf]['general']['str_unique_chains']:
                            counters[method]['True Positives'].append(inputf)
                            if len(counters[method]['True Negatives'])+len(counters[method]['True Positives'])!=len(counters[method]['Correct']):
                                raise SystemError

                if len(counters[method]['True Negatives']) + len(counters[method]['True Positives']) != len(counters[method]['Correct']):
                    raise SystemError


            print('benchmark_set: ', symm[inputf]['benchmark_set'])
            print('cesymm: ', symm[inputf]['CE-Symm'])
            print('analysis: ', symm[inputf]['EncoMPASS']['analysis'])
            print('transfer: ', symm[inputf]['EncoMPASS']['transfer'], '\n\n')
            print('Number of benchmark symmetries: ', str(len(counters['Benchmark']['Unique Symmetries'])), counters['Benchmark']['Unique Symmetries'], '\n')
            for method in methods:
                if method == 'SymD':
                    print('SymD z-score Threshold: ', str(zScoreThreshold), '\n')
                print(method, 'Correct:', str(len(counters[method]['Correct'])), '\n')
                print(method, 'Completely Correct:', str(len(counters[method]['Completely Correct'])), '\n')
                print(method, 'True Positives:', str(len(counters[method]['True Positives'])), '\n')
                print(method, 'False Positives:', str(len(counters[method]['False Positives'])), '\n')
                print(method, 'True Negatives:', str(len(counters[method]['True Negatives'])),'\n')
                print(method, 'False Negatives:', str(len(counters[method]['False Negatives'])), '\n')


    out.write('Number of benchmark symmetries: ' + str(len(counters['Benchmark']['Unique Symmetries'])) + ': ' + ', '.join(counters['Benchmark']['Unique Symmetries']) + '\n')
    for method in methods:
        if method == 'SymD':
            out.write('SymD z-score Threshold: ' + str(zScoreThreshold) + '\n')
        out.write(method + ' Correct: ' + str(len(counters[method]['Correct'])) + '\n')
        out.write(method + ' Completely Correct: ' + str(len(counters[method]['Completely Correct'])) + '\n')
        out.write(method + ' True Positives: ' + str(len(counters[method]['True Positives'])) + '\n')
        out.write(method + ' False Positives: ' + str(len(counters[method]['False Positives'])) + '\n')
        out.write(method + ' True Negatives: ' + str(len(counters[method]['True Negatives'])) + '\n')
        out.write(method + ' False Negatives: ' + str(len(counters[method]['False Negatives'])) + ': ' +
                  ', '.join(counters[method]['False Negatives']) + '\n')
    out.write('\n')

    print('Total number of benchmark symmetries', len(counters['Benchmark']['Unique Symmetries']))
    for method in methods:
        for p in sorted(counters[method]):
            print(method, p, len(counters[method][p]))

    print(counters)

    if omit_type == 'beta':
        type_name = 'alpha'
    else:
        type_name = 'beta'
    writeout_benchmark_results(benchmark_dic, out_dir + 'MemSTATS-' + tab + '_' + type_name +'_quat_zscore' + str(zScoreThreshold) +'_' + str(datetime.date.today()) + '_results.tsv', bs, zScoreThreshold, locations['FSYSPATH']['quatsymm'])
    return counters

########## Main #################
if __name__ == "__main__":

    # Get the benchmark set
    tab = '20-May-19'
    memstats_file = '../MemSTATS_dataset.xlsx'
    data_xlsx = pd.read_excel(memstats_file, tab + '-internal')
    data_xlsx.to_csv('MemSTATS_dataset_' + tab + '-internal.tsv', encoding='utf-8', index=False, sep='\t')
    data_xlsx = pd.read_excel(memstats_file, tab + '-quaternary')
    data_xlsx.to_csv('MemSTATS_dataset_' + tab + '-quaternary.tsv', encoding='utf-8', index=False, sep='\t')

    # Define paths to key directories and files
    symmetry_dir = 'symmetry/'

    locations = {}
    locations['FSYSPATH'] = {}
    locations['FSYSPATH']['cesymm'] = symmetry_dir + 'cesymm/'  # path to CE-Symm output files
    locations['FSYSPATH']['symd'] = symmetry_dir + 'symd/'  # path to SymD output files
    locations['FSYSPATH']['ananas'] = symmetry_dir + 'ananas/'  # path to AnAnaS output files
    locations['FSYSPATH']['transfer'] = symmetry_dir + 'encompass/transfer/' # path to inferred symmetry output files (EncoMPASS)
    locations['FSYSPATH']['quatsymm'] =  symmetry_dir + 'quatsymm/memstats_with_quatsymm.pkl'  # dictionary of QuatSymm results
    locations['FSYSPATH']['whole'] = '../MemSTATS_pdbs/'  # path to pdb files

    analysis_whole_file = open(symmetry_dir + 'encompass/analysis_symmetries_whole.txt', 'r')
    analysis_chain_file = open(symmetry_dir + 'encompass/analysis_symmetries_chains.txt', 'r')
    tm_archive = pkl.load(open('.tm_archive.pkl', 'rb')) # includes information extracted from inserting proteins in the membrane with OPM
    pdb_suff = '_sb' # coordinate files are named XXXX + pdb_suff + '.pdb', eg. '1okc_sb.pdb'
    out_dir = os.getcwd() + '/results/'
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    locations['OPT'] = {}
    locations['OPT']['temppath'] = out_dir + 'tmp/'  # path to working directory for storing intermediary outputs


    # Define SymD significance threshold. Results with higher z-score will be considered significant.
    zScoreThreshold = 10 # authors recommend 8 or 10

    # Define results file with overall performance statistics
    out_file = open(out_dir + "MemSTATS-" + tab + "_stats_zscore" + str(zScoreThreshold) + '_' + str(datetime.date.today()) + ".txt", "w")

    # Initiate benchmarking
    bs_int = open('MemSTATS_dataset_' + tab + '-internal.tsv', 'r')
    bs_quat = open('MemSTATS_dataset_' + tab + '-quaternary.tsv', 'r')
    results = {}
    types = ['beta', 'alpha']
    for i,omit_type in enumerate(types):
        results[types[(i + 1) % len(types)]] = {}
        results[types[(i + 1) % len(types)]]['quaternary'] = count_quat_symmetry(locations, bs_quat, tm_archive, analysis_whole_file, out_file, omit_type, pdb_suff, tab)
        results[types[(i + 1) % len(types)]]['internal'] = count_internal_symmetry(locations, bs_int, tm_archive, analysis_chain_file, out_file, omit_type, 'no', pdb_suff, tab)
        results[types[(i + 1) % len(types)]]['internal with transfer'] = count_internal_symmetry(locations, bs_int, tm_archive, analysis_chain_file, out_file, omit_type, 'yes',pdb_suff, tab)

    pkl.dump(results, open(out_dir + 'benchmark_results_' + tab + '_zscore' + str(zScoreThreshold) + '_' + str(datetime.date.today()) +'.pkl', 'wb'))


    bs_int.close()
    bs_quat.close()
    analysis_whole_file.close()
    analysis_chain_file.close()
    out_file.close()

