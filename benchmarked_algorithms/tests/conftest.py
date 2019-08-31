import pytest
import os
import sys
import pickle as pkl
sys.path.append(os.getcwd() + '/')

symmetry_dir = 'symmetry/'
tm_archive_file = '.tm_archive.pkl'
pdb_suff = '_sb'
pdb_list_to_test = ['1okc', '4hea', '2a65']

@pytest.fixture
def test_locations_dic():
    locations = {}
    locations['FSYSPATH'] = {}
    locations['FSYSPATH']['cesymm'] = symmetry_dir + 'cesymm/'  # path to CE-Symm output files
    locations['FSYSPATH']['symd'] = symmetry_dir + 'symd/'  # path to SymD output files
    locations['FSYSPATH']['ananas'] = symmetry_dir + 'ananas/'  # path to AnAnaS output files
    locations['FSYSPATH']['transfer'] = symmetry_dir + 'encompass/transfer/' # path to inferred symmetry output files (EncoMPASS)
    locations['FSYSPATH']['quatsymm'] =  symmetry_dir + 'quatsymm/memstats_with_quatsymm.pkl'  # dictionary of QuatSymm results
    locations['FSYSPATH']['whole'] = '../MemSTATS_pdbs/'  # path to pdb files
    locations['OPT'] = {}
    locations['OPT']['temppath'] = os.getcwd() + 'tmptest/'  # path to working directory for storing intermediary outputs
    if os.path.isdir(locations['OPT']['temppath']) == False:
        os.mkdir(locations['OPT']['temppath'])
    locations['out_dir'] = os.getcwd() + 'tmptest/'

    return locations

@pytest.fixture
def test_tm_archive():
    tm_archive = pkl.load(open(tm_archive_file, 'rb'))
    assert '1okc' in tm_archive
    return tm_archive


@pytest.mark.parametrize("pdb", pdb_list_to_test)
@pytest.fixture
def test_tm_chains(test_tm_archive, pdb):
    tm_chains_list = ';'.join(sorted(test_tm_archive[pdb]['tmchains']))
    return tm_chains_list
