import pytest
import os
import sys
sys.path.append(os.getcwd() + '/')


symmetry_dir = 'symmetry/'
tm_archive_file = '.tm_archive.pkl'
pdb_suff = '_sb'
pdb_list_to_test = ['1okc', '4hea', '2a65']



@pytest.mark.parametrize('pdb,order', [('1okc', None), ('4f35', 2), ('4u1w', 2), ('5aji', 7), ('1k4c', 4)])
def test_ananas_symmetry(test_tm_chains, test_locations_dic, pdb, order):
    from benchmarking_with_MemSTATS import ananas_symmetry
    dic = ananas_symmetry(pdb, test_tm_chains, test_locations_dic, pdb_suff)
    print(dic)
    assert dic[0]['order'] == order


@pytest.mark.parametrize('pdb,order', [('1okc', 3), ('4f35', 2), ('4u1w', None), ('5aji', 7), ('1k4c', None)])
def test_cesymm_symmetry(test_tm_chains, test_locations_dic, pdb, order):
    from benchmarking_with_MemSTATS import cesymm_symmetry
    dic = cesymm_symmetry(pdb, test_tm_chains, 'no', test_locations_dic, pdb_suff)
    print(dic)
    assert dic[0]['order'] == order

@pytest.mark.parametrize('pdb,order', [('1okc', 3), ('4f35', 2), ('4u1w', 2), ('5aji', 7), ('1k4c', 4)])
def test_symd_symmetry(test_tm_chains, test_locations_dic, pdb, order):
    from benchmarking_with_MemSTATS import symd_symmetry
    dic = symd_symmetry(pdb, test_tm_chains, 'no', test_locations_dic, pdb_suff, zScoreThreshold = 8)
    print(dic)
    assert dic[0]['order'] == order