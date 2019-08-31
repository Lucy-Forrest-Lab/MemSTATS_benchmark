import pytest
import os
import sys
import tempfile
sys.path.append(os.getcwd() + '/')


symmetry_dir = 'symmetry/'
tm_archive_file = '.tm_archive.pkl'
pdb_suff = '_sb'
pdb_list_to_test = ['1okc', '4hea', '2a65']

@pytest.mark.unit
def test_file_system_case_sensitive():
    path = os.getcwd()
    with tempfile.NamedTemporaryFile(prefix='TmP', dir=path) as tmp_file:
        case_sensitive = (not os.path.exists(tmp_file.name.lower()))
    assert case_sensitive == True


@pytest.mark.data
@pytest.mark.parametrize('pdb,order', [('1okc', None), ('4f35', 2), ('4u1w', 2), ('5aji', 7), ('1k4c', 4)])
def test_ananas_symmetry(test_tm_chains, test_locations_dic, pdb, order):
    from benchmarking_with_MemSTATS import ananas_symmetry
    dic = ananas_symmetry(pdb, test_tm_chains, test_locations_dic, pdb_suff)
    print(dic)
    assert dic[0]['order'] == order


@pytest.mark.data
@pytest.mark.parametrize('pdb,order', [('1okc', 3), ('4f35', 2), ('4u1w', 2), ('5aji', 7), ('1k4c', 4)])
def test_symd_symmetry(test_tm_chains, test_locations_dic, pdb, order):
    from benchmarking_with_MemSTATS import symd_symmetry
    dic = symd_symmetry(pdb, test_tm_chains, 'no', test_locations_dic, pdb_suff, zScoreThreshold = 8)
    assert dic[0]['order'] == order


@pytest.mark.data
@pytest.mark.parametrize('pdb,dir,topology', [('4dx5_analysis', 'symmetry/encompass/selected/', 'Parallel;'),
                                              ('1kqf_I_df_mlev', 'symmetry/encompass/cesymm_minlen/', 'Parallel;')])
def test_read_analysis_file(pdb, dir, topology):
    from benchmarking_with_MemSTATS import read_analysis_files
    dic = read_analysis_files(dir, pdb)
    print(dic)
    assert dic['topology'] == topology


@pytest.mark.data
@pytest.mark.parametrize('pdb,fname,dir,order', [('1okc','1okc','symmetry/cesymm/', 3), ('4f35','4f35','symmetry/cesymm/', 2),
                                           ('4u1w','4u1w','symmetry/cesymm/', None), ('5aji','5aji','symmetry/cesymm/', 7),
                                           ('1k4c','1k4c','symmetry/cesymm/', None),
                                           ('4dx5', '4dx5_analysis', 'symmetry/encompass/selected/', 3),
                                                 ('1j4n', '1j4n_A', 'symmetry/cesymm/', 2)])
def test_cesymm_symmetry(test_tm_chains, test_locations_dic, pdb, fname, dir, order):
    from benchmarking_with_MemSTATS import cesymm_symmetry
    dic = cesymm_symmetry(dir, fname, test_tm_chains, 'no', test_locations_dic, pdb_suff)
    print(dic)
    assert dic[0]['order'] == order


@pytest.mark.data
@pytest.mark.parametrize('pdb,dir,order', [('2nwx_B_12-255_transfer', 'symmetry/encompass/transfer/', 2)])
def test_read_transfer_symmetry(pdb, dir, order):
    from benchmarking_with_MemSTATS import read_transfer_symmetry
    dic = read_transfer_symmetry(dir, pdb)
    print(dic)
    assert dic[0]['order'] == order