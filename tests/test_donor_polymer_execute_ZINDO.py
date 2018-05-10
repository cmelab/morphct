from morphct.definitions import TEST_ROOT
from morphct.code import helper_functions as hf
from morphct import run_MorphCT
from testing_tools import TestCommand
import os
import shutil
import sys
import pytest

@pytest.fixture(scope='module')
def run_simulation():
    # ---==============================================---
    # ---======== Directory and File Structure ========---
    # ---==============================================---

    input_morph_dir = TEST_ROOT + '/assets/donor_polymer'
    output_morph_dir = TEST_ROOT + '/output_EZ'
    input_device_dir = TEST_ROOT + '/assets/donor_polymer'
    output_device_dir = TEST_ROOT + '/output_EZ'

    # ---==============================================---
    # ---========== Input Morphology Details ==========---
    # ---==============================================---

    morphology = 'donor_polymer.xml'
    input_sigma = 3.0
    device_morphology = None
    device_components = {
    }
    overwrite_current_data = True
    random_seed_override = 929292929

    # ---==============================================---
    # ---============= Execution Modules ==============---
    # ---==============================================---

    execute_fine_graining = False                 # Requires: None
    execute_molecular_dynamics = False            # Requires: fine_graining
    execute_obtain_chromophores = False           # Requires: Atomistic morphology, or molecular_dynamics
    execute_zindo = True                         # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = False  # Requires: execute_zindo
    execute_calculate_mobility = False            # Requires: calculate_transfer_integrals
    execute_device_simulation = False              # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---============ Chromophore Parameters ==========---
    # ---==============================================---
    molecule_terminating_connections = {
        'C1': [[2, 1]],
        'C10': [[2, 1]]
    }
    remove_orca_inputs = True
    remove_orca_outputs = True

    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [i for i in dir() if (not i.startswith('_')) and (not i.startswith('@'))\
                      and (i not in ['run_MorphCT', 'helper_functions', 'hf', 'os', 'shutil', 'TestCommand',
                                    'TEST_ROOT', 'setup_module', 'teardown_module', 'testing_tools', 'sys',
                                    'pytest', 'request'])]
    parameters = {}
    for name in parameter_names:
        parameters[name] = locals()[name]


    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---
    try:
        shutil.rmtree(output_morph_dir)
    except OSError:
        pass
    os.makedirs(os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code'))
    shutil.copy(os.path.join(TEST_ROOT, 'assets', os.path.splitext(morphology)[0], 'OC',
                             morphology.replace('.xml', '_post_obtain_chromophores_voronoi.pickle')),
                os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code',
                             morphology.replace('.xml', '.pickle')))
    run_MorphCT.simulation(**parameters)  # Execute MorphCT using these simulation parameters
    # The output dictionary from this fixing
    fix_dict = {}
    # Load the output pickle
    output_pickle_data = hf.load_pickle(os.path.join(output_morph_dir, os.path.splitext(morphology)[0],
                                                          'code', morphology.replace('.xml', '.pickle')))
    fix_dict['output_AA_morphology_dict'] = output_pickle_data[0]
    fix_dict['output_CG_morphology_dict'] = output_pickle_data[1]
    fix_dict['output_CG_to_AAID_master'] = output_pickle_data[2]
    fix_dict['output_parameter_dict'] = output_pickle_data[3]
    fix_dict['output_chromophore_list'] = output_pickle_data[4]
    # Load the correct expected pickle
    expected_pickle_data = hf.load_pickle(os.path.join(input_morph_dir, 'EZ',
                                                       morphology.replace('.xml', '_post_execute_ZINDO.pickle')))
    fix_dict['expected_AA_morphology_dict'] = expected_pickle_data[0]
    fix_dict['expected_CG_morphology_dict'] = expected_pickle_data[1]
    fix_dict['expected_CG_to_AAID_master'] = expected_pickle_data[2]
    fix_dict['expected_parameter_dict'] = expected_pickle_data[3]
    fix_dict['expected_chromophore_list'] = expected_pickle_data[4]
    return fix_dict


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_AA_morphology_dict(self, run_simulation):
        for key in run_simulation['expected_AA_morphology_dict']:
            try:
                self.compare_equal(len(run_simulation['output_AA_morphology_dict'][key]),
                                   len(run_simulation['expected_AA_morphology_dict'][key]))
            except TypeError:
                self.compare_equal(run_simulation['output_AA_morphology_dict'][key],
                                   run_simulation['expected_AA_morphology_dict'][key])

    def test_check_CG_morphology_dict(self, run_simulation):
        self.compare_equal(run_simulation['output_CG_morphology_dict'],
                           run_simulation['expected_CG_morphology_dict'])

    def test_check_CG_to_AAID_master(self, run_simulation):
        self.compare_equal(run_simulation['output_CG_to_AAID_master'],
                           run_simulation['expected_CG_to_AAID_master'])

    def test_check_parameter_dict(self, run_simulation):
        # Pop the system-dependent keys, such as the input and output dirs since this will
        # always be system-dependent
        output_pars = {}
        expected_pars = {}
        for key in run_simulation['expected_parameter_dict']:
            if key  in ['parameter_file', 'output_morph_dir', 'CG_to_template_dirs', 'output_morphology_directory',
                        'input_device_dir', 'input_morphology_file', 'output_device_dir', 'input_morph_dir']:
                continue
            output_pars = run_simulation['output_parameter_dict'][key]
            expected_pars = run_simulation['expected_parameter_dict'][key]
        self.compare_equal(output_pars, expected_pars)

    def test_check_chromophore_list(self, run_simulation):
        for index, chromophore in enumerate(run_simulation['expected_chromophore_list']):
            for key in chromophore.__dict__.keys():
                self.compare_equal(run_simulation['output_chromophore_list'][index].__dict__[key],
                                   chromophore.__dict__[key])

    def test_check_input_single_orca_files_created(self, run_simulation):
        input_morph_dir = run_simulation['output_parameter_dict']['input_morph_dir']
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        chromo_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                                  'input_orca', 'single')
        asset_dir = os.path.join(input_morph_dir, 'EZ', 'input_orca', 'single')
        for file_name in os.listdir(asset_dir):
            self.confirm_file_exists(os.path.join(chromo_dir, file_name))

    def test_check_input_pair_orca_files_created(self, run_simulation):
        input_morph_dir = run_simulation['output_parameter_dict']['input_morph_dir']
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        chromo_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                                  'input_orca', 'pair')
        asset_dir = os.path.join(input_morph_dir, 'EZ', 'input_orca', 'pair')
        for file_name in os.listdir(asset_dir):
            self.confirm_file_exists(os.path.join(chromo_dir, file_name))

    def test_check_output_single_orca_files_created(self, run_simulation):
        input_morph_dir = run_simulation['output_parameter_dict']['input_morph_dir']
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        chromo_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                                  'output_orca', 'single')
        asset_dir = os.path.join(input_morph_dir, 'EZ', 'output_orca', 'single')
        for file_name in os.listdir(asset_dir):
            self.confirm_file_exists(os.path.join(chromo_dir, file_name))

    def test_check_output_pair_orca_files_created(self, run_simulation):
        input_morph_dir = run_simulation['output_parameter_dict']['input_morph_dir']
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        chromo_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                                  'output_orca', 'pair')
        asset_dir = os.path.join(input_morph_dir, 'EZ', 'output_orca', 'pair')
        for file_name in os.listdir(asset_dir):
            self.confirm_file_exists(os.path.join(chromo_dir, file_name))


# TODO: Parse the output files to check we get the same energy levels (should be deterministic)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, 'output_EZ'))


if __name__ == "__main__":
    run_simulation()