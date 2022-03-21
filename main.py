import os
import math
import sys
from collections import Counter
from random import shuffle
from system_files.convert_results_dump_to_species import convert_dump_to_species
from system_files.complexes_post_processor import complexes_post_process
from datetime import datetime
from shutil import move
from joblib import Parallel, delayed
import json
from pathlib import Path
from system_files.shared_classes import (convert_link_address,
                                         read_file,
                                         write_file,
                                         delete_temp_files)

with open('simulation_parameters.json') as f:
    parameters = json.load(f)

# initialize parameters from simulation_parameters.json file
number_of_parallel_threads = parameters["number_of_parallel_threads"]
simulation_time = parameters["simulation_time"]
number_of_test_suites = parameters["number_of_test_suites"]
input_species_file = parameters["input_species_file"]
save_results_directory = parameters["save_results_directory"]
perl_interpreter = parameters["perl_interpreter"]
nfsim_perl_interface = parameters["nfsim_perl_interface"]
nfsim_simulator = parameters["nfsim_simulator"]
delete_temporary_files = parameters["delete_temporary_files"]

# number of splits to be performed during the simulations
number_of_splits = math.ceil(simulation_time / 0.1)

# k1 (kinetic rate) to be calculated using this value
k1_coefficient = 1.225

# the length of time to perform the simulation
run_time = round((simulation_time / number_of_splits), 6)

# simulation duration for a single step
run_time_per_step = number_of_splits / simulation_time  # only for time stamp purpose

if number_of_splits <= 1:
    raise Exception('Too short simulation time for distributed processing.')


def run_bngl_on_threads():

    # maximum possible number of threads, i.e max splittable
    alternative_n_threads = number_of_parallel_threads
    input_species_file_name = Path(input_species_file).stem
    sum_nucleotides = None

    # formulate NFSim simulation run command
    def get_nfsim_run_command(current_step_folder, step_bngl_file, step_xml_file, result_dump_folder, thread_run_time):

        simulation_command = 'START CMD /C "CD "{}" && ' \
                             '"{}" "{}" -xml "{}" && ' \
                             '"{}" -utl 1000 -xml "{}" -dump "[0:{}:{}]-^>{}/" ' \
                             '-oSteps 1 -sim {}"'.format(current_step_folder,
                                                         convert_link_address(perl_interpreter),
                                                         convert_link_address(nfsim_perl_interface),
                                                         step_bngl_file,
                                                         convert_link_address(nfsim_simulator),
                                                         step_xml_file, thread_run_time,
                                                         thread_run_time,
                                                         result_dump_folder,
                                                         thread_run_time)

        return simulation_command

    global step_session_data, step_session_folder
    step_session_data, step_session_folder = {}, None

    # initialize main session folder name
    session_folder_name = '{}---simulation_results---{}'.format(input_species_file_name,
                                                                datetime.now().strftime('%d-%m-%Y--%H%M%S'))
    # initialize and create main session directory
    session_directory = os.path.join(save_results_directory,
                                     session_folder_name)
    os.mkdir(session_directory)

    # declare essential static bngl syntax
    begin_molecule_syntax = ['begin molecule types\n',
                             '',
                             'end molecule types']

    begin_species_syntax = ['begin species\n',
                            '',
                            'end species']

    # update k1 kinetic value on the bngl file to be used in each simulation step's thread
    bngl_parameters_k1_updated = ['k1 '
                                  + str(round((number_of_parallel_threads * k1_coefficient)
                                              * float(p.split(' ')[1]), 6))
                                  if 'k1 ' in p else p for p in read_file('bngl_script_files/parameters.bngl')]
    bngl_parameters = read_file('bngl_script_files/parameters.bngl')
    bngl_species = read_file('bngl_script_files/species.bngl')
    bngl_observables = read_file('bngl_script_files/observables.bngl')
    bngl_functions = read_file('bngl_script_files/functions.bngl')
    bngl_reaction_rules = read_file('bngl_script_files/reaction_rules.bngl')

    # setup all session variables including bngl and xml files to relevant a dictionaries
    def setup_session_variables(complexes_list_thread, n_runs):
        global step_session_data, step_session_folder

        # current step/round directory of the simulation
        current_step_directory = os.path.join(session_directory, 'step---{}'.format(n_runs))
        os.mkdir(current_step_directory)

        session_dictionary = {}

        # loop through the complexes splits and include them in the bngl and save then in the relevant directory
        for thread, complex, begin_molecule_state in zip(range(1, alternative_n_threads + 1),
                                                         complexes_list_thread['complexes_all_split'],
                                                         complexes_list_thread['begin_molecule_state']):

            thread_data_directory = os.path.join(current_step_directory, 'thread---{}'.format(thread))
            os.mkdir(thread_data_directory)

            thread_bngl_file_name = '{}---{}---{}.bngl'.format(input_species_file_name, n_runs, thread)
            thread_bngl = os.path.join(thread_data_directory, thread_bngl_file_name)
            thread_xml_file_name = '{}---{}---{}.xml'.format(input_species_file_name, n_runs, thread)
            thread_xml = os.path.join(thread_data_directory, thread_xml_file_name)

            begin_molecule_syntax[1] = begin_molecule_state + '\n'
            species_script = [begin_species_syntax[0], *complex, '\n', begin_species_syntax[-1]]

            # formulate bngl file syntax
            bngl_file_syntax = bngl_parameters_k1_updated \
                               + begin_molecule_syntax \
                               + species_script \
                               + bngl_observables \
                               + bngl_functions \
                               + bngl_reaction_rules

            write_file(thread_bngl, bngl_file_syntax)

            session_dictionary.update({thread: {'thread_dir': thread_data_directory,
                                                'thread_bngl': thread_bngl,
                                                'thread_xml': thread_xml}})

        # assign current_step and session_dictionary to global step_session_folder
        step_session_folder = current_step_directory
        step_session_data = session_dictionary

    # fetch all complexes from all simulated threads and attach them together
    def attach_complexes():
        all_complexes_attached = []
        n_nucleotides_fetched = 0

        # read through the dumped files and keep fetching output data until all fetched
        while sum_nucleotides != n_nucleotides_fetched:
            all_complexes_attached = []
            n_nucleotides_fetched = 0
            for thread_dir in step_session_data.values():

                dump_file_links = list(next(os.walk(thread_dir['thread_dir'])))
                dump_file_link = os.path.join(dump_file_links[0],
                                              str([i for i in dump_file_links[2]
                                                   if i.endswith('.0') and not i.endswith('.0.dump.0')][0]))

                species_list = convert_dump_to_species(dump_file_link, '', '', 'read_dump')

                if type(species_list) == list:
                    try:
                        n_nucleotides_fetched += sum([len(i[0].split('.')) * int(i[1])
                                                      for i in [[l.rsplit('  ', 1)[0], l.rsplit('  ', 1)[1]]
                                                                for l in species_list]])
                    except:
                        pass

                    for single_complex in species_list:
                        all_complexes_attached.append(single_complex)

        return all_complexes_attached

    # split complexes to given number of baskets which are to be processed by parallel threads
    def get_split_complexes(all_complexes):
        complex_index, complexes_dictionary, sum_of_nucleotides = [], {}, 0

        # add index to all complexes and count their weight in (number of nucleotides per complex)
        n = 1
        for complex in all_complexes:
            item = [int(complex[1]), len(complex[0].split('.')), n]
            complex_index.append(item)
            complexes_dictionary.update({n: complex[0]})
            sum_of_nucleotides += item[0] * item[1]
            n += 1

        # reinitialize for sorting purpose
        expanded = []
        for c_index in complex_index:
            for _ in range(c_index[0]):
                expanded.append([c_index[2], c_index[1], 1])

        # sort complexes by weight descending order and get n larger ones seperated (one per basket)
        # then initialize the rest of them anf shuffle
        baskets = {}
        sorted_by_comp_size = list(reversed(sorted(expanded, key=lambda x: x[1])))
        biggest_comps_per_n_thread = sorted_by_comp_size[:alternative_n_threads]
        rest_of_comps = sorted_by_comp_size[alternative_n_threads:]

        shuffle(rest_of_comps)

        expanded_rearranged = biggest_comps_per_n_thread + rest_of_comps

        for thread in range(1, alternative_n_threads + 1):
            baskets.update({thread: [[], []]})

        # put available complexes to baskets. put the upcoming complex on a loop to the smallest basket by weight basis
        total_weight = 0
        current_basket = 1
        for complex_e_r in expanded_rearranged:
            baskets[current_basket][0].append(complex_e_r[0])
            baskets[current_basket][1].append(complex_e_r[1])
            total_weight += complex_e_r[1]
            basket_sizes = {i[0]: sum(i[1][1]) for i in baskets.items()}
            smallest_basket = min(basket_sizes, key=basket_sizes.get)
            current_basket = smallest_basket

        # check the minimum possible threads to be utilized on the next split/round
        possible_thread_count = sum([1 for i in baskets.items() if i[1] != [[], []]])

        complexes_for_threads = []
        for b_item in baskets.items():
            if b_item[1] != [[], []]:
                get_count = Counter(b_item[1][0])
                item_set = [[complexes_dictionary[f[0]], f[1]] for f in get_count.items()]
                complexes_for_threads.append(item_set)

        return {'complexes_for_threads': complexes_for_threads,
                'possible_thread_count': possible_thread_count}

    # convert file links to command line acceptable format
    def convert_to_run_formats(thread_data):

        bngl_file = convert_link_address(thread_data[1]['thread_bngl'])
        xml_file = convert_link_address(thread_data[1]['thread_xml'])
        dump_dir = convert_link_address(thread_data[1]['thread_dir'])

        run_d_dic = {'bngl_file': bngl_file,
                     'xml_file': xml_file,
                     'dump_dir': dump_dir}

        return run_d_dic

    # pause program until simulation processes in parallel are finished, return true when done
    def wait_for_process(current_session_folder, n_threads):
        check_items_directory = {}

        # loop through and check if the results dumps are created at each thread's directory
        while sum([i for i in check_items_directory.values()]) != n_threads:
            dump_p = list(os.walk(current_session_folder))[1:]
            dump_f = list(os.walk(current_session_folder))[0][1]

            for d, f in zip(dump_p, dump_f):

                check_dump = len([i for i in d[2] if i.endswith('.0') and not i.endswith('.0.dump.0')])

                if check_dump != 0:
                    if f not in check_items_directory:
                        check_items_directory.update({f: check_dump})

        return True

    # set fg~ state to species complexes
    def complexes_set_fg_state(all_complexes):

        comps_and_bngl_info = {'complexes_all_split': [], 'begin_molecule_state': []}

        for complex in all_complexes:

            comps_with_state, begin_mol_line = [], 'N(b~A~T~C~G,5,3,W,fg'
            for ssdna, c_len in zip(complex, range(0, len(complex))):
                ssdna_with_state = '.'.join(['{},fg~{})'.format(nuc[:-1], c_len)
                                             for nuc in ssdna[0].split('.')]) + '  ' + str(ssdna[1])

                comps_with_state.append(ssdna_with_state)
                begin_mol_line += '~' + str(c_len)

            comps_and_bngl_info['complexes_all_split'].append(comps_with_state)
            comps_and_bngl_info['begin_molecule_state'].append(begin_mol_line + ')')

        return comps_and_bngl_info

    # run simulation by call with command
    def run_simulation(cmd_command):
        os.system(cmd_command)

    # initialize input species as a list
    species_set = [['.'.join([','.join(n.split(',')[:-1]) + ')'
                              for n in l.split('  ')[0].split('.')]), l.split('  ')[1]] for l in
                   read_file(input_species_file) if l.startswith('N')]

    # calculate and save number of nucleotides in the input file
    sum_nucleotides = sum([len(i[0].split('.')) * int(i[1]) for i in species_set])

    pct = 100 / number_of_splits
    print('\rSimulation progress...0%, 0/{} steps completed.'.format(number_of_splits), end="")

    for run_step in range(1, number_of_splits + 1):
        progress_pct = run_step * pct
        step_model_time = run_step / run_time_per_step

        if run_step % run_time_per_step == 0:
            step_model_time = int(step_model_time)

        # split complexes saved at species_set
        split_complexes = get_split_complexes(species_set)

        # give fg~ state
        complexes_state_given = complexes_set_fg_state(split_complexes['complexes_for_threads'])

        # declare number of parallel threads to be utilized
        alternative_n_threads = split_complexes['possible_thread_count']

        setup_session_variables(complexes_state_given, run_step)

        # setup list of command line callable commands for the list of jobs (parallel simulations)
        job_list = []
        for thread in step_session_data.items():
            run_data = convert_to_run_formats(thread)
            nfsim_run_command = get_nfsim_run_command(run_data['dump_dir'], 
                                                      run_data['bngl_file'], 
                                                      run_data['xml_file'], 
                                                      run_data['dump_dir'], run_time)
            job_list.append(nfsim_run_command)

        # run simulations in parallel
        Parallel(n_jobs=alternative_n_threads)(delayed(run_simulation)(inputTuple) for inputTuple in job_list)

        # make a pause until all simulations are done
        wait_for_process(step_session_folder, alternative_n_threads)

        # results to save as species file name
        save_species_file_name = '{}_(step-{})_(threads-{})_nf.{}_step_result.species'.format(input_species_file_name,
                                                                                              run_step,
                                                                                              alternative_n_threads,
                                                                                              round(step_model_time, 6))
        # make a path link to the file to be saved
        save_species_path = os.path.join(step_session_folder, save_species_file_name)

        # get all complexes formed by each thread as a single bunch
        # so this will be saved to a single file as the step's results
        all_complexes_attached = attach_complexes()

        # run post process, which is necessary to reduce identical complexes
        species_set = complexes_post_process(all_complexes_attached, '', '', '')

        # convert post processed complexes list to species standard format
        complexes_list_st_format = [''.join(['# ' + ''.join([c[4] for c in e[0].split('.')]) + ", 5' - 3'\n",
                                             str(e[0] + '  ' + e[1])]) + '\n' for e in species_set]

        # write results file at the relative step directory
        write_file(save_species_path, complexes_list_st_format)

        # put number of parallel threads back to user desired number
        alternative_n_threads = number_of_parallel_threads

        print('\rSimulation progress...{}%, {}/{} steps completed.'.format(round(progress_pct), 
                                                                           run_step, 
                                                                           number_of_splits), end="")

    return session_directory


# run given number of test suites
for _ in range(number_of_test_suites):
    time_start = datetime.now()
    run_and_get_address = run_bngl_on_threads()
    time_end = datetime.now()
    sim_duration = str(time_end - time_start).rsplit('.', 1)[0].replace(':', '.')

    # delete temporary files for saving hard drive space, given by user "delete_temporary_files": true/false
    if delete_temporary_files:
        delete_temp_files(run_and_get_address)

    # rename the session directory with run time duration
    move(run_and_get_address, run_and_get_address + '---' + sim_duration)
