# Distributed_DNA_Multi-Strand_Simulator

A software tool simulates DNA multi-strands on a distributed model using NFsim simulator 

by Michael Sneddon http://michaelsneddon.net/nfsim/.  


Program developed by Frankie Spencer,

Project led by Dr. Eugen Czeizler,


At: Algorithms and Software in Bio-Engineering (AlgoBio), 

    Dept. of Information Technologies, 
    
    Åbo Akademi University, 
    
    Åbo, Finland
    

Runs only on Windows 7 or later versions


Requirements

Python 3.7.0 installed
joblib==1.1.0

Active perl 64 - https://www.activestate.com/products/perl/downloads/

NFsim_v1.11 - http://michaelsneddon.net/nfsim/download/


Run instructions:

Please fill the “simulation_parameters.json” according to following instructions. 
Then run main.py. There will be number of command line windows will pop up (number will depend on how many CPU threads are used). 

Note: It may be not possible to use the system during these simulations since they don't run on background. 

--------------------------------------------------------------------------

{

Desired number of CPU threads to be utilized for the simulation, e.g. 8

   "number_of_parallel_threads": 8,

Simulation model time in seconds, e.g. 10

   "simulation_time": 10,

Desired number of test suites to run. Here you can set multiple simulation sets to run, e.g. 10

   "number_of_test_suites": 10,

The input file location on the system. You can find our samples ".species" files at "test_samples"

   "input_species_file": "C:/my_species_file_directory/dx_tile_small-scale_example.species",

Desired directory to collect results. Make sure this directory exists.

   "save_results_directory": "C:/my_save_results_directory/my_test_simulation",

Directory where "perl.exe" located. 

   "perl_interpreter": "C:/my_Perl_directory/Perl64/bin/perl.exe",

Directory where "bng2.pl" located at NFsim simulator package files.

   "nfsim_perl_interface": "C:/my_NFsim_directory/NFsim_v1.11/bng2.pl",

Directory where "NFsim_MSWin32.exe" located at NFsim simulators package files.

   "nfsim_simulator": "C:/my_NFsim_directory/NFsim_v1.11/bin/NFsim_MSWin32.exe",

Delete temporary files for saving system space. We recommend set it to "true" unless you decide otherwise.

   "delete_temporary_files": true
 
}
