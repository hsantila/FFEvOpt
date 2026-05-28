# FFEvOpt
Evolutionary algorithm used to optimize a force field to NMR order parameters and area per lipid

The optimization altgorithm is cma-es (wikipedia article here https://en.wikipedia.org/wiki/CMA-ES) based on the deep python package (https://deap.readthedocs.io/en/master/). Algorithm figure attached for a reminder.

At the moment you can either change the diherals (angle, barrier or just barrier) or charges/LJ sigmas. To go around the restrictions imposed by atom typing, the itp file where these edits happen in an optimization round has to be formatted in a spesific way. The original force field is currently read in from lipid_anc.itp

The atoms for which parameters are optimized are read in from a file called opt_file.txt

For optimizing charges, atoms have to be grouped in neutral charge groups to retain the correct total charge. What is then tuned by the algorithm is the scaling for the charges/LJ sigmas in that charge group. For dihedrals there are options to retain symmetry, randomize angles, zero barries etc which are unfortunately hard coded.

The code forms a checkpoint file once per optimization cycle from which it can be restarted.

For every round one has to run an MD run for the force field candidates (individuals) being tested. It is advisable to optimize sampling vs computational time for this part, aka exactly the type of run you wanna run. The input files for MD are stored in a folder called packet. This is then copied and modified for each individual to a folder named optimizationround_individualnumber. The starting configuration is selected from the individual in the previous round closest to the new candidate.

The analysis of target values from the MD is run in parallel (within the job script sent to the cluster) and the results is stored to the individual's folder. This is then read by the optimization algorithm after all individuals at that round have finished running. The algorithm uses this information to rank the force fields and generates a new set of candidates.

The python package slurm_tools.py sends the jobs for the individuals, and queries their status from the cluster, allowing the optimization only to proceed once the jobs have run.

The file fitness.py contains the script calculating the fitness of the individual. Currently this targets the area per lipid and order parameters but you can customize it to your own needs.

There are quite some things, filenames etc, hard coded for now.


By far the messiest part of the code is the itp handling (cmaes_pre_post.py) and the handling of crashes (in the main run file optimization.py). We can go over this in person.
