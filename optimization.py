


"""
Implementation of CMA-ES optimization using the Deap python library by Hanne Antila. The implementation is constructed based on the example codes given on the DEAP github https://github.com/DEAP/deap/tree/82f774d9be6bad4b9d88272ba70ed6f1fca39fcf/examples.
"""

from collections import deque
from cmaes_pre_post import *
from slurm_tools import *
from fitness import *
import multiprocessing
import subprocess
import time
import numpy
import os
import re
import random
import pickle
import random
import argparse


from deap import algorithms
from deap import base
from deap import benchmarks
from deap import cma
from deap import creator
from deap import tools



creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

def parse_args():
    parser=argparse.ArgumentParser()
    
    parser.add_argument("--opt_file", type=str, default="opt_file.txt", help="A file that contains instructions on which parameters to optimize and how")
    parser.add_argument("--anc_itp", type=str, default="lipid_anc.itp", help="A file that contains a gromacs itp file for the force field used as a starting point")
    parser.add_argument("--checkpoint", type=str,default=None ,help="A file that contains a checkpoint from previous round. If not given, writes to a file checkpoint.pkl")
    parser.add_argument("--logbook", type=str, default="logbook.txt",help="A file that contains a checkpoint from previous round. If not given, writes to a file checkpoint.pkl")
    parser.add_argument("--zero", type=str, default="no", help="No or yes to indicate if the optimization should start from zero for the dihedral values being optimized (or random number in case of dihedral angle). Will not affect optimization of charges/LJ sigmas")
    parser.add_argument("--sigma", type=float, default=0.15, help="Sets the sigma---aka the strength of mutation---for the optimization")
    parser.add_argument("--sym", type=bool, default=True, help="In case dihedrals are optimized, this sets the symmetries equal to the force field used as a starting point. Two diherdals identical in the starting point will then be also identical in the optimized force field")
    args=parser.parse_args()
    return args
	
def checkBounds(min_ang, max_ang, min_b, max_b, dNindep, Ngroup, dopttype):
    """
    This function sets boundaries for optimized values: 
    min_ange= minimun dihedral angle, max_angle = max dihedral angle, min_b = minimum barrier height, max_b = max barrier height. dNindep = number of independent dihedrals to be optimized, Ngroup = number of scaling groups (neutral atom groups to be optimized) for sigma and charges. dopttype = type of dihedral optimization.
    """
    def decorator(func):
       def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                 if dNindep!=0 and dopttype==0:
                    for i in xrange(2*dNindep):
                       if i%2==0:
                          if child[i] > max_ang:
                             child[i] = child[i]-max_ang
                          elif child[i] < min_ang:
                             child[i] = child[i]+max_ang
                       else:
                          if child[i] > max_b:
                             child[i] = max_b
                          elif child[i] < min_b:
                             child[i] = min_b	

                 if dNindep!=0 and dopttype==1:
                    for i in xrange(dNindep):
                       if child[i] > max_b:
                          child[i] = max_b
                       elif child[i] < min_b:
                          child[i] = min_b	
			#else: Add boundaries here for sigma and q scaling
			#prevent negative prefactors
			#For now, the last else makes sure that none of the scalings for charges or Lennard-Jones sigmas are negative.
                 else:
                    for i in range(2*dNindep,2*dNindep+2*Ngroup):
                       if child[i]<0:
                          child[i]=0														
            return offspring
       return wrapper
    return decorator

def calc_boundarypenalty(ind,dNindep, dopptype):
#This function penalizes large barriers in the dihedral optimization.  Inputs the individual, number of independent dihedrals optimized and the dihedral optimization type. Penalty is applied beyon the boundary value set by ref for each barrier over ref or under -ref. 
	fit=0
	ref=7
	if dopptype==1:
		for i in xrange(dNindep):
			if ind[i]<-ref or ind[i]>ref:
                                fit=fit+0.05

	if dopptype==0:
		for i in range(0,2*dNindep):
			if ind[2*i+1]<-ref or ind[2*i+1]>ref:
				fit=fit+0.05	
	
	return fit	

def main():

	args = parse_args()
	anc_itp=args.anc_itp
	opt_file=args.opt_file
	checkpoint=args.checkpoint
	logbook_fileN=args.logbook
	zero=args.zero
	symflag=arg.sym
	
	#this extracts information from the ancestor itp file and the optimization file to start the optimization
	ancestor, dihedtot,dNindep,daddres, atoms, sigmas, scaling_ndx, Ngrps, dopttype=read_itp(anc_itp,opt_file,symflag)
	
	
	#initialize logbook file
	lb_file=open(logbook_fileN,'a',0)	
	
	
	
	# Problem size
	N = len(ancestor)
	
	#here modify ancestor if neccesary
	random.seed(0)
	if dopttype==0 and zero=="yes":
		#this sets the angle to random number for dihedral and barrier to zero 
		for j in range(0,dNindep):
			ancestor[j*2]=random.random()*2*numpy.pi
			ancestor[j*2+1]=0
	if dopttype==1 and zero=="yes":
                #this sets the dihedral barrier to zero comment of if not needed
                ancestor=[0]*N

	#this sets the scalings for lj sigmas and charges to 1
	if Ngrps!=0:
                ancestor[-2*Ngrps:]=[1 for i in range(2*Ngrps)]
           
			
	sigma = args.sigma    # 1/5th of the domain is a solid choice usually
	ind_id=0
	lambda_ = 4 + int(3 * numpy.log(N))
	
	toolbox = base.Toolbox()
	#lean_evaluate is the function used to evaluate individuals. You can write your own here. Should take in an individual (np vector) and return a fitness value
	toolbox.register("evaluate", lean_evaluate)
	halloffame = tools.HallOfFame(1)
	stats = tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)

	logbooks = list()
	#this is the variable for iteration step number
	t = 0

	#The stopping criteria are hard coded in
	#this will make the optimization stop at 750, could be implemented as command line input
	MAXITER = 750
	TOLHISTFUN = 10**-12
	TOLHISTFUN_ITER = 10 + int(numpy.ceil(30. * N / lambda_))
	EQUALFUNVALS = 1. / 3.
	EQUALFUNVALS_K = int(numpy.ceil(0.1 + lambda_ / 4.))
	TOLX = 10**-12
	#if factors weighting small barriers are not added the minimun fitness is 0
	TOLFIT = 0.01
	TOLUPSIGMA = 10**20
	CONDITIONCOV = 10**14
	STAGNATION_ITER = int(numpy.ceil(0.2 * t + 120 + 30. * N / lambda_))
	NOEFFECTAXIS_INDEX = t % N
	equalfunvalues = list()
	bestvalues = list()
	medianvalues = list()
	mins = deque(maxlen=TOLHISTFUN_ITER)
	
	#Min and max values for dihedral barriers/angles are hard coded in for now
	MAX_ang = 2*numpy.pi
	MIN_ang = 0
	MAX_b=12
	MIN_b=-12

	pop_tmp=[]
	pop_chkpnt=[]
	

	# We start with the original force field as centroid 
	strategy = cma.Strategy(centroid=ancestor, sigma=sigma, lambda_=lambda_)
	
	
	toolbox.register("generate", strategy.generate, creator.Individual)
	toolbox.decorate("generate", checkBounds(MIN_ang, MAX_ang, MIN_b, MAX_b, dNindep, Ngrps,dopttype))
	toolbox.register("update", strategy.update)

	logbooks.append(tools.Logbook())
	logbooks[-1].header = "gen", "evals","sigma_val", "std", "min", "avg", "max"

	conditions = {"MaxIter" : False, "TolHistFun" : False, "EqualFunVals" : False,
				  "TolX" : False, "TolUpSigma" : False, "Stagnation" : False,
				  "ConditionCov" : False, "NoEffectAxis" : False, "NoEffectCoor" : False, "TolFit":False}	
		  		
	if checkpoint:
		# A file name has been given, then load the data from the file
		with open(checkpoint, "r") as cp_file:
			cp = pickle.load(cp_file)
		population = cp["population"]
		t = cp["generation"]
		dopttype = cp["dopt_type"]
		halloffame = cp["halloffame"]
		bestvalues = cp["best"]
		medianvalues = cp["median"]
		mins = cp["mins"]
		equalfunvalues = cp["eqfval"]
		logbooks = cp["logbook"]

		strategy=cp["strat"]
		toolbox.register("generate", strategy.generate, creator.Individual)
		toolbox.decorate("generate", checkBounds(MIN_ang, MAX_ang, MIN_b, MAX_b, dNindep, Ngrps,dopttype))
		toolbox.register("update", strategy.update)
		for ind_id, ind in enumerate(population):
			pop_chkpnt.append(ind)	
		toolbox.update(population)
	
		t=t+1
	


	## Note that the algorithm won't stop by itself on the optimum (0.0 on rastrigin).
	while not any(conditions.values()):
		

		# Generate a new population
		population = toolbox.generate()
		
		pop_chkpnt.append(ancestor)
		
		#sends jobs to cluster
		#this is likely something you will have to rewrite to work on your own cluster. 
		procesids=send_jobs_bal(t, population,pop_chkpnt, daddres, dihedtot, dNindep,atoms,sigmas,scaling_ndx,Ngrps,dopttype)
		exitcodes=wait_jobs(t, procesids)
		
		# if all individuals exited with nonzero code, something went very wrong
		if numpy.count_nonzero(exitcodes)==lambda_:
			sys.exit("Whole generation exited with nonzero code, breaking the iteration")	

		pop_tmp=[]
		pop_chkpnt=[]
		#once jobs have finished, individuals can be evaluated
		for ind_id, ind in enumerate(population):
		  
			idval=str(t)+'_'+str(ind_id)			
			
			# assume that if the run exited, it crashed and the proposed force field is faulty. Assign random, large fitness value which will push the individual to the end of the sorted population
			if exitcodes[ind_id]!=0:
				ind.fitness.values = (random.random()*15+100,)
				pop_tmp.append(ind)
				print("individual "+str(ind_id)+" exited with non-zero code")
			else:	
				ind.fitness.values = toolbox.evaluate(idval)
				# A boundary penalty could be applied like this
				#ind.fitness.values[0]=ind.fitness.values[0]+calc_boundarypenalty(ind,dNindep,dopttype)
				
				if ind.fitness.values[0]<=100:
					pop_tmp.append(ind)	
			#this is to keep unsorted population for chekpointing and for new simulations using old indexing
			pop_chkpnt.append(ind)			
		l_tmp=len(pop_tmp)
			
		halloffame.update(pop_tmp)
		#statistics calculated for valid individuals only
		record = stats.compile(pop_tmp)
		logbooks[-1].record(gen=t, evals=l_tmp,sigma_val=strategy.sigma, **record)
		
		#write log to file
		lb_file.write(logbooks[-1].stream+"\n")
		
		if l_tmp<=int(strategy.lambda_ / 2):
			print("Faulty individuals included in the production of next generation")
			toolbox.update(population)
			# Count the number of times the k'th best solution is equal to the best solution
			# At this point the population is sorted (happens in above toolbox.update)
			if population[-1].fitness == population[-EQUALFUNVALS_K].fitness:
				equalfunvalues.append(1)
		
			# Log the best and median value of this population
			bestvalues.append(population[-1].fitness.values)
			medianvalues.append(population[int(round(len(population)/2.))].fitness.values)
		else:	
			toolbox.update(pop_tmp)
			# Count the number of times the k'th best solution is equal to the best solution
			# At this point the population is sorted (method update)
			if pop_tmp[-1].fitness == pop_tmp[-EQUALFUNVALS_K].fitness:
				equalfunvalues.append(1)
		
			# Log the best and median value of this population
			bestvalues.append(pop_tmp[-1].fitness.values)
			medianvalues.append(pop_tmp[int(round(len(pop_tmp)/2.))].fitness.values)


		#checkpointing
		if t!=0:
			cp = dict(population=pop_chkpnt, generation=t, halloffame=halloffame,best=bestvalues,median=medianvalues,mins=mins, eqfval=equalfunvalues, logbook=logbooks,strat=strategy, dopt_type=dopttype)
			if checkpoint==None:
				checkpoint="checkpoint.pkl"
			with open(checkpoint, "wb") as cp_file:
				pickle.dump(cp, cp_file)
				# Update the strategy with the evaluated individuals
				
                #clean directories from last round
		prev_dir=str(int(t-1))+'_{0..'+str(int(lambda_-1))+'}'
		command='find '+prev_dir+' -name "*out*" -delete'
		os.system(command)
		command='find '+prev_dir+' -name "*nc" -delete'
		os.system(command)
		command='find '+prev_dir+' -name "*trr" -delete'
		os.system(command)
		command='find '+prev_dir+' -name "*xvg" -delete'
		os.system(command)
		command='find '+prev_dir+' -name "*mdinfo" -delete'
		os.system(command)
		


		t += 1
		STAGNATION_ITER = int(numpy.ceil(0.2 * t + 120 + 30. * N / lambda_))
		NOEFFECTAXIS_INDEX = t % N
		
		

		if t >= MAXITER:
			# The maximum number of iteration per CMA-ES ran
			conditions["MaxIter"] = True
		
		mins.append(record["min"])
		if abs(population[-1].fitness.values[0]) < TOLFIT:
			# The optimization finds good enough minimum
			conditions["TolFit"] = True
		if (len(mins) == mins.maxlen) and max(mins) - min(mins) < TOLHISTFUN:
			# The range of the best values is smaller than the threshold
			conditions["TolHistFun"] = True

		if t > N and sum(equalfunvalues[-N:]) / float(N) > EQUALFUNVALS:
			# In 1/3rd of the last N iterations the best and k'th best solutions are equal
			conditions["EqualFunVals"] = True

		if all(strategy.pc < TOLX) and all(numpy.sqrt(numpy.diag(strategy.C)) < TOLX):
			# All components of pc and sqrt(diag(C)) are smaller than the threshold
			conditions["TolX"] = True
		
		# Need to transfor strategy.diagD[-1]**2 from pyp/numpy.float64 to python
		# float to avoid OverflowError
		if strategy.sigma / sigma > float(strategy.diagD[-1]**2) * TOLUPSIGMA:
			# The sigma ratio is bigger than a threshold
			conditions["TolUpSigma"] = True
		
		if len(bestvalues) > STAGNATION_ITER and len(medianvalues) > STAGNATION_ITER and \
		   numpy.median(bestvalues[-20:]) >= numpy.median(bestvalues[-STAGNATION_ITER:-STAGNATION_ITER + 20]) and \
		   numpy.median(medianvalues[-20:]) >= numpy.median(medianvalues[-STAGNATION_ITER:-STAGNATION_ITER + 20]):
			# Stagnation occured
		   conditions["Stagnation"] = True

		if strategy.cond > 10**14:
			# The condition number is bigger than a threshold
			conditions["ConditionCov"] = True

		if all(strategy.centroid == strategy.centroid + 0.1 * strategy.sigma * strategy.diagD[-NOEFFECTAXIS_INDEX] * strategy.B[-NOEFFECTAXIS_INDEX]):
			# The coordinate axis std is too low
			conditions["NoEffectAxis"] = True

		if any(strategy.centroid == strategy.centroid + 0.2 * strategy.sigma * numpy.diag(strategy.C)):
			# The main axis std has no effect
			conditions["NoEffectCoor"] = True

	stop_causes = [k for k, v in conditions.items() if v]
	print("Stopped because of condition%s %s" % ((":" if len(stop_causes) == 1 else "s:"), ",".join(stop_causes)))
	
	for ind in population:
		print(ind,  ind.fitness.values)
	return halloffame

if __name__ == "__main__":
    main()
