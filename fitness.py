import numpy as np
import os 
import subprocess 
from cmaes_pre_post import *
from OrderParameter import *
#from amber_interface import *
import random
import time
def evaluate(myid,Fdata,exp_op, q_min):
	"""
	Old standalone analysis of fitness based on the trajectory which was done within the optimization code. Left here as an example.	
	"""

	#decend into the individual subfolder
	wd = os.getcwd()
	if not os.path.isdir(myid):
		return 100+random.random()*15,
	os.chdir(myid)
	Fsimufile=open('Form_Factor_From_Simulation.dat','w')
	OPsimufile=open('OPs_From_Simulation.dat','w')
	#converting amber trajectory back
	amber2grotraj()
	
	if not os.path.exists("prod.trr"):
		os.chdir(wd)
		print("prod.trr was not produced for individual "+str(myid)) 
		return 100+random.random()*15,
	

		
	traj_fname="prod.trr"
	top_fname="prod.tpr"
	mapping_fname="../mappingPOPCcharmm.txt"
	centered_fname="centered.xtc"
	
	#contains order parameter definitions in order : headgroup, sn1, sn2 	
	inp_Sname="CHARMM_SimulationPOPC.def"
		
	fweight=1
	fit=0	
	
	
	print("Evaluating individual"+str(myid))
	fitnessf=open('fitness.txt','w')
	subprocess.call(['bash','../form_preprocess.sh', mapping_fname, top_fname,traj_fname,centered_fname])
	#subprocess.call([['bash','order_preprocess.sh']])

	# Calculate form factor part of the fitness evaluation
	
	lines = subprocess.check_output(['bash','../form_calc.sh']).splitlines()
	Nform=len(lines)
	Fsimu=np.zeros((Nform, 2))
	
	i=0
	for line in lines:
		Fsimufile.write(line+' \n')
		parts=line.split()
		if float(parts[0])<q_min:
			continue
		else:
			Fsimu[i][0]=float(parts[0])
			Fsimu[i][1]=float(parts[1])
			i=i+1
			
	(maxs, mins) = find_extrema(Fsimu)
	hrat=[]
	

	hrat.append(Fsimu[maxs[0]][1]/Fsimu[maxs[1]][1])
	hrat.append(Fsimu[maxs[1]][1]/Fsimu[maxs[2]][1])			


	fit=fit+fweight*(abs(hrat[0]-Fdata[0]) + abs(hrat[1]-Fdata[1]) + abs(Fsimu[maxs[0]][0]-Fdata[2])+abs(Fsimu[maxs[1]][0]-Fdata[3])+abs(Fsimu[mins[0]][0]-Fdata[4])+abs(Fsimu[mins[1]][0]-Fdata[5]))
	fitnessf.write("fitness from FF "+str(fit)+'\n')
	OrdParam=find_OP(inp_Sname, top_fname, centered_fname)
	if len(OrdParam)!=len(exp_op):
		print("Number of experimental order parameters "+str(len(exp_op))+" does not match the simulated ones" +str(len(OrdParam)))
	
		sys.exit(1)
	for i, op in enumerate(OrdParam.values()):
        	
		(op.avg, op.std, op.stem) =op.get_avg_std_stem_OP
		#consider weighting with error here
		fitnessf.write(op.name+" "+str(exp_op[i])+" "+str(op.avg)+" "+str(op.stem)+" "+str(op.avg-exp_op[i])+"\n")
		fit=fit+abs(op.avg-exp_op[i])
			
	
	fitnessf.write("final fitness "+str(fit)+'\n')
	fitnessf.close()
	#return to main folder
	os.chdir(wd)
				            
	return fit,       	
	
#-----------------------------------------------------------------------------	
def standalone_analysis():
	"""
	Old standalone analysis of fitness based on the trajectory. Left here as an example.
	"""
	starttime=time.time()
	#load and process experimental data
	Fsimufile=open('Form_Factor_From_Simulation.dat','w')
	OPsimufile=open('OPs_From_Simulation.dat','w')
	fexp=read_expdata("../Form_Factor_From_Experiments.dat")
	(maxs, mins) = find_extrema(fexp)
	Fdata=[]

	Fdata.append(fexp[maxs[0]][1]/fexp[maxs[1]][1])
	Fdata.append(fexp[maxs[1]][1]/fexp[maxs[2]][1])	
	Fdata.append(fexp[maxs[0]][0])
	Fdata.append(fexp[maxs[1]][0])
	Fdata.append(fexp[mins[0]][0])
	Fdata.append(fexp[mins[1]][0])

	q_min=fexp[0][0]
	
	exp_op=read_OP("../OP_From_Experiments.dat")
	
	amber2grotraj()
	traj_fname="prod.trr"
	top_fname="prod.tpr"
	mapping_fname="../mappingPOPCcharmm.txt"
	centered_fname="centered.xtc"
	
	#contains order parameter definitions in order : headgroup, sn1, sn2 	
	inp_Sname="CHARMM_SimulationPOPC.def"
		
	fweight=1
	fit=0		
	fitnessf=open('fitness.txt','w')
	subprocess.call(['bash','../form_preprocess.sh', mapping_fname, top_fname,traj_fname,centered_fname])
	lines = subprocess.check_output(['bash','../form_calc.sh']).splitlines()
	Nform=len(lines)
	Fsimu=np.zeros((Nform, 2))
	
	i=0
	for line in lines:
		Fsimufile.write(line+' \n')
		parts=line.split()
		if float(parts[0])<q_min:
			continue
		else:
			Fsimu[i][0]=float(parts[0])
			Fsimu[i][1]=float(parts[1])
			i=i+1
			
	(maxs, mins) = find_extrema(Fsimu)
	hrat=[]
	

	hrat.append(Fsimu[maxs[0]][1]/Fsimu[maxs[1]][1])
	hrat.append(Fsimu[maxs[1]][1]/Fsimu[maxs[2]][1])			


	fit=fit+fweight*(abs(hrat[0]-Fdata[0]) + abs(hrat[1]-Fdata[1]) + abs(Fsimu[maxs[0]][0]-Fdata[2])+abs(Fsimu[maxs[1]][0]-Fdata[3])+abs(Fsimu[mins[0]][0]-Fdata[4])+abs(Fsimu[mins[1]][0]-Fdata[5]))
	fitnessf.write("fitness from FF "+str(fit)+'\n')
	OrdParam=find_OP(inp_Sname, top_fname, centered_fname)
	
	if len(OrdParam)!=len(exp_op):
		print("Number of experimental order parameters "+str(len(exp_op))+" does not match the simulated ones "+str(len(OrdParam)))
	
		sys.exit(1)
	for i, op in enumerate(OrdParam.values()):
        	
		(op.avg, op.std, op.stem) =op.get_avg_std_stem_OP
		#consider weighting with error here
		fitnessf.write(op.name+" "+str(exp_op[i])+" "+str(op.avg)+" "+str(op.stem)+" "+str(op.avg-exp_op[i])+"\n")
		fit=fit+abs(op.avg-exp_op[i])
			
	
	fitnessf.write("final fitness "+str(fit)+'\n')
	fitnessf.close()
	print(str(time.time()-starttime))
#-----------------------------------------------------------------------------
def lean_evaluate(myid):
	"""
	The evaluation function takes the jobid (myid) and decends into the folder where the simulation for that individual is. It then reads the fitness file therein. The fitness in the fitnessfile is pre-calculated at the end of the simulation based on the analysis of the trajectory, hence the optimization code only reads it. If the fitness file is faulty/missing, this function assumes the simulation/individual was faulty and assignes a random, large number for the fitness to penalize. 
	"""
	wd = os.getcwd()
	
	if not os.path.isdir(myid):
		return 100+random.random()*15,
	os.chdir(myid)
		
	print("Evaluating individual"+str(myid))
	if not os.path.exists("prod.trr"):
		os.chdir(wd)
		print("prod.trr was not produced for individual "+str(myid)) 
		return 100+random.random()*15,
		
	if not os.path.exists("fitness.txt"):
		os.chdir(wd)
		print("fitness.txt was not produced for individual "+str(myid)) 
		return 100+random.random()*15,	
		
	
	fo = open("fitness.txt", "r")
	lines=fo.readlines()
	
	
	if len(lines)==0:
		os.chdir(wd)
		print("fitness.txt was not complete "+str(myid)) 
		return 100+random.random()*15,
		
	parts=lines[-1].split()			
	if len(parts)!=3 or parts[0]!="final":
		os.chdir(wd)
		print("fitness.txt was not complete "+str(myid)) 
		return 100+random.random()*15,	
	
	fit=float(parts[2])
	os.chdir(wd)
				            
	return fit,  
