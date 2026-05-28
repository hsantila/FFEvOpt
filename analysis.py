import numpy as np
import os 
import subprocess 
from cmaes_pre_post import *
from OrderParameter import *
from amber_interface import *
import random
import time
starttime=time.time()
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
mapping_fname="../mappingPOPCslipids.txt"
centered_fname="centered.xtc"
	
#contains order parameter definitions in order : headgroup, sn1, sn2 	
inp_Sname="../Slipids_POPC.def"
		
fweight=1
fit=0		
fitnessf=open('fitness.txt','w')
subprocess.call(['bash','../form_preprocess.sh', mapping_fname, top_fname,traj_fname,centered_fname])
lines = subprocess.check_output(['bash','../analyse_area.sh']).splitlines()

lines=[float(i) for i in lines[1:]]

area=np.mean(lines)
area=area/17.0

if math.sqrt((area-0.64)**2)<0.04:
	fit=0
else:
	fit=fit+15*(area-0.64)**2

fitnessf.write("fitness from area "+str(fit)+'\n')

OrdParam=find_OP(inp_Sname, top_fname, centered_fname)
	
if len(OrdParam)!=len(exp_op):
	print "Number of experimental order parameters",len(exp_op), "does not match the simulated ones", len(OrdParam)
	
	sys.exit(1)
for i, op in enumerate(OrdParam.values()):
        	
	(op.avg, op.std, op.stem) =op.get_avg_std_stem_OP
	#consider weighting with error here
	fitnessf.write(op.name+" "+str(exp_op[i])+" "+str(op.avg)+" "+str(op.stem)+" "+str(op.avg-exp_op[i])+"\n")
	
	if math.sqrt(abs(op.avg-exp_op[i])**2)<0.01:
		continue				
	fit=fit+10*abs(op.avg-exp_op[i])**2
		
fitnessf.write("final fitness "+str(fit)+'\n')
fitnessf.close()
print "analysis took"+str(starttime-time.time()) 
