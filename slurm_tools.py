from cmaes_pre_post import *
#from amber_interface import *
import multiprocessing
import subprocess
import time
import numpy
import sys
import os
import re

#----------------------------------------------------------------------
def send_jobs_bal(Iterround, population, prev_population, daddres, dihedtot, dNindep,atoms,sigmas,scaling_ndx,Ngrps,dopttype):

	"""
	Sends jobs to be simulated, each individual as separate job. Takes the input configuration from the last frame of the individual closest to the individual being submitted.
	Input: iterround = optimization round number, population = set of individuals to be converted into simulation inputs, prev_population = set of individuals at the previous round. Used to extract the starting configuration. Daddress = index conversion between dihedral indexing in the individual and the itp. dihedtot = str info read from the ancestor. DNindep = number of independently optimized dihedrals, atoms = atom/charge info, sigmas = LJ sigmas in the ancestor. scaling_ndx = scaling group indexing, Ngrps = number of scaling groups, dopttype = dihedral optimization type
	Output: submits the jobs to the cluster and returns an array containing their processids in the cluster.
	NOTE: 
	-this could be properly rewritten using drmaa or other python package to control Slurm function
	-this will need to be rewritten for each spesific cluster

	"""
	 	
	base='packet'
	arg1='-Tr'
	arg2='..'
	itpname='lipid.itp'
	wd = os.getcwd()
	jobids=[]
	Nold=len(prev_population)

	mindist=1000
	minind=500
	print("Sending jobs for round: ",Iterround)
	for i,ind in enumerate(population):
		
		ind_id=str(Iterround)+'_'+str(i)
		#subprocess.call(['pwd'])
		#copies the "packet" template folder new folder for the individual simulation. Then the job submission file and the itps in that folder are modified accordingly
		subprocess.call(['cp',arg1,base,ind_id])		
		os.chdir(ind_id)
		replace='s/Jobname/qs3nodh_'+str(ind_id)+'/g'
		
		subprocess.call(['sed','-i',replace,'prod.job' ])
		write_itp(ind,dihedtot,dNindep,daddres,atoms,sigmas,scaling_ndx,Ngrps,"lipid_anc.itp",itpname,dopttype)
		#finds the starting structure from the previous optimization round
		for j in range(0,Nold):
			d=numpy.linalg.norm(numpy.subtract(ind,prev_population[j]))
			if d<mindist:
				mindist=d	
				minind=j
		if minind==(Nold-1) or Nold==1:

				amber_inputs("topol.top","equilibrated.gro")
				
		else:
				file_path="../"+str(Iterround-1)+"_"+str(minind)+"/frame.gro"
							
				if os.path.exists(file_path):

					subprocess.call(["cp",file_path,"equilibrated.gro"])
					amber_inputs("topol.top","equilibrated.gro")
					
				else:	
					amber_inputs("topol.top","equilibrated.gro")
		#job submission happends here				
		parts=subprocess.check_output(['sbatch','prod.job']).split()
		
		jobids.append(int(parts[3]))
	
		os.chdir(wd)
	
	return jobids	
#----------------------------------------------------------------------		
	
def wait_jobs(Iterround, procesids):
	"""
	"Listens" to the jobs at the cluster based on their process ids and outputs info to stdout
	NOTE:
	- this could be properly rewritten using drmaa or other python package to control slurm functions
	- Will likely have to be redone in some for to port to new clusters
	- extract_squeue uses username (hard_coded below)
	"""
	 
	username=hantila
	Np=len(procesids)
	exitcodes=[0]*Np
	starttime=time.time()
	exitflag=0

	while exitflag==0:
		finished=[]
		# this extracts the queue information based on username. Mine is hard coded in atm
		runstate=extract_squeue(username)
		
		print("Round: ", Iterround)
		print("Running:")
		print("--------------------------------------")
		print("job-ID\tJob name\tstatus")
		
		for pid in procesids:
			
			p_state = filter(lambda x: x['JOBID'] == str(pid), runstate)
			
			if len(p_state)==0:
				finished.append(pid)
			else:	
				print(str(pid)+" "+p_state[0]['NAME']+"\t"+p_state[0]['ST'])
		print("\n")		
		print("Finished: "+str(len(finished))+" out of "+str(Np))
		print("--------------------------------------")
		print("job-ID \t Job name \t exit_status \t exit_code")
	
		for pid in finished:
			p_state=extract_sacct(pid)
			print(pid, p_state[2],"\t", p_state[0],"\t", p_state[1])
			try:
				indx=re.match('.*?([0-9]+)$',p_state[2]).group(1)
			except AttributeError:
				print("Can't retrieve individual number for job " + str(pid)+"\n")
				continue	
			if p_state[0]=='TIMEOUT':
				exitcodes[int(indx)]==0
			#elif:
				#here crash detection based on smth else than exit codes
			else:
				exitcodes[int(indx)]==int(p_state[1][0])

		if len(finished)==Np:
			exitflag=1
			continue
		#if (time.time()-starttime)/(60*60)>60:
		#	exitflag=2
		#	continue
		sys.stdout.flush()
		
		time.sleep(1*60)

		print("\n\n")	
	return exitcodes
#------------------------------------------------------------------------			

def extract_squeue(username):
	"""
	extracts queue information based on username.
	"""
	pattern = re.compile(r"^\s+(?P<JOBID>\d+)\s+(?P<PARTITION>\w+)\s+(?P<NAME>\w+)\s+(?P<USER>\w+)\s+(?P<ST>\w+)\s+(?P<TIME>[\d:]+)\s+(?P<NODES>\d+)\s+(?P<NODELIST>[\d\w\(\)\-:,\[\]]*)$", re.M)
	
	lines = subprocess.check_output(['squeue','-u',username])
	fields = lines[:lines.find("\n")].split()
	
	records=pattern.findall(lines)
	for j in range(len(records)):
		records[j] = list(records[j])

	runinfo=[dict(zip(fields,records[i])) for i in range(len(records))]
	return runinfo        
    

#---------------------------------------------------------------------------

def extract_sacct(pid):
	"""
	Extracts the the exit status of jobs based on their id using sacct
	"""
	try:
		lines = subprocess.check_output(['sacct','-j',str(pid),"--format=JobId,JobName%20,State,Exitcode"], stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError:
		status="not found"
		exitcode="0:0"
		jobname="not found"
		return [status, exitcode, jobname]
	
	pattern="^"+str(pid)+"\s+(?P<NAME>[\S]+)\s+(?P<State>[\w\+]+)\s+(?P<ExitCode>[\d:]+)\s*$"
	info=re.findall(pattern, lines, re.MULTILINE)
	
	if len(info)==0:
		status="infolength zero"
		exitcode="0:0"
		jobname="not found"
		return [status, exitcode, jobname]
		
	status=info[0][1]
	exitcode=info[0][2]
	jobname=info[0][0]
	return [status, exitcode, jobname]
