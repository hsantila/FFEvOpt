import numpy as np
import subprocess
			

#--------------------------------------------------------
def read_itp(itpname, optfile, symflag):
	"""
	Extracts information for optimization from the ancestor itp file (itpname) and the optimization file (optfile). Determines the symmetry based on symflag.
Returns:
individual = the vector to be optimized, dihedtot = all info read in ffrom dihedrals,dNindep = number of independent dihedrals to be optimized,daddres = conversion between individual indexing and the index of the dihedral in dihedtot, atoms = infor read in for atoms/charges, sigmas = info read in for LJ sigmas, 	scaling_ndx, Ngrps = number of (neutral) scaling groups for charges/sigmas, dopttype = dihedral optimization type (for now, only one type per optimization)
	"""
	objectitp=open(itpname, 'r')

	
	Ndihed=0
	dihedtot=[]
	dihedral_vals=[]
	dihed_ndxs=[]
	atoms={}
	#contains lj sigmas, not to be confused with the optimization sigma
	sigmas={}
	pairtot=[]
	typestot=[]
	atomstot=[]	
	while True:
		line=objectitp.readline()
	
		if not line:
			break	
		if 'dihedrals' in line:
			while True:

				line=objectitp.readline()
	
				if line == '':
					break
				if line == '\n':
					break
				if '[' in line:
					break
				if ';'==line[0]:
					continue
				parts=line.split()
			
				if parts[4]=='2' or parts[4]=='4':
					continue
				dihedtot.append(parts)
				dihedral_vals.append(parts[5:8])
				dihed_ndxs.append(parts[0:4])
				Ndihed=Ndihed+1
				
		if 'pairtypes' in line:
			while True:

				line=objectitp.readline()
	
				if line == '':
					break
				if line == '\n':
					break
				if '[' in line:
					break
				if ';'==line[0]:
					continue
				parts=line.split()
			
				pairtot.append(parts)
		

		if 'atoms' in line:
			while True:

				line=objectitp.readline()
	
				if line == '':
					break
				if line == '\n':
					break
				if '[' in line:
					break
				if ';'==line[0]:
					continue
				parts=line.split()
				atomstot.append(parts)
				atoms[parts[1]]=[int(parts[0]), float(parts[6])] #add sigma here when it is known?

		if 'atomtypes' in line:
			while True:

				line=objectitp.readline()
	
				if line == '':
					break
				if line == '\n':
					break
				if '[' in line:
					break
				if ';'==line[0]:
					continue
				parts=line.split()
				typestot.append(parts)
				sigmas[parts[0]]=[float(parts[5])]
				#note, if it has "special sigmas" for pairtypes section, this finds and reads them in
				if len(parts)>7:
					sigmas[parts[0]]=[float(parts[5]),float(parts[-2])]	

	#retrieves information from the opt file
	dihed_ndx, scaling_ndx, Ngrps,dopttype = read_optfile(optfile)		
	dvals=np.asarray(dihedral_vals)
	dndx=np.asarray(dihed_ndxs)
	
	#here you can choose if dihedrals of same parameters in the ancestor are optimized individually (find_dihed) or optimized as one value so that they will have identical parameters (find_dihed_sym) also at the end.
	if symflag==True:
		individual_tmp, daddres, dNindep=find_dihed_sym(dihed_ndx, dvals, dndx)
	else:
		individual_tmp, daddres, dNindep=find_dihed(dihed_ndx, dvals, dndx)
	individual=[]
	
	#extracts the dihedral parameter optimized based on the optimization type: barrier only or both barrier and equilibrium angle
	if dopttype==1:
		for i in range(0,dNindep):
			individual.append(individual_tmp[2*i+1])

	if dopttype==0:
		individual = list(individual_tmp)
	
	for i in range(0,Ngrps):
		individual.append(0)
		individual.append(0)
	
	

	print(len(individual))
	return individual, dihedtot,dNindep,daddres, atoms, sigmas, scaling_ndx, Ngrps, dopttype
	
#--------------------------------------------------------
def write_itp(individual,dihedtot,dNindep,daddres,atoms,sigmas,scaling_ndx,Ngroup,anc_itpname, indv_itpname, dopttype):	
	"""
	Writes the itp files so that a simulation can be run for each individual in the optimization. 
	Input: individua = vector of the parameters for this individual, dihedtot = parsed str information for the itp. To be replaced with optimized parameters. dNindep = number of independent dihedrals, daddess = conversion between dihedral itp position and the indexing in the individual, atoms = charges of atoms, sigmas = LJ sigmas for atom numbers, scaling_ndx = scaling groups for charges/sigmas, Ngroup = number of scaling groups, anc_itpname = ancestor itp filename, to be read and used as a base for writing the new itps, indv_itpname = name of the itp file to be written for this individual, dopttype = dihedral optimization type.
	Output = None (writes the itp)
	"""
	if dopttype==0:
		for i in range(0,dNindep):
		
			indx=daddres[i]
			for j in indx:
				dihedtot[j][5]=str(individual[2*i]/np.pi*180)	
				dihedtot[j][6]=str(individual[2*i+1])
		noptdih=2*dNindep	
	if dopttype==1:
		for i in range(0,dNindep):
			indx=daddres[i]	
			for j in indx:
				dihedtot[j][6]=str(individual[i])
		noptdih=dNindep

	anc_itp=open(anc_itpname,"r")
	indv_itp=open('tmp.itp', "w")

	lines=indv_itp.readlines

	newline=""	
	while True:

		line=anc_itp.readline()
		
		if not line:
			break	
		if line == '':
			break
		if "dihedrals" in line and dNindep!=0:
			
			indv_itp.write(line)
			i=0	
			Ndihed=len(dihedtot[0])
			newline=""
			while True:
				line=anc_itp.readline()
				parts=line.split()
				
				
				if line == '':
					break
				if ';'==line[0]:
					indv_itp.write(line)
					continue
				if line == ' ':
					indv_itp.write(line)
					
					break
				if line == '\n':
					indv_itp.write(line)
					
					break
				if '[' in line:
					indv_itp.write(line)	
								
					break
				if len(parts)<8:
					indv_itp.write(line)
					break
			
				if parts[4]=='2' or parts[4]=='4':
					indv_itp.write(line)
					continue				
			
				newparts=dihedtot[i]
				for k in range(0,len(newparts)):
					newline=newline+newparts[k]+'\t'
					#print k, parts[k], newline, dihedtot[i]
			
				newline=newline+"\n"
				indv_itp.write(str(newline))	
				i=i+1
				newline=""
			continue
				
		if "atomtypes" in line:
			
			indv_itp.write(line)
				
			while True:
				line=anc_itp.readline()
				parts=line.split()
				
				
				if line == '':
					break
				if ';'==line[0]:
					indv_itp.write(line)
					continue
				if line == ' ':
					indv_itp.write(line)
					
					break
				if line == '\n':
					indv_itp.write(line)
					
					break
				if '[' in line:
					indv_itp.write(line)	
								
					break
		
				if parts[0] not in atoms:
					indv_itp.write(line)
					continue

				type_inx=atoms[parts[0]][0]
			
				if type_inx in scaling_ndx:
					group=scaling_ndx[type_inx]
					newline=parts[0]+" "+parts[1]+" "+parts[2]+" "+parts[3]+" "+parts[4]+" "+str(float(parts[5])*individual[noptdih+Ngroup+group])+" "+parts[6]+"\n"
					indv_itp.write(newline)
					
				else:
					indv_itp.write(line)	
			continue
			
		if "atoms" in line:
			
			indv_itp.write(line)
				
			while True:
				line=anc_itp.readline()
				parts=line.split()
				
				
				if line == '':
					break
				if ';'==line[0]:
					indv_itp.write(line)
					continue
				if line == ' ':
					indv_itp.write(line)
					
					break
				if line == '\n':
					indv_itp.write(line)
					
					break
				if '[' in line:
					indv_itp.write(line)	
								
					break
				if parts[1] not in atoms:
					indv_itp.write(line)
					continue

				type_inx=atoms[parts[1]][0]
				if type_inx in scaling_ndx:
					group=scaling_ndx[type_inx]
					newline=parts[0]+" "+parts[1]+" "+parts[2]+" "+parts[3]+" "+parts[4]+" "+parts[5]+" "+str(float(parts[6])*individual[noptdih+group])+" "+parts[7]+"\n"
					indv_itp.write(newline)
					
				else:
					indv_itp.write(line)	
			continue
#this section writes the pairtypes (special pair interactions) for the optimized sigmas according to the mixing rules. Uses special types of (CHARMM-style) sigmas fo this if they were present in the original itp.
		if "pairtypes" in line:
			
			indv_itp.write(line)
				
			while True:
				line=anc_itp.readline()
				parts=line.split()
				
				
				if line == '':
					break
				if ';'==line[0]:
					indv_itp.write(line)
					continue
				if line == ' ':
					indv_itp.write(line)
					
					break
				if line == '\n':
					indv_itp.write(line)
					
					break
				if '[' in line:
					indv_itp.write(line)	
								
					break
				
				if (parts[0] not in atoms) and (parts[1] not in atoms):
					indv_itp.write(line)
					continue
				if parts[0] not in atoms:
					type_inx1=-1
					type_inx2=atoms[parts[1]][0]	

				elif parts[1] not in atoms:
					type_inx1=atoms[parts[0]][0]
					type_inx2=-1	
				
				else:
					type_inx1=atoms[parts[0]][0]
					type_inx2=atoms[parts[1]][0]		

				if type_inx2 in scaling_ndx and type_inx1 in scaling_ndx:

					group1=scaling_ndx[type_inx1]
					group2=scaling_ndx[type_inx2]
					sigma1=sigmas[parts[0]]
					sigma2=sigmas[parts[1]]
					pair_sigma=0.5*(individual[noptdih+Ngroup+group1]*sigma1[-1]+individual[noptdih+Ngroup+group2]*sigma2[-1])
					newline=parts[0]+" "+parts[1]+" "+parts[2]+" "+str(pair_sigma)+" "+parts[4]+"\n"
					indv_itp.write(newline)	
					continue
										
				elif type_inx2 in scaling_ndx:	

					group=scaling_ndx[type_inx2]
					sigma1=sigmas[parts[0]]
					sigma2=sigmas[parts[1]]
					pair_sigma=0.5*(sigma1[-1]+individual[noptdih+Ngroup+group]*sigma2[-1])
					newline=parts[0]+" "+parts[1]+" "+parts[2]+" "+str(pair_sigma)+" "+parts[4]+"\n"
					indv_itp.write(newline)
					
				elif type_inx1 in scaling_ndx:	
	
					group=scaling_ndx[type_inx1]
					sigma1=sigmas[parts[0]]
					sigma2=sigmas[parts[1]]
					pair_sigma=0.5*(sigma1[-1]*individual[noptdih+Ngroup+group]+sigma2[-1])
					newline=parts[0]+" "+parts[1]+" "+parts[2]+" "+str(pair_sigma)+" "+parts[4]+"\n"
					indv_itp.write(newline)					
					
				else:
					indv_itp.write(line)	
			continue				
		else:
			indv_itp.write(line)		
		


	indv_itp.close()
	anc_itp.close()
	subprocess.call(['cp', "tmp.itp", indv_itpname ])	
	
#--------------------------------------------------------
def find_dihed_sym(dihed_ndx, dvals, dndx):
	"""
	finds the ancestor to be optimized based on the atom indices (dihed_ndx) read in from the opt file. Dvals are the diheral values read in from the ancestor itp.
	Dndx contains gorups of 4 atoms making up the diherals read in from the itp. This function is the version that forces symmetry in the diherals.
	Output: ancestor individual, dihedral parameter to output index conversion, Nindep = number of independent diherals to be optimized.
	"""
	Nindep=0
	individual=[]
	addres=[]
	#print dvals[0:3]
	indxsym={}
	indxd=[]
	
	for i in dihed_ndx:

		indx=np.where(dndx==str(i))
		
		for d in indx[0]:
			indxd.append(d)
	
	for ind in indxd:
		if ind not in indxsym:
			indx=np.where((dvals==(dvals[ind][0],dvals[ind][1],dvals[ind][2])).all(axis=1))
			indxsym[ind]=indx[0]


	for ind in indxd:
		toaddres=[]
		if dvals[ind][0]=='-1':
			continue	
		individual.append(float(dvals[ind][0])/180*np.pi)
		#individual.append(float(dvals[ind][0]))
		individual.append(float(dvals[ind][1]))	
		if ind in indxsym:
	
			for s in indxsym[ind]:
				if s==ind:
					continue	
				if s in indxd:
					toaddres.append(s)
					dvals[s][0]='-1'
		
		dvals[ind][0]='-1'
		toaddres.append(ind)
		addres.append(toaddres)	
		
				
		Nindep=Nindep+1

	return individual, addres, Nindep	
#--------------------------------------------------------
def find_dihed(dihed_ndx, dvals, dndx):
	"""
	finds the ancestor to be optimized based on the atom indices (dihed_ndx) read in from the opt file. Dvals are the diheral values read in from the ancestor itp.
	Dndx contains gorups of 4 atoms making up the diherals read in from the itp. This function does not force symmetry.
	Output: ancestor individual, dihedral parameter to output index conversion, Nindep = number of independent diherals to be optimized.
	"""

	Nindep=0
	individual=[]
	addres=[]
	#print dvals[0:3]	
	for i in dihed_ndx:
		#if barrier is zero, assume no need to optimize

		#if abs(float(dvals[i][1]))<10**-5:
		#	continue


		#raw_input("Press the <ENTER> key to continue...")	
		#indx=np.where((dvals==(dvals[i][0],dvals[i][1],dvals[i][2])).all(axis=1))
		indx=np.where(dndx==str(i))
		indx=indx[0]


		
		
		for ind in indx:
			if dvals[ind][0]=='-1':
				continue	
			individual.append(float(dvals[ind][0])/180*np.pi)
			#individual.append(float(dvals[ind][0]))
			individual.append(float(dvals[ind][1]))	
			addres.append(ind)
			Nindep=Nindep+1
			dvals[ind][0]='-1'
	return individual, addres, Nindep

		
#----------------------------------------------------------------------------	
def read_optfile(optfile):
	"""
	Reads in: optimization file
	Outputs: dihed_ndx = indexes for atoms whose diherals will be optimized, Nqs = number of neutral charge groups to sigma/charge optimization. scaling_ndx = charge group index for each atom. dopttype = diheral optimization type (barier or both barrier and equilibrium angle)
	"""	
	fopt=open(optfile, 'r')
	lines=fopt.readlines()
	Nqs=0
	scaling_ndx={}
	dihed_ndx=[]
	dopttype=0
	for line in lines:
		parts=line.split()

		if line == '':
			break
		if line == '\n':
			break
		if parts[0]=="#":
			continue
		if parts[0]=='yes':		
			Nqs=Nqs+1	
			for i in parts[2:]:
				scaling_ndx[int(i)]=Nqs-1

		if parts[1]=='yes' or parts[1]=='bar':	
			for p in parts[2:]:
				dihed_ndx.append(int(p))
		#now the last line of the opt file determines the dihedral optimization type for all. This should be fixed		
		if parts[1]=='bar':
			dopttype=1
	return dihed_ndx, scaling_ndx, Nqs, dopttype	
		
#-----------------------------------------------------------------------------

def read_expdata(datafile):
	opf=open(datafile,'r')
	lines=opf.readlines()
	data=[]
	for line in lines:
		if '#' in line:
			continue
		if 'label' in line:
			continue
		parts=line.split()
	
		data.append([float(p) for p in parts])
	n=len(data)
	m=len(data[1])
	#data_out=np.empty((n, m))
	data_out=np.array(data)

	return data_out

#---------------------------------------------------------------------------------

def read_OP(datafile):
	opf=open(datafile,'r')
	lines=opf.readlines()
	data=[]
	for line in lines:
		if '#' in line:
			continue
		if 'label' in line:
			continue
		parts=line.split()
	
		data.append(float(parts[1]))
	data_out=np.array(data)

	return data_out	

#--------------------------------------------------------------------------------
def find_extrema(data):
	maxs=[]	
	mins=[]	
	Npoints=len(data)

	data_grad=np.gradient(data[:,1])
	for i in range(0,Npoints-2):

		v1=data_grad[i]
		v2=data_grad[i+1]

		if v1>=0 and v2 < 0:
			maxs.append(i)
		if v1<=0 and v2 > 0:
			mins.append(i)

	return maxs, mins

#--------------------------------------------------------------------------------
def traj_preprosess():

	subprocess.call(['bash','form_preprocess.sh'])
	subprocess.call([['bash','order_preprocess.sh']])

#--------------------------------------------------------------------------------
def is_number(s):
#this is obsolete for now
    try:
        float(s)
        return True
    except ValueError:
        return False
#--------------------------------------------------------------------------------
def read_OPdef(op_def):
	op_names=[]	
	OPs={"gamma_C13a",
	"gamma_C13b",
	"gamma_C13c",
	"gamma_C14a",
	"gamma_C14b",
	"gamma_C14c",
	"gamma_C15a",
	"gamma_C15b",
	"gamma_C15c",
	"beta1",
	"beta2",
	"alpha1",
	"alpha2",
	"g3_1",
	"g3_2",
	"g2_1",	
	"g1_1",
	"g1_2"}

	opdeff=open(op_def, 'r')
	lines=opdeff.readlines()
	for line in lines:
		for op in OPs:
			if op in line:
				parts=line.split()	
				op_names.append(parts)		
	return op_names	
#------------------------------------------------------------------------------------	

#---------------------------------------------------------------------------------
def find_op_ind(def_itp, op_names):
	objectitp=open(def_itp, 'r')
	names=[]
	while True:
		line=objectitp.readline()
		
		if 'atoms' in line:
			while True:
				line=objectitp.readline()
					
				if line == '\n':
					break
				if '[' in line:
					break
				parts=line.split()	
				if ';' in parts[0]:
					continue
					
				for i in range(0,len(op_names)):
					
					if op_names[i][2]==parts[4] or op_names[i][3]==parts[4]:
			#			print line, names[i][2], names[i][3]	
						parts=line.split()
						names.append(parts[0])	
			#			print names, parts, op_names[i][2], op_names[i][3]	
						break
										
						
			break

	return names



	
