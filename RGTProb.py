from decimal import *
from sys import exit
import math

# Preliminary data acquisition, and preference settings, we produce a list of species trees, a list of gene trees and possibly a dictionary of time lenghts

InputFormat=input('Select the desired input format: type 1 for standard format or 2 for newick format\n')
InputType=input('Select how to submit inputs: type 1 to submit files or 2 to type trees manually\n')

# Checking for valid inputs

if InputFormat not in ['1','2'] or InputType not in ['1','2']:
   exit('Error, invalid input format request!\nPlease restart the program.')

if InputType=='1': # processing files to extract the trees on which to perform the calculations
   fileST=input('Type the name of the file containing SPECIES trees as filename.txt\n')
   openST=open(fileST)
   dataST=[]
   fileGT=input('Type the name of the file containing GENE trees as filename.txt\n')
   openGT=open(fileGT)
   dataGT=[]
   for st in openST:
      st=st.strip()
      dataST.append(st)
   for gt in openGT:
      gt=gt.strip()
      dataGT.append(gt)
      
if InputType=='2': # Allows for the manual insertion of GTs and STs
   if InputFormat=='1':
	   print('Please write trees without blank spaces and characters that are not letters, numbers or parenthesis: an example of correct typing is ((AB)2C)1')
   else:
	   print('Please write trees without blank spaces and characters that are not letters, numbers, parenthesis, commas or colons: an example of correct typing is ((A:0.2,B:0.2):0.2,C:0.4)')
   speciestree=input('Type the SPECIES tree:\n')
   dataST=[speciestree]
   genetree=input('Type the GENE tree:\n')
   dataGT=[genetree]

# Only standard format needs other information to create a dictionary of time lengths, newick format already contains the time information

if InputFormat=='1': 
   IntervalMode=input("Select the desired time interval mode: type 1 for Yule's model, 2 to insert time intervals manually or 3 for a polynomial output (readable by Wolfram Mathematica)\n")
   if IntervalMode=='1':
      SpecRatePar=input('Insert the parameter for the speciation rate\n')
   elif IntervalMode=='2':
      times_raw=input('Write the desired time intervals from T(2) to T(N-1) only separated by a comma\n')
      times=times_raw.split(',')
      for n in range (2,NumSpec):
         timeDict[n]=float(times[j])
   elif IntervalMode!='3':
      exit('ERROR: invalid input. Please restart the program.')

def Time_Yule(par):# function that sets the time intervals according tu Yule's model
   global timeDict
   timeDict=dict()
   for i in range(2,NumSpec):
      timeDict[i]=1/(i*par)
   return timeDict

# Now we have our data in a way that allows us to use the same functions for both cases. Now the desired output format must be chosen

if IntervalMode!='3':
	OutputFormat=input('Choose the desired output format: type 1 for desktop print (not optimal for big trees or a collection of trees), 2 for 1 file for each ST-GT couple or 3 for 1 file for each ST.\n')
else:
	OutputFormat=input('Choose the desired output format: type 2 for 1 file for each ST-GT couple, 3 for 1 file for each ST.\n')
OutputProbab=input('Type 1 if you only want total GT|ST probabilities, 2 if you also want GT&RH|ST probabilities\n')

# This function calculates the vector S_species of the positions of the species names in the input string st, it also calculates the number of species
 
def S_calc(st):
	global S_species,NumSpec
	NumSpec=0
	S_species=[]
	for i in range(len(st)):
		if st[i] not in ['1','2','3','4','5','6','7','8','9','0','.',':',',','(',')']:
			S_species.append(i)
			NumSpec+=1

# We could say that each caracter of the input string represents a part of the tree (parentheses represent arches, number represents internal nodes and
# letters represents leaves). The following function calculates the vector Rd (Root distance) s.t. Rd[i] is the number of nodes separating from the root 
# the part of the tree represented by the i^th caracter of the input string

def T_calc(Tree): 
	n=0
	Rd=[]
	for i in range(len(Tree)):
		if Tree[i]=='(':
			n+=1
		if Tree[i]==')':
			n-=1
		Rd.append(n)
	return Rd
	
# Given a ranked tree, the following two functions calculate the matrix M_ij containing the rank of the "last common ancestor" (from now on LCA) of the
# leaves i and j. This matrix contains all the information about his ranked topology. 
# M is a symmetrical matrix and M_ii=0, so we only need to calculate M_ij where i>j. Our final output will be a list where M[i][j]=M_{i,j+i+1}.
# To calculate M we observe that the part of the string referring to the LCA of A and B is the one with the minimum distance from the root that we can 
# find between the letters A and B. So we can find that part finding the minimum value of the vector Rd between the positions of A and B.

# The following function analyzes species trees provided in newick format, it also extrapolates the dictionary of times. It first calculates the vector
# BranchL s.t. BranchL[i] isn't 0 only if the i^th character of the string is ':', in that case BranchL[i] is the branch lenght written immediately after
# the ':' in the input string.

def newick_analysis(ST):
	global timeDict
	BranchL=[]
	Rd=T_calc(ST)
	for i in range(len(ST)):
		if ST[i] in [':']:
			j=i+1
			string=''
			while ST[j] in ['1','2','3','4','5','6','7','8','9','0','.']:
				string=string+ST[j]
				j=j+1
			BranchL.append(Decimal(string))
		else:
			BranchL.append(0)

# Now we calculate the matrix preM s.t. preM_ij is the sum from T(N) to T(s+1) where s is the rank of the LCA of i and j. To calculate preM_ij we sum
# the lenghts of all the branches separating i from the LCA of i and j

	preM=[]
	times=set() # This will be the set of values appearing in preM
	for i in range(NumSpec-1):
		preM.append([])
		for j in range(i+1,NumSpec):
			n=min(Rd[S_species[i]:S_species[j]])
			h=S_species[i]+1
			t=BranchL[h]
			for k in range(Rd[S_species[i]]-n):
				while Rd[h]!=Rd[S_species[i]]-k-1:
					h+=1
				h+=1
				t=t+BranchL[h]
			preM[i].append(t)
			times.add(t)

# Now we calculate the ranking of the internal nodes by sorting the values if the times in preM, then we produce the matrix M with these ranks. 
# The dictionary of times contains the difference between consecutive values in the ordered list of times.

	tim=list(times)
	tim.sort(reverse=True)
	tim.append(0)
	timeDict=dict()
	M=[]
	for i in range(NumSpec-1):
		M.append([])
		timeDict[i+2]=float(tim[i]-tim[i+1])
		for j in range(NumSpec-1-i):
			h=0
			while tim[h]!=preM[i][j]:
				h+=1
			M[i].append(h+1)
	return M
	
# The following function analyzes species trees provided in standard format and gene trees. The vector S indicates the position of the species names in
# the input string s.t. string[S[i]]=st[S_specie[i]], so if string is the species tree S=S_species, but if string is a gene tree they can be different. 
# The procedure is analogous to the one in the previous function.

def read(string,S): #Processing of a species tree when provided in standard format
   Rd=T_calc(string)
   M=[]
   for i in range(NumSpec-1):
      M.append([])
      for j in range(i+1,NumSpec):
         if S[i]<S[j]:
            n=min(Rd[S[i]:S[j]])
         else: #this is useful to process GTs in which the order of the species is different from the ST
            n=min(Rd[S[j]:S[i]])
         k=max(S[i],S[j])
         while Rd[k]>=n:
            k+=1
         if k+2==len(string):
            M[i].append(1)
         elif string[k+2] in ['1','2','3','4','5','6','7','8','9','0']: #this is useful for big trees with more than 10 leaves
            M[i].append(int(string[k+1:k+3]))
         else:
            M[i].append(int(string[k+1]))
   return M

#Now we can calculate the Maximal Ranked History (MRH) of a GT in a ST, the vector nodes and the matrices topology and Top.

def findaMRH(ST,GT):
	global nodes
	S_gene=[] #list containing the position of letters in GT in a way s.t. GT[S_gene[i]]=ST[S_species[i]]
	for i in range(NumSpec):
		n=0
		while GT[n]!=ST[S_species[i]]:
			n+=1
		S_gene.append(n)
	gene=read(GT,S_gene)
	L=[] #L[i] will contain the ranks of the LCAs (according to the ST) of the couples of species whose LCA in the GT has rank i+1
	nodes=[]
	for i in range(NumSpec-1):
		L.append([])
	for i in range(NumSpec-1):
		for j in range(NumSpec-i-1):
			L[gene[i][j]-1].append(specie[i][j])
		k=0
		while i+1 not in gene[k]:
			k+=1
		nodes.append(k)
	MRH=[L[NumSpec-2][0]]
	for i in range(NumSpec-2):
		MRH.append(min(MRH[i],min(L[NumSpec-3-i])))
	MRH.sort()
	return MRH
   
def topocalc():
	global Top,topology
	topology=[[]] #topology[i][j] is the branch containing the (i+1) phyletic line during the (j+1)^th time interval
	Top=[] #Top[i][j] is the branch in the (i+1)^th time interval containing the phyletic line of the branch j+1 in the (i+2)^th time interval
	for i in range(NumSpec-1):
		topology.append([])
		topology[0].append(1)
		Top.append([])
	for i in range(NumSpec-1):
		n=0
		for j in range(NumSpec-1):
			if specie[i][0]>j:
				topology[i+1].append(topology[i][j])
			else:
				topology[i+1].append(topology[i][j]+1)
			if specie[j][0]==i+1:
				n+=1
				Top[i].append(n)
				Top[i].append(n)
			elif specie[j][0]<i+1:
				n+=1
				Top[i].append(n)

def rhcalc(tup): #input: MRH. Output: a list with all possible rhs. Uses update_rh
   D=[(1,)]
   for j in range(1,len(tup)):
      D=update_rh(D,tup[j])
   return D

def update_rh(ins,Z): #Input: a list of RHs. Output: a list of RHs which are all those obtained by elongating of 1 element the
   d=[]#a list of all ranked histories
   for key in ins:
      for n in range(max(key),Z+1):
         newkey=key+(n,)
         d.append(newkey)
   return d
   
def k_calc(RH): #Output: a dictionary of k[ijz], it is useful to define m[NumSpec]=0 and k(NumSpec,0,z)=1, it also calculates the dictionary of m
	global m
	Pos=[] #a matrix s.t. Pos[i][j] is the branch containing the (j+1)^th coalescence in tau_(i+1)
	for i in range(NumSpec-1):
		Pos.append([])
		Pos[RH[i]-1].append(topology[nodes[i]][RH[i]-1])
	m=dict() #m[i] is the number of nodes of the gene tree in the i-th interval of the species tree
	k=dict()
	for i in range(1,NumSpec+1):
		m[i]=RH.count(i)
		k[(NumSpec,0,i)]=1
	for i in range(1,NumSpec):
		for z in range (NumSpec-i):
			k[(NumSpec-i,m[NumSpec-i],z+1)]=0
		for j in range(m[NumSpec-i]+1):
			if j==0:
				for z in range(NumSpec-i+1):
					k[(NumSpec-i,m[NumSpec-i],Top[NumSpec-1-i][z])]+=k[(NumSpec-i+1,0,z+1)]
			else:
				for z in range(NumSpec-i):
					if z+1==Pos[NumSpec-1-i][m[NumSpec-i]-j]:
						k[(NumSpec-i,m[NumSpec-i]-j,z+1)]=k[(NumSpec-i,m[NumSpec-i]-j+1,z+1)]-1
					else:
						k[(NumSpec-i,m[NumSpec-i]-j,z+1)]=k[(NumSpec-i,m[NumSpec-i]-j+1,z+1)]
	for i in range(NumSpec):
		del(k[(NumSpec,0,i+1)])
	return k
   
def Lambda_calc(RH): #Input: a dictionary of k[ijz]. Output: a dictionary of lambda(i,j). Uses k_calc
   global l
   l=dict()
   k=k_calc(RH)
   for chiave in k:
      ij=chiave[0:2]
      l[ij]=l.get(ij,0)+k[chiave]*(k[chiave]-1)//2

def den(rh,i,j):
	d1=1
	for k in range(m[i]+1):
		if k!=j:
			d1=d1*(l[(i,k)]-l[(i,j)])
	return d1

def probab(ST,GT,RH): #Input: ST, GT and a RH. Output: P(G and RH|ST). Uses Lambda_calc
   Lambda_calc(RH)
   P=dict()
   if IntervalMode=='3':
      P[NumSpec-1]='1'
      for i in range(NumSpec-1,1,-1):
         summ='('
         for j in range(m[i]+1):
            summ+='Exp['+str(-l[(i,j)])+'*T['+str(i)+']]/('+str(den(RH,i,j))+')+'
         summ=summ[0:len(summ)-1]
         summ+=')'
         P[i-1]=P[i]+'*'+summ
      P[0]=P[1]+'*'+str(2**m[1])+'/'+str(math.factorial(m[1]+1)*math.factorial(m[1]))
   else:
      P[NumSpec-1]=1
      for i in range(NumSpec-1,1,-1):
         summ=0
         for j in range(m[i]+1):
            summ+=math.exp(-l[(i,j)])**timeDict[i]/den(RH,i,j)
         P[i-1]=P[i]*summ
      P[0]=P[1]*2**m[1]/(math.factorial(m[1]+1)*math.factorial(m[1]))
   return(P[0])

for st in dataST:
	S_calc(st)
	if InputFormat=='1':
		specie=read(st,S_species)
	if IntervalMode=='1':
		Time_Yule(float(SpecRatePar))
	elif InputFormat=='2':
		specie=newick_analysis(st)
	topocalc()
	if IntervalMode!='3':
		strTimeDict='('+str(timeDict[2])
		for i in range(3,NumSpec):
			strTimeDict+=','+str(timeDict[i])
			strTimeDict+=')'
	if OutputFormat=='1':
		print('\nST = '+st+'\nTimes = '+strTimeDict+'\n')
	if OutputFormat=='3':
		if IntervalMode=='3':
			filename=st+'_prob.nb'
		elif IntervalMode in ['1','2']:
			filename=st+'_'+strTimeDict+'_prob.txt'
		else:
			filename=st+'_prob.txt'
		fout=open(filename, 'w')
		if IntervalMode!='3':
			fout.write('Times= '+strTimeDict+'\n\n')
	for gt in dataGT:
		MRH=findaMRH(st,gt)
		rhlist=rhcalc(MRH)
		if OutputFormat=='2':
			if IntervalMode=='3':
				filename=st+'_'+gt+'_prob.nb'
			elif IntervalMode in ['1','2']:
				filename=st+'_'+gt+'_'+strTimeDict+'_prob.txt'
			else:
				filename=st+'_'+gt+'_prob.txt'
			fout=open(filename, 'w')
			if IntervalMode!='3':
				fout.write('Times = '+strTimeDict+'\n\n')
		if IntervalMode=='3':
			pt=''
		else:
			pt=0
		if OutputProbab=='2':
			probs=''
		for rh in rhlist:
			p=probab(st,gt,rh)
			if IntervalMode=='3':
				pt+='+'+p
				if OutputProbab=='2':
					probs+='Print["P'+str(rh)+' = ", Simplify['+str(p)+']]\n'
			else:
				pt+=p
				if OutputProbab=='2':
					probs+='P'+str(rh)+' = '+str(p)+'\n'
		if IntervalMode=='3':
			pt=pt[1:]
			fout.write('Print["P('+gt+'|'+st+') = ", Simplify['+str(pt)+']]\n')
		elif OutputFormat=='1':
			print('P('+gt+')= '+str(pt)+'\n')
		else:
			fout.write('P('+gt+') = '+str(pt))
		if OutputProbab=='2':
			if OutputFormat=='1':
				print(probs)
			else:
				fout.write(probs)
		if OutputFormat=='2':
			fout.close()
	if OutputFormat=='3':
		fout.close()

r=input('to close the program type "0" and submit\n')
if r=='0':
	pass
