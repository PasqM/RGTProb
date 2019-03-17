#Project developed by Guido Narduzzi and Pasquale Miglionico, undergraduate fellows of the Scuola Normale Superiore of Pisa and students of 
#the Department of Biology of the University of Pisa, under the direction of Dr. Filippo Di Santo, University of Pisa, Department of 
#Mathematics

'''Commenti su cose da fare/migliorare:
1- La funzione newick analysis va resa più leggibile, le variabili interne devono essere più leggibili. Boh, posso tradurre i nomi italiani, ma certe a certe variabili non riesco a dare un nome self-explaining
2- Ma col newick non possiamo fare i polinomi??? Boh, perchè no? La domanda è perchè qualcuno vorrebe usare newick per farsi calcolare dei polinomi
3- Ma newick che formato indica alla fine??? Newick è quello con le lunghezze dei rami esplicite, se serve cambiamo i nomi
La revisione arriva fino a update_rh (circa metà programma
'''

print('''Welcome to RGTProb, the gene tree probability calculator. This tool will allow you to calculate Ranked Gene Tree
         probabilities conditioning on a species tree in multiple formats. You can modify the program as you like, but please always mention the
         original project. If you develop something interesting please let us know, so that we can add it to our software and distribute
         your content to other users.''')

#Main code

from decimal import *
import math

#Initialization of useful variables

Top=''
topology=''
nodes=''
timeDict=dict()#A dictionary containing the lengths of the branches
l=dict()
IntervalMode='000' #Specifies the mode in which the interval lengths should be determined
S_species=''#Vector of the positions of the species names in the input string

InputFormat=input('Select the desired input format: type 1 for standard format or 2 for newick format\n')
InputType=input('Select how to submit inputs: type 1 to submit files or 2 to type trees manually\n')
NumSpec=int(input('Select the number of species\n'))

#Checking for valid inputs

if InputFormat not in ['1','2']:
   print('Error, invalid input format request!\nPlease restart the program.')

if InputType not in ['1','2']:
   print('Error, invalid input submission mode request!\nPlease restart the program.')

### Preliminary data acquisition, and preference settings

if InputFormat=='1': #standard format: needs to create a dictionary of time lengths 
   IntervalMode=input("Select the desired time interval mode: type 1 for Yule's model, 2 to insert time intervals manually or 3 for a polynomial output (readable by Wolfram Mathematica)\n")
   if IntervalMode=='1':#creates a standard set of times according to the function Time_Yule
      SpecRatePar=input('Insert the parameter for the speciation rate\n')
      Time_Yule(float(SpecRatePar))
   elif IntervalMode=='2': # insert intervals manually
      times_raw=input('Write the desired time intervals from T(2) to T(N-1) only separated by a comma\n')
      times=times_raw.split(',')
      for j in range (NumSpec-2):
         n=j+2
         timeDict[n]=float(times[j])
   elif IntervalMode=='3': # polynomials
      pass
   else:
      print('ERROR: invalid input. Please restart the program.')


if InputType=='1': #processing files to extract the trees on which to perform the calculations
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

if InputType=='2': #Allows for the manual insertion of GTs and STs
   print('Please write trees without blank spaces and characters that are not letters, numbers parenthesis: an example of correct typing is ((AB)2C)1')
   speciestree=input('Type the SPECIES tree:\n')
   dataST=[speciestree]
   genetree=input('Type the GENE tree:\n')
   dataGT=[genetree]

def Time_Yule(par):#function that sets the time intervals according tu Yule's model
   global timeDict, NumSpec
   for i in range(NumSpec-2):
      timeDict[i+2]=1/((i+2)*par)
   return timeDict

#now we have our data in a way that allows us to use the same functions for both cases. Now the desired output format must be chosen

OutputFormat=input('''Choose the desired output format: type 1 for desktop print (not optimal for big trees or a collection of trees),
             2 for 1 file for each ST-GT couple or 3 for 1 file for each ST.\n''')
OutputProbab=input('Type 1 if you only want total GT|ST probabilities, 2 if you also want GT&RH|ST probabilities\n')

#The following functions calculate the matrix M[i][j] containing the rank of the "last common ancestor" (from now on LCA) of i and j

def S_calc(st): #calculates the vector S_species of the positions of the species names in the input string st
	global S_species
	S_species=[]
	for i in range(len(st)):
		if st[i] not in ['1','2','3','4','5','6','7','8','9','0','.',':',',','(',')']:
			S_species.append(i)
	return

def T_calc(Tree): #calculates the vector LRd (Leaves-Root distance) vector indicating the number of nodes separating the leaves from the root
	n=0
	LdR=[]
	for i in range(len(Tree)):
		if Tree[i]=='(':
			n=n+1
		if Tree[i]==')':
			n=n-1
		LdR.append(n)
	return LdR

def newick_analysis(ST): #Processing of the species tree when given in newick format, it also extrapolates the dictionary of times
	global timeDict
	K=[]
	BranchL=[]
	times=set()
	M=[]
	LdR=T_calc(ST)
	for i in range(len(ST)): #the vector BranchL contains the branch lenghts
		if ST[i] in [':']:
			j=i+1
			string=''
			while ST[j] in ['1','2','3','4','5','6','7','8','9','0','.']:
				string=string+ST[j]
				j=j+1
			BranchL.append(Decimal(string))
		else:
			BranchL.append(0)
	for i in range(NumSpec-1): #the matrix K[i][j] contains the sum from T(N) to T(s+1) where s is the rank of the LCA of i and j
		K.append([])
		for j in range(NumSpec-1-i):
			n=min(LdR[S_species[i]:S_species[j+i+1]])
			h=S_species[i]+1
			t=BranchL[h]
			for k in range(LdR[S_species[i]]-n):
				while LdR[h]!=LdR[S_species[i]]-k-1:
					h+=1
				h+=1
				t=t+BranchL[h]
			K[i].append(t)
			times.add(t)
	tim=list(times)
	tim.sort(reverse=True)
	tim.append(0)
	timeDict=dict()
	for i in range(NumSpec-1):
		M.append([])
		timeDict[i+2]=float(tim[i]-tim[i+1])
		for j in range(NumSpec-1-i):
			h=0
			while tim[h]!=K[i][j]:
				h+=1
			M[i].append(h+1)
	return M
	
def read(string,S): #Processing of a species tree when provided in standard format
   LdR=T_calc(string)
   M=[]
   for i in range(NumSpec-1):
      M.append([])
      for j in range(NumSpec-1-i):
         if S[i]<S[i+j+1]:
            n=min(LdR[S[i]:S[j+i+1]])
         else: #this is useful to process GTs in which the order of the species is different from the ST
            n=min(LdR[S[j+i+1]:S[i]])
         k=max(S[i],S[j+i+1])
         while LdR[k]>=n:
            k+=1
         if k+2==len(string):
            M[i].append(1)
         elif string[k+2] in ['1','2','3','4','5','6','7','8','9','0']:
            M[i].append(int(string[k+1:k+3]))
         else:
            M[i].append(int(string[k+1]))
   return M

#Now we can calculate the Maximal Ranked History (MRH) of a GT in a ST and the vector nodes and the matrices topology and Top.

def findaMRH(ST,GT):
   global Top,gene,topology,nodes
   if InputFormat=='1':
      specie=read(ST,S_species)
   elif InputFormat=='2':
      specie=newick_analysis(ST)
   S_gene=[] #list containing the position of letters in GT in a way s.t. GT[S_gene[i]]=ST[S_species[i]]
   for i in range(NumSpec):
      n=0
      while GT[n]!=ST[S_species[i]]:
         n+=1
      S_gene.append(n)
   gene=read(GT,S_gene)
   L=[] #L[i] will contain the ranks of the LCAs (according to the ST) of the couples of species whose LCA in the GT has rank i
   topology=[[]] #topology[i][j] is the branch containing the (i+1) phyletic line during the (j+1)^th time interval
   Top=[]
   nodes=[]#nodes[i] indicates one of the species descending from the node with rank i+1
   for i in range(NumSpec-1):
      L.append([])
      topology.append([])
      topology[0].append(1)
      Top.append([])
   for i in range(NumSpec-1):
      for j in range(NumSpec-i-1):
         L[gene[i][j]-1].append(specie[i][j])
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
      k=0
      while i+1 not in gene[k]:
         k+=1
      nodes.append(k)
   MRH=[L[NumSpec-2][0]]
   for i in range(NumSpec-2):
      MRH.append(min(MRH[i],min(L[NumSpec-3-i])))
   MRH.sort()
   return MRH

def calcolarh(tup): #input: MRH. Output: a list with all possible rhs. Uses update_rh
   D=[(1,)]
   j=1
   while j<len(tup):
      D=update_rh(D,tup[j])
      j+=1
   return D

def update_rh(ins,Z): #Input: a list of RHs. Output: a list of RHs which are all those obtained by elongating of 1 element the
   #input with valid values
   d=[]#a list of all ranked histories
   for key in ins:
      n=max(key)
      while n<=Z:
         newkey=key+(n,)
         d.append(newkey)
         n+=1
   return d

def probab(ST,GT,RH): #Input: ST, GT and a RH. Output: P(G and RH|ST). Uses trovam, vall
   global l, timeDict, NumSpec, InputFormat, IntervalMode
   calcoLambda(RH)
   P=dict()
   m=trovam(RH)
   P[NumSpec-1]=1#it is necessarily so
   nums=list(range(NumSpec-1))
   if InputFormat=='1' and IntervalMode=='3':
      tau=dict()
      for t in range(NumSpec-2):
         tau[t+2]='T['+str(t+2)+']'
      polystr=''
      P[NumSpec-1]='1'
      for i in nums[::-1]:
         if i==0:
            break
         summ='('
         for j in range(m[i+1]+1):
            summ+='('+str(num(RH,i,j))+'*'+tau[i+1]+'])*('+str(den(RH,i,j))+')+'
         summ=summ[0:len(summ)-1]
         summ+=')'
         P[i]=P[i+1]+'*'+summ
         if i==1: #it stops at the P from s(NumSpec-1) to s(1), because P in tau_1 is taken care of in the next part of the function
            break
      P[0]=P[1]+'*'+'(2)^('+str(m[1])+')'+'/('+str(m[1]+1)+'!*'+str(m[1])+'!)'
   else:
      for i in nums[::-1]:
         if i==0:
            break
         summ=0
         for j in range(m[i+1]+1):
            summ+=(num(RH,i,j))**(timeDict[i+1])*den(RH,i,j)
         P[i]=P[i+1]*summ
         if i==1: #it stops at the P from s(NumSpec-1) to s(1), because P in tau_1 is taken care of in the next part of the function
            break
      P[0]=P[1]*2**(m[1])/(math.factorial(m[1]+1)*math.factorial(m[1]))
   return(P[0])

def den(rh,i,j): #uses trovam, calcoLambda
   global l, IntervalMode
   m=trovam(rh)
   d1=1
   if IntervalMode=='3':
      d1='('
   for k in range(m[i+1]+1):
      if k!=j:
         if IntervalMode=='3':
            d1+=str(l[(i+1,k)]-l[(i+1,j)])+'*'
         else:
            d1=d1*(l[(i+1,k)]-l[(i+1,j)])
   if IntervalMode=='3':
      if d1=='(':
         den='1'
      else:
         d1=d1[:len(d1)-1]
         d1+=')'
         den='(1/'+d1+')'
   else:
      den=1/d1
   return den

def num(rh,i,j):
   global l, IntervalMode
   if IntervalMode=='3':
      num='Exp['+str(-l[(i+1,j)])
   else:
      num=math.exp(-l[(i+1,j)])
   return num

def trovam(RH): #Output: a dictionary of m, where m[i] is the number of nodes of the gene tree in the i-th interval of the species tree
   m=dict()
   for n in range(len(RH)+1):
      m[n]=RH.count(n)
   return m

def calcoLambda(RH): #Input: a dictionary of k[ijz]. Output: a dictionary of lambda(i,j). Uses valori, binomiale
   global l
   l=dict()
   k=valori(RH)
   for chiave in k:
      ij=chiave[0:2]
      l[ij]=l.get(ij,0)+binomiale(k[chiave],2)

def binomiale(n,k): #calculates the binomial n over k
   if n<k:
      return 0
   return math.factorial(n)//(math.factorial(k)*math.factorial(n-k))

def valori(RH): #Output: a dictionary of k[ijz], it is useful to define m[NumSpec]=0 and k(NumSpec,0,z)=1
   global NumSpec,Top,topology,nodes
   Pos=[] #a matrix s.t. Pos[i][j] is the branch containing the (j+1)^th coalescence in tau_(i+1)
   for i in range(NumSpec-1):
      Pos.append([])
      Pos[RH[i]-1].append(topology[nodes[i]][RH[i]-1])
   m=[]
   k=dict()
   for i in range(NumSpec):
      m.append(RH.count(i+1))
      k[(NumSpec,0,i+1)]=1
   for i in range(NumSpec-1):
      for z in range (NumSpec-1-i):
         k[(NumSpec-1-i,m[NumSpec-2-i],z+1)]=0
      for j in range(m[NumSpec-2-i]+1):
         if j==0:
            for z in range(NumSpec-i):
               k[(NumSpec-1-i,m[NumSpec-2-i],Top[NumSpec-2-i][z])]+=k[(NumSpec-i,0,z+1)]
         else:
            for z in range(NumSpec-1-i):
               if z+1==Pos[NumSpec-2-i][m[NumSpec-2-i]-j]:
                  k[(NumSpec-1-i,m[NumSpec-2-i]-j,z+1)]=k[(NumSpec-1-i,m[NumSpec-2-i]-j+1,z+1)]-1
               else:
                  k[(NumSpec-1-i,m[NumSpec-2-i]-j,z+1)]=k[(NumSpec-1-i,m[NumSpec-2-i]-j+1,z+1)]
   for i in range(NumSpec):
      del(k[(NumSpec,0,i+1)])
   return k

if OutputFormat=='1': #at the moment it works with standard format, symbolic must be implemented, newick corrected
   for st in dataST:
      S_calc(st)
      if InputFormat=='2':
         newick_analysis(st) #così abbiamo timeDict e M
      strTimeDict='('
      if IntervalMode!='3':
         strTimeDict+=str(timeDict[2])
         for i in range(NumSpec-3):
            strTimeDict+=','+str(timeDict[i+3])
         strTimeDict+=')'
      if strTimeDict!='(':
         print('Times= '+strTimeDict)
      for gt in dataGT:
         MRH=findaMRH(st,gt)
         rhlist=calcolarh(MRH)
         if IntervalMode=='3':
            pt=''
         else:
            pt=0
         probs=''
         print('ST= '+st+' GT= '+gt)
         for rh in rhlist:
            p=probab(st,gt,rh)
            if IntervalMode=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+=str(rh)+': '+str(p)+'; '
         if IntervalMode=='3':
            pt=pt[1:]
         print('P(GT|ST)= '+str(pt))
         if OutputProbab=='2':
            print(probs)
         print()
      print()

if OutputFormat=='2': #for each gt-st couple
   for st in dataST:
      S_calc(st)
      if InputFormat=='2':
         newick_analysis(st) #così abbiamo timeDict
      strTimeDict=''
      if IntervalMode!='3':
         strTimeDict+='('
         strTimeDict+=str(timeDict[2])
         for i in range(NumSpec-3):
            strTimeDict+=','+str(timeDict[i+3])
         strTimeDict+=')'
      for gt in dataGT:
         MRH=findaMRH(st,gt)
         rhlist=calcolarh(MRH)
         if InputFormat=='1':
            filename=st+'_'+gt+'_'+strTimeDict+'_prob.nb'
         else:
            stname=''
            for letter in st:
               if letter not in '0123456789,.:':
                  stname+=letter
                  stname+='newick'
            filename=stname+'_'+gt+'_'+strTimeDict+'_prob.nb'
         fout=open(filename, 'w')
         if IntervalMode!='3':
            fout.write('Times='+strTimeDict+'\n')
            pt=0
         if IntervalMode=='3':
            pt=''
         probs=''
         for rh in rhlist:
            p=probab(st,gt,rh)
            if IntervalMode=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+='Print["P('+str(rh)+')=", Expand['+str(p)+']]\n'
         if IntervalMode=='3':
            pt=pt[1:]
         fout.write('Print["P('+gt+'|'+st+')= ", Expand['+str(pt)+']]\n')
         if OutputProbab=="2":
            fout.write(probs)
         fout.close()
         
if OutputFormat=='3':
   for st in dataST:
      S_calc(st)
      if InputFormat=='2':
         newick_analysis(st)
      strTimeDict=''
      if IntervalMode!='3':
         strTimeDict+='('
         strTimeDict+=str(timeDict[2])
         for i in range(NumSpec-3):
            strTimeDict+=','+str(timeDict[i+3])
         strTimeDict+=')'
      if InputFormat=='1':
         filename=st+'_'+strTimeDict+'_prob.nb'
      else:
         stname=''
         for letter in st:
            if letter not in '0123456789,.:':
               stname+=letter
               stname+='newick'
         filename=stname+'_'+strTimeDict+'_prob.nb'
      fout=open(filename, 'w')
      if IntervalMode!='3':
         fout.write('Times= '+strTimeDict+'\n\n')
      for gt in dataGT:
         MRH=findaMRH(st,gt)
         rhlist=calcolarh(MRH)
         if IntervalMode=='3':
            pt=''
         else:
            pt=0
            fout.write('Times='+strTimeDict+'\n')
         probs=''
         for rh in rhlist:
            p=probab(st,gt,rh)
            if IntervalMode=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+='Print["P('+str(rh)+')=", Expand['+str(p)+']]\n'
         if IntervalMode=='3':
            pt=pt[1:]
         fout.write('Print["P('+gt+'|'+st+')= ", Expand['+str(pt)+']]\n')
         if OutputProbab=='2':
            fout.write(probs)
         fout.write('\n')
      fout.close()


r=input('to close the program type "0" and submit\n')
if r=='0':
	pass
