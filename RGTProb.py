#Project developed by Guido Narduzzi and Pasquale Miglionico, ordinary students of Biology at the Scuola Normale Superiore of Pisa
#under the direction of Prof. Filippo Di Santo, University of Pisa, Department of Mathematics


print('''Welcome to RGTProb, the gene tree probability calculator. This tool will allow you to calculate Ranked Gene Tree
         probabilities conditioning on a species tree in multiple formats. You can modify the program as you like, but please always mention the
         original project. If you develop something interesting please let us know, so that we can add it to our software and distribute
         your content to other users.''')

#Main code

from decimal import *
import math

Top=''
gene=''
topologia=''
nodi=''
dt=dict()
l=dict()
IntervalMode='000'
S_specie=''

InputFormat=input('Select the desired input format: type 1 for standard format or 2 for newick format\n') #Inp1
InputType=input('Select how to submit inputs: type 1 to submit files or 2 to type trees manually\n') #Inp2
NumSpec=int(input('Select the number of species\n'))

#Functions for time models

def Time_Yule(par):
   global dt, N
   for i in range(NumSpec-2):
      dt[i+2]=1/((i+2)*par)
   return dt

###

if InputFormat=='1': #standard format
   IntervalMode=input("Select the desired time interval mode: type 1 for Yule's model, 2 to insert time intervals manually or 3 for a polynomial output (readable by Wolfram Mathematica)\n")#Inp4
   if IntervalMode=='1':#creates a standard set of times according to the function Time_Yule
      SpecRatePar=input('Insert the parameter for the speciation rate\n')
      Time_Yule(float(SpecRatePar))
   elif IntervalMode=='2': # insert intervals
      times_raw=input('Write the desired time intervals from T(2) to T(N-1) only separated by a comma\n')
      times=times_raw.split(',')
      for j in range (N-2):
         n=j+2
         dt[n]=float(times[j])
   elif IntervalMode=='3': # polynomials
      pass
   else:
      print('ERROR: invalid input')


if InputType=='1': #process files
   fileST=input('Type the name of the file containing STs as "filename.txt"\n')
   openST=open(fileST)
   dataST=[]
   fileGT=input('Type the name of the file containing GTs as "filename.txt"\n')
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
   speciestree=input('Type the ST:\n')
   dataST=[speciestree]
   genetree=input('Type the GT:\n')
   dataGT=[genetree]

#now we have our data in a way that allows us to use only one script for both cases. Now the desired output format must be chosen

OutputFormat=input('''Choose the desired output format: type 1 for desktop print (this might not work if it has too much to print),
             2 for 1 file for each ST-GT couple or 3 for 1 file for each ST.\n''') #notare che potrebbe funzionare bene in coll diretto a mathematica solo il 2
#inp6
OutputProbabs=input('Type 1 if you want only total GT|ST probabilities, 2 if you want also GT&RH|ST probabilities\n')#input7

#Now all necessary information has been collected, we can write functions and the rest of the program

def T_calc(Tree): #calculates a vector indicating the number of nodes which separe the leafs from the root
	n=0
	T=[]
	for i in range(len(Tree)):
		if Tree[i]=='(':
			n=n+1
		if Tree[i]==')':
			n=n-1
		T.append(n)
	return T

def newick_analysis(ST): #newick processing, it extrapolates the dictionary of times and the matrix M
	global dt
	K=[]
	Val=[]
	tempi=set()
	M=[]
	T=T_calc(ST)
	for i in range(len(ST)): #the vector Val contains the branch lenghts
		if ST[i] in [':']:
			j=i+1
			stringa=''
			while ST[j] in ['1','2','3','4','5','6','7','8','9','0','.']:
				stringa=stringa+ST[j]
				j=j+1
			Val.append(Decimal(stringa))
		else:
			Val.append(0)
	for i in range(N-1):
		K.append([])
		for j in range(N-1-i):
			n=min(T[S_specie[i]:S_specie[j+i+1]])
			h=S_specie[i]+1
			t=Val[h]
			for k in range(T[S_specie[i]]-n):
				while T[h]!=T[S_specie[i]]-k-1:
					h+=1
				h+=1
				t=t+Val[h]
			K[i].append(t)
			tempi.add(t)
	tem=list(tempi)
	tem.sort(reverse=True)
	tem.append(0)
	dt=dict()
	for i in range(N-1):
		M.append([])
		dt[i+2]=float(tem[i]-tem[i+1])
		for j in range(N-1-i):
			h=0
			while tem[h]!=K[i][j]:
				h+=1
			M[i].append(h+1)
	return M

def trovaMRH(ST,GT):#its inputs are a RST and a RGT. the function returns the MRH of th RGT in the RST as a list
   global N,Top,gene,topologia,nodi,S_specie,Inp1
   #this function only uses the function lettura
   if InputFormat=='1':	   
      specie=lettura(ST,S_specie)
   elif Inp1=='2':
      specie=newick_analysis(ST)
   S_gene=[] #list containing the position of letters in GT in a way s.t. GT[S_gene[i]]=ST[S_specie[i]]
   for i in range(N):
      n=0
      while GT[n]!=ST[S_specie[i]]:
         n+=1
      S_gene.append(n)
   gene=lettura(GT,S_gene)
   #now we calculate the MRH. L, Top, topologia, nodi are created at the same time
   L=[] #L[i] will be the time of the (i+1)th node of the gene tree
   topologia=[[]] #topologia[i][j] is the branch containing the (i+1) phyletic line during the (j+1)^th time interval
   Top=[]
   nodi=[]
   for i in range(N-1):
      L.append([])
      topologia.append([])
      topologia[0].append(1)
   for i in range(N-1):
      for j in range(N-i-1):
         L[gene[i][j]-1].append(specie[i][j])
      Top.append([])
      n=0
      for j in range(N-1):
         if specie[i][0]>j:
            topologia[i+1].append(topologia[i][j])
         else:
            topologia[i+1].append(topologia[i][j]+1)
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
      nodi.append(k)
   MRH=[L[N-2][0]]
   for i in range(N-2):
      MRH.append(min(MRH[i],min(L[N-3-i])))
   MRH.sort()
   return MRH

def lettura(stringa,S):
   T=T_calc(stringa)
   M=[]
   for i in range(N-1):
      M.append([])
      for j in range(N-1-i):
         if S[i]<S[i+j+1]:
            n=min(T[S[i]:S[j+i+1]])
         else:
            n=min(T[S[j+i+1]:S[i]])
         k=max(S[i],S[j+i+1])
         while T[k]>=n:
            k+=1
         if k+2==len(stringa):
            M[i].append(1)
         elif stringa[k+2] in ['1','2','3','4','5','6','7','8','9','0']:
            M[i].append(int(stringa[k+1:k+3]))
         else:
            M[i].append(int(stringa[k+1]))
   return M

def calcolarh(tup): #input: MRH. Output: a list with all possible rhs. Uses aggiorna_rh
   D=[(1,)]
   j=1
   while j<len(tup):
      D=aggiorna_rh(ins=D,Z=tup[j])
      j+=1
   return D

def aggiorna_rh(ins,Z=1000): #Input: a list of RHs. Output: a list of RHs which are all those obtained by elongating of 1 element the
   #input with valid values
   d=[]
   for chiave in ins:
      N=len(chiave)
      break
   T=min(N+1,Z)
   for chiave in ins:
      n=0
      M=max(chiave)
      while n<M:
         n+=1
      while n<=T:
         newkey=chiave+(n,)
         d.append(newkey)
         n+=1
   return d

def probab(ST,GT,RH): #Input: ST, GT and a RH. Output: P(G and RH|ST). Uses trovam, vall
   global l, dt, N, Inp1, Inp4
   calcoLambda(RH)
   P=dict()
   m=trovam(RH)
   P[N-1]=1#it is necessarily so
   nums=list(range(N-1))
   if InputFormat=='1' and IntervalMode=='3':
      tau=dict()
      for t in range(N-2):
         tau[t+2]='T['+str(t+2)+']'
      polystr=''
      P[N-1]='1'
      for i in nums[::-1]:
         if i==0:
            break
         summ='('
         for j in range(m[i+1]+1):
            summ+='('+str(num(RH,i,j))+'*'+tau[i+1]+'])*('+str(den(RH,i,j))+')+'
         summ=summ[0:len(summ)-1]
         summ+=')'
         P[i]=P[i+1]+'*'+summ
         if i==1: #it stops at the P from s(N-1) to s(1), because P in tau_1 is taken care of in the next part of the function
            break
      P[0]=P[1]+'*'+'(2)^('+str(m[1])+')'+'/('+str(m[1]+1)+'!*'+str(m[1])+'!)'
   else:
      for i in nums[::-1]:
         if i==0:
            break
         summ=0
         for j in range(m[i+1]+1):
            summ+=(num(RH,i,j))**(dt[i+1])*den(RH,i,j)
         P[i]=P[i+1]*summ
         if i==1: #it stops at the P from s(N-1) to s(1), because P in tau_1 is taken care of in the next part of the function
            break
      P[0]=P[1]*2**(m[1])/(math.factorial(m[1]+1)*math.factorial(m[1]))
   return(P[0])

def den(rh,i,j): #uses trovam, calcoLambda
   global l, Inp4
   m=trovam(rh)
   d1=1
   if Inp4=='3':
      d1='('
   for k in range(m[i+1]+1):
      if k!=j:
         if Inp4=='3':
            d1+=str(l[(i+1,k)]-l[(i+1,j)])+'*'
         else:
            d1=d1*(l[(i+1,k)]-l[(i+1,j)])
   if Inp4=='3':
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
   global l, Inp4
   if Inp4=='3':
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

def valori(RH): #Output: a dictionary of k[ijz], it is useful to define m[N]=0 and k(N,0,z)=1
   global N,Top,gene,topologia,nodi
   Pos=[] #a matrix s.t. Pos[i][j] is the branch in which rhe (j+1)^th coalescence happens in tau_(i+1)
   for i in range(N-1):
      Pos.append([])
      Pos[RH[i]-1].append(topologia[nodi[i]][RH[i]-1])
   m=[]
   k=dict()
   for i in range(N):
      m.append(RH.count(i+1))
      k[(N,0,i+1)]=1
   for i in range(N-1):
      for z in range (N-1-i):
         k[(N-1-i,m[N-2-i],z+1)]=0
      for j in range(m[N-2-i]+1):
         if j==0:
            for z in range(N-i):
               k[(N-1-i,m[N-2-i],Top[N-2-i][z])]+=k[(N-i,0,z+1)]
         else:
            for z in range(N-1-i):
               if z+1==Pos[N-2-i][m[N-2-i]-j]:
                  k[(N-1-i,m[N-2-i]-j,z+1)]=k[(N-1-i,m[N-2-i]-j+1,z+1)]-1
               else:
                  k[(N-1-i,m[N-2-i]-j,z+1)]=k[(N-1-i,m[N-2-i]-j+1,z+1)]
   for i in range(N):
      del(k[(N,0,i+1)])
   return k

def S_calc(st): #finds the positions of the species names in the input string
	global S_specie
	S_specie=[]
	for i in range(len(st)):
		if st[i] not in ['1','2','3','4','5','6','7','8','9','0','.',':',',','(',')']:
			S_specie.append(i)
	return

if Inp6=='1': #at the moment it works with standard format, symbolic must be implemented, newick corrected
   for st in dataST:
      S_calc(st)
      if Inp1=='2':
         newick_analysis(st) #così abbiamo dt e M
      strdt='('
      if Inp4!='3':
         strdt+=str(dt[2])
         for i in range(N-3):
            strdt+=','+str(dt[i+3])
         strdt+=')'
      if strdt!='(':
         print('Times= '+strdt)
      for gt in dataGT:
         MRH=trovaMRH(st,gt)
         rhlist=calcolarh(MRH)
         if Inp4=='3':
            pt=''
         else:
            pt=0
         probs=''
         print('ST= '+st+' GT= '+gt)
         for rh in rhlist:
            p=probab(st,gt,rh)
            if Inp4=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+=str(rh)+': '+str(p)+'; '
         if Inp4=='3':
            pt=pt[1:]
         print('P(GT|ST)= '+str(pt))
         if Inp7=='2':
            print(probs)
         print()
      print()

if Inp6=='2': #for each gt-st couple
   for st in dataST:
      S_calc(st)
      if Inp1=='2':
         newick_analysis(st) #così abbiamo dt
      strdt=''
      if Inp4!='3':
         strdt+='('
         strdt+=str(dt[2])
         for i in range(N-3):
            strdt+=','+str(dt[i+3])
         strdt+=')'
      for gt in dataGT:
         MRH=trovaMRH(st,gt)
         rhlist=calcolarh(MRH)
         if Inp1=='1':
            filename=st+'_'+gt+'_'+strdt+'_prob.txt'
         else:
            stname=''
            for letter in st:
               if letter not in '0123456789,.:':
                  stname+=letter
                  stname+='newick'
            filename=stname+'_'+gt+'_'+strdt+'_prob.txt'
         fout=open(filename, 'w')
         fout.write('GT='+gt+'\nST='+st+'\n')
         if Inp4!='3':
            fout.write('Times='+strdt+'\n')
            pt=0
         if Inp4=='3':
            pt=''
         probs=''
         for rh in rhlist:
            p=probab(st,gt,rh)
            if Inp4=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+=str(rh)+': '+str(p)+'\n'
         if Inp4=='3':
            pt=pt[1:]
         fout.write('P(GT|ST)= '+str(pt)+'\n')
         if Inp7=="2":
            fout.write(probs)
         fout.close()
         
if Inp6=='3':
   for st in dataST:
      S_calc(st)
      if Inp1=='2':
         newick_analysis(st)
      strdt=''
      if Inp4!='3':
         strdt+='('
         strdt+=str(dt[2])
         for i in range(N-3):
            strdt+=','+str(dt[i+3])
         strdt+=')'
      if Inp1=='1':
         filename=st+'_'+strdt+'_prob.txt'
      else:
         stname=''
         for letter in st:
            if letter not in '0123456789,.:':
               stname+=letter
               stname+='newick'
         filename=stname+'_'+strdt+'_prob.txt'
      fout=open(filename, 'w')
      fout.write('ST= '+st+'\n')
      if Inp4!='3':
         fout.write('Times= '+strdt+'\n\n')
      for gt in dataGT:
         MRH=trovaMRH(st,gt)
         rhlist=calcolarh(MRH)
         fout.write('GT='+gt+'\n')
         if Inp4=='3':
            pt=''
         else:
            pt=0
            fout.write('Times='+strdt+'\n')
         probs=''
         for rh in rhlist:
            p=probab(st,gt,rh)
            if Inp4=='3':
               pt+='+'+p
            else:
               pt+=p
            probs+=str(rh)+': '+str(p)+'\n'
         if Inp4=='3':
            pt=pt[1:]
         fout.write('P(GT|ST)= '+str(pt)+'\n')
         if Inp7=='2':
            fout.write(probs)
         fout.write('\n')
      fout.close()


r=input('to close the program type "0" and submit\n')
if r=='0':
   pass
