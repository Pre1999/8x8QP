import numpy as np

##################################################################

#Function to Solve for X and Y
def X_Y(Nlist,Plist):
	cmatrix=np.zeros((gates,gates),dtype=float)
	bx=np.zeros((gates,1),dtype=float)
	by=np.zeros((gates,1),dtype=float)
	X=np.zeros((gates,1),dtype=float)
	Y=np.zeros((gates,1),dtype=float)
	
	#Forming Connectivity Matrix
	for i in range(1,len(Nlist)):
		for j in range(1,len(Nlist[i])):
			for k in range(j+1,len(Nlist[i])):
				if (len(Nlist[i])>2):
					cmatrix[Nlist[i][j]-1][Nlist[i][k]-1]+=1/(len(Nlist[i])+Nlist[i][0]-2)
					cmatrix[Nlist[i][k]-1][Nlist[i][j]-1]+=1/(len(Nlist[i])+Nlist[i][0]-2)
	
	#Forming A,Bx and By Matrix
	for i in range(len(cmatrix)):
		sumc=0
		for j in range(len(cmatrix)):
			sumc+=cmatrix[i][j]
			cmatrix[i][j]=-cmatrix[i][j]
		cmatrix[i][i]=sumc
	
	for i in range(1,len(Plist)):
		for j in range(1,len(Nlist[Plist[i][1]])):
			cmatrix[Nlist[Plist[i][1]][j]-1][Nlist[Plist[i][1]][j]-1]+=1/(len(Nlist[Plist[i][1]])+Nlist[Plist[i][1]][0]-2)
			bx[Nlist[Plist[i][1]][j]-1]+=Plist[i][2]*(1/(len(Nlist[Plist[i][1]])+Nlist[Plist[i][1]][0]-2))
			by[Nlist[Plist[i][1]][j]-1]+=Plist[i][3]*(1/(len(Nlist[Plist[i][1]])+Nlist[Plist[i][1]][0]-2))

	#Solving For X and Y
	Xc=np.linalg.lstsq(cmatrix,bx,rcond=None)
	Yc=np.linalg.lstsq(cmatrix,by,rcond=None)
	
	for i in range(gates):
		X[i][0]=round(Xc[0][i][0],8)
		Y[i][0]=round(Yc[0][i][0],8)
	return(X,Y)

##################################################################

#Containment Function
#	!!	(xmin & xmax , ymin & ymax) -> Values in which contained gates are present
def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax,shift):
	for i in range(len(contgates)):
		flagc=1
		for j in range(3,len(Glist[contgates[i]-1])):
			for m in range(1,len(Nlist[Glist[contgates[i]-1][j]])):
				for k in range(len(uncontgates)):
					if(Nlist[Glist[contgates[i]-1][j]][0]==1 and flagc==1):
						for n in range(1,len(Plist)):
							if(Plist[n][1]==Glist[contgates[i]-1][j]):
								#print("Gate",contgates[i],"is connected to Pad",n)
								if(Plist[n][3]>=ymax and Plist[n][2]<xmin):				#  __|
									temp=[0,Plist[n][1],xmin,ymax]
								elif(Plist[n][3]>=ymax and Plist[n][2]>=xmin and Plist[n][2]<xmax):	# |__|
									temp=[0,Plist[n][1],Plist[n][2],ymax]
								elif(Plist[n][3]>=ymax and Plist[n][2]>=xmax):			# |__
									temp=[0,Plist[n][1],xmax,ymax]
								elif(Plist[n][3]<ymax and Plist[n][3]>=ymin and Plist[n][2]>=xmax):	# | =
									temp=[0,Plist[n][1],xmax,Plist[n][3]]
								elif(Plist[n][3]<ymin and Plist[n][2]>=xmax):				# |""
									temp=[0,Plist[n][1],xmax,ymin]
								elif(Plist[n][3]<=ymin and Plist[n][2]>=xmin and Plist[n][2]<xmax):	# |""|
									temp=[0,Plist[n][1],Plist[n][2],ymin]
								elif(Plist[n][3]<ymin and Plist[n][2]<xmin):				#  ""|
									temp=[0,Plist[n][1],xmin,ymin]
								elif(Plist[n][3]>=ymin and Plist[n][3]<ymax and Plist[n][2]<=xmin):	#  = |
									temp=[0,Plist[n][1],xmin,Plist[n][3]]
								Plistn.append(temp)
								flagc=0
								break
					if(Nlist[Glist[contgates[i]-1][j]][m]==uncontgates[k]):
						#print("Gate",contgates[i],"is connected to Gate",uncontgates[k],"via Net",Glist[contgates[i]-1][j]) 
						if(Y[uncontgates[k]-1][0]>=ymax and X[uncontgates[k]-1][0]<xmin):			#  __|
							temp=[0,Glist[contgates[i]-1][j],xmin,ymax]
						elif(Y[uncontgates[k]-1][0]>=ymax and X[uncontgates[k]-1][0]>=xmin and X[uncontgates[k]-1][0]<xmax):																# |__|
							temp=[0,Glist[contgates[i]-1][j],X[uncontgates[k]-1][0],ymax]
						elif(Y[uncontgates[k]-1][0]>=ymax and X[uncontgates[k]-1][0]>=xmax):			# |__
							temp=[0,Glist[contgates[i]-1][j],xmax,ymax]
						elif(Y[uncontgates[k]-1][0]<ymax and Y[uncontgates[k]-1][0]>=ymin and X[uncontgates[k]-1][0]>=xmax):																# | =
							temp=[0,Glist[contgates[i]-1][j],xmax,Y[uncontgates[k]-1][0]]
						elif(Y[uncontgates[k]-1][0]<ymin and X[uncontgates[k]-1][0]>=xmax):			# |""
							temp=[0,Glist[contgates[i]-1][j],xmax,ymin] 
						elif(Y[uncontgates[k]-1][0]<=ymin and X[uncontgates[k]-1][0]>=xmin and X[uncontgates[k]-1][0]<xmax):																# |""|
							temp=[0,Glist[contgates[i]-1][j],X[uncontgates[k]-1][0],ymin]
						elif(Y[uncontgates[k]-1][0]<ymin and X[uncontgates[k]-1][0]<xmin):			#  ""|
							temp=[0,Glist[contgates[i]-1][j],xmin,ymin]
						elif(Y[uncontgates[k]-1][0]>=ymin and Y[uncontgates[k]-1][0]<ymax and X[uncontgates[k]-1][0]<=xmin):																#  = |
							temp=[0,Glist[contgates[i]-1][j],xmin,Y[uncontgates[k]-1][0]]
						else:
							if(shift==1):
								temp=[0,Glist[contgates[i]-1][j],xmax,Y[uncontgates[k]-1][0]]
							elif(shift==2):
								temp=[0,Glist[contgates[i]-1][j],xmin,Y[uncontgates[k]-1][0]]
							elif(shift==3):
								temp=[0,Glist[contgates[i]-1][j],X[uncontgates[k]-1][0],ymin]
							elif(shift==4):
								temp=[0,Glist[contgates[i]-1][j],X[uncontgates[k]-1][0],ymax]
						Plistn.append(temp)
						break
	return Plistn
	
##################################################################

#Bubble Sort Function to Sort Variables with Key 100000*X+Y
def bubble_sort(sort,X,Y):
	for i in range(len(sort)-1):
		for j in range(len(sort)-1-i):
			if((100000*X[sort[j]-1]+Y[sort[j]-1])>(100000*X[sort[j+1]-1]+Y[sort[j+1]-1])):
				temp=sort[j]
				sort[j]=sort[j+1]
				sort[j+1]=temp
	return sort

##################################################################

#Initialising Pads For Containment
def init_pads(Plist_append,Plist,xmin,xmax,ymin,ymax):
	for i in range(1,len(Plist)):
		if(Plist[i][2]<=xmax and Plist[i][2]>=xmin and Plist[i][3]>=ymin and Plist[i][3]<=ymax):
			Plist_append.append(Plist[i])
	return Plist_append
	
##################################################################

#Function to Create Uncontained Gate List
def uncont_list(contgates,Glist,gates):
	k=0
	uncontgates=np.zeros(gates-len(contgates),dtype=int)
	uncontgates=uncontgates.tolist()
	for i in range(len(Glist)):
		flag=1
		for j in range(len(contgates)):
			if(Glist[i][1]==contgates[j]):
				flag=1
				break
			else:
				flag=0
		if(flag==0):
			uncontgates[k]=Glist[i][1]
			k+=1
	return uncontgates

##################################################################

#Recursive Function to Partition the Gates
# Data Needed => gate array(sort), X, Y, xlb, xub, ylb, yub, Plist, Nlist, Glist, nets, gates

def recursive(Nlist,Glist,Plist,gates,nets,xlb,xub,ylb,yub,sort,X,Y):
	##################################################################
	#        !!!!!!!!!!        Starting Vertical Cut        !!!!!!!!!!
	##################################################################
	
	global rec
	rec=rec+1
	sort=bubble_sort(sort,X,Y)

	#Assigning Left and Right Sets of Gates
	leftg=np.zeros(int((len(sort))/2),dtype=int)
	if((len(sort))%2==0):
		rightg=np.zeros((int((len(sort))/2)),dtype=int)
	else:
		rightg=np.zeros((int((len(sort))/2)+1),dtype=int)
	for i in range(int((len(sort))/2)):
		leftg[i]=sort[i]

	for i in range(int((len(sort))/2),(len(sort))):
		rightg[i-int((len(sort))/2)]=sort[i]
	print("Left:\n ", leftg)
	print("Right:\n ", rightg)
	print()

	#Initialising Pads for Left Containment
	Plistl=np.zeros(1,dtype=int)
	Plistl=Plistl.tolist()
	Plistl=init_pads(Plistl,Plist,xlb,(xlb+xub)/2,ylb,yub)
	'''print("Pads for Left Containment")
	print(Plistl)'''

	#Calling Containment Function for Left Side
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax):
	Plistl=Containment(Plistl,Nlist,Glist,Plist,leftg,rightg,xlb,(xlb+xub)/2,ylb,yub,1)
	print("!!	Pads for Left Containment  :")
	print(Plistl)
	print()

	#Forming Net List for Left-Side Gates
	Nlist_left=np.zeros([nets+1,1],dtype=int)
	Nlist_left=Nlist_left.tolist()
	for i in range(len(leftg)):
		for j in range(Glist[leftg[i]-1][2]):
			Nlist_left[Glist[leftg[i]-1][j+3]].append(Glist[leftg[i]-1][1])
	'''print("Net-List for Left-Side Gates:")
	print(Nlist_left)
	print()'''

	#Forming Modified Net List for Left-Side Gates
	for i in range(1,len(Plistl)):
		Nlist_left[Plistl[i][1]][0]=1
	'''print("Modified Net-List for Left-Side Gates:")
	print(Nlist_left)
	print()'''

	#Calling X_Y Function
	X_left,Y_left=X_Y(Nlist_left,Plistl)

	#Updating X and Y values
	for i in range(len(leftg)):
		X[leftg[i]-1][0]=X_left[leftg[i]-1][0]
		Y[leftg[i]-1][0]=Y_left[leftg[i]-1][0]

	'''print("Updated X:")
	print(X)
	print()
	print("Updated Y:")
	print(Y)
	print()'''

	#Initialising Pads for Right Containment
	Plistr=np.zeros(1,dtype=int)
	Plistr=Plistr.tolist()
	Plistr=init_pads(Plistr,Plist,(xlb+xub)/2,xub,ylb,yub)
	'''print("Pads for Right Containment")
	print(Plistr)
	print()'''

	#Calling Containment Function for Right Side
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax):
	Plistr=Containment(Plistr,Nlist,Glist,Plist,rightg,leftg,(xlb+xub)/2,xub,ylb,yub,2)
	print("!!	Pads for Right Containment  :")
	print(Plistr)
	print()

	#Forming Net List for Right-Side Gates
	Nlist_right=np.zeros([nets+1,1],dtype=int)
	Nlist_right=Nlist_right.tolist()
	for i in range(len(rightg)):
		for j in range(Glist[rightg[i]-1][2]):
			Nlist_right[Glist[rightg[i]-1][j+3]].append(Glist[rightg[i]-1][1])

	#Forming Modified Net List for Right-Side Gates
	for i in range(1,len(Plistr)):
		Nlist_right[Plistr[i][1]][0]=1
	'''print("Modified Net-List for Right-Side Gates:")
	print(Nlist_right)
	print()'''

	#Calling X_Y Function
	X_right,Y_right=X_Y(Nlist_right,Plistr)

	#Updating X and Y values
	for i in range(len(rightg)):
		X[rightg[i]-1][0]=X_right[rightg[i]-1][0]
		Y[rightg[i]-1][0]=Y_right[rightg[i]-1][0]
	print("Final X:")
	print(X)
	print()
	print("Final Y:")
	print(Y)
	print()

	##################################################################
	#     !!!!!!!!!!        Starting Horizontal Cuts        !!!!!!!!!!
	##################################################################

	#Sorting 2(Left)
	sort_left=bubble_sort(leftg,X,Y)

	#Sorting 3(Right)
	sort_right=bubble_sort(rightg,X,Y)

	'''print(sort_left)
	print(sort_right)'''

	#Assigning Top and Bottom Sets of Gates on Left Side [b1,b2]
	b1=np.zeros((int(len(sort_left)/2)),dtype=int)
	if(len(sort_left)%2==0):
		b2=np.zeros((int(len(sort_left)/2)),dtype=int)
	else:
		b2=np.zeros((int(len(sort_left)/2)+1),dtype=int)

	for i in range(int(len(sort_left)/2)):
		b1[i]=sort_left[i]

	for i in range(int(len(sort_left)/2),len(sort_left)):
		b2[i-int(len(sort_left)/2)]=sort_left[i]

	'''print(b1)
	print(b2)'''

	#Assigning Top and Bottom Sets of Gates on Right Side [b3,b4]
	b3=np.zeros((int(len(sort_right)/2)),dtype=int)
	if(len(sort_right)%2==0):
		b4=np.zeros((int(len(sort_right)/2)),dtype=int)
	else:
		b4=np.zeros((int(len(sort_right)/2)+1),dtype=int)

	for i in range(int(len(sort_right)/2)):
		b3[i]=sort_right[i]

	for i in range(int(len(sort_right)/2),len(sort_right)):
		b4[i-int(len(sort_right)/2)]=sort_right[i]

	'''print(b3)
	print(b4)'''

	#	!!!!!!!!!!!!!!!         Starting B1         !!!!!!!!!!!!!!!

	#Initialising Pads for B1 Containment
	Plistb1=np.zeros(1,dtype=int)
	Plistb1=Plistb1.tolist()
	Plistb1=init_pads(Plistb1,Plist,xlb,(xlb+xub)/2,(ylb+yub)/2,yub)
	'''print("Pads for B1 Containment")
	print(Plistb1)'''

	#Creating Uncontained Gate List for B1
	#def uncont_list(contgates,Glist,gates):
	uncontb1=uncont_list(b1,Glist,gates)

	#Calling Containment Function for B1
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax):
	Plistb1=Containment(Plistb1,Nlist,Glist,Plist,b1,uncontb1,xlb,(xlb+xub)/2,(ylb+yub)/2,yub,3)
	'''print("!!	Pads for B1 Containment  :")
	print(Plistb1)
	print()'''

	#Forming Net List for b1 Gates
	Nlist_b1=np.zeros([nets+1,1],dtype=int)
	Nlist_b1=Nlist_b1.tolist()
	for i in range(len(b1)):
		for j in range(Glist[b1[i]-1][2]):
			Nlist_b1[Glist[b1[i]-1][j+3]].append(Glist[b1[i]-1][1])

	#Forming Modified Net List for b1 Gates
	for i in range(1,len(Plistb1)):
		Nlist_b1[Plistb1[i][1]][0]=1
	'''print("Modified Net-List for b1 Gates:")
	print(Nlist_b1)
	print()'''

	#Calling X_Y Function
	X_b1,Y_b1=X_Y(Nlist_b1,Plistb1)

	#Updating X and Y values
	for i in range(len(b1)):
		X[b1[i]-1][0]=X_b1[b1[i]-1][0]
		Y[b1[i]-1][0]=Y_b1[b1[i]-1][0]
		
	#	!!!!!!!!!!!!!!!         Starting B2         !!!!!!!!!!!!!!!

	#Initialising Pads for B2 Containment
	Plistb2=np.zeros(1,dtype=int)
	Plistb2=Plistb2.tolist()
	Plistb2=init_pads(Plistb2,Plist,xlb,(xlb+xub)/2,ylb,(ylb+yub)/2)
	'''print("Pads for B2 Containment")
	print(Plistb2)'''

	#Creating Uncontained Gate List for B2
	#def uncont_list(contgates,Glist,gates):
	uncontb2=uncont_list(b2,Glist,gates)

	#Calling Containment Function for B2
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax):
	Plistb2=Containment(Plistb2,Nlist,Glist,Plist,b2,uncontb2,xlb,(xlb+xub)/2,ylb,(ylb+yub)/2,4)
	'''print("!!	Pads for B2 Containment  :")
	print(Plistb2)
	print()'''

	#Forming Net List for b2 Gates
	Nlist_b2=np.zeros([nets+1,1],dtype=int)
	Nlist_b2=Nlist_b2.tolist()
	for i in range(len(b2)):
		for j in range(Glist[b2[i]-1][2]):
			Nlist_b2[Glist[b2[i]-1][j+3]].append(Glist[b2[i]-1][1])

	#Forming Modified Net List for b2 Gates
	for i in range(1,len(Plistb2)):
		Nlist_b2[Plistb2[i][1]][0]=1
	'''print("Modified Net-List for b2 Gates:")
	print(Nlist_b2)
	print()'''

	#Calling X_Y Function
	X_b2,Y_b2=X_Y(Nlist_b2,Plistb2)

	#Updating X and Y values
	for i in range(len(b2)):
		X[b2[i]-1][0]=X_b2[b2[i]-1][0]
		Y[b2[i]-1][0]=Y_b2[b2[i]-1][0]

	#	!!!!!!!!!!!!!!!         Starting B3         !!!!!!!!!!!!!!!

	#Initialising Pads for B3 Containment
	Plistb3=np.zeros(1,dtype=int)
	Plistb3=Plistb3.tolist()
	Plistb3=init_pads(Plistb3,Plist,(xlb+xub)/2,xub,(ylb+yub)/2,yub)
	'''print("Pads for B3 Containment")
	print(Plistb3)'''

	#Creating Uncontained Gate List for B3
	#def uncont_list(contgates,Glist,gates):
	uncontb3=uncont_list(b3,Glist,gates)

	#Calling Containment Function for B3
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax,shift,xshiftval,yshiftval):
	Plistb3=Containment(Plistb3,Nlist,Glist,Plist,b3,uncontb3,(xlb+xub)/2,xub,(ylb+yub)/2,yub,3)
	'''print("!!	Pads for B3 Containment  :")
	print(Plistb3)
	print()'''

	#Forming Net List for b3 Gates
	Nlist_b3=np.zeros([nets+1,1],dtype=int)
	Nlist_b3=Nlist_b3.tolist()
	for i in range(len(b3)):
		for j in range(Glist[b3[i]-1][2]):
			Nlist_b3[Glist[b3[i]-1][j+3]].append(Glist[b3[i]-1][1])

	#Forming Modified Net List for b3 Gates
	for i in range(1,len(Plistb3)):
		Nlist_b3[Plistb3[i][1]][0]=1
	'''print("Modified Net-List for b3 Gates:")
	print(Nlist_b3)
	print()'''

	#Calling X_Y Function
	X_b3,Y_b3=X_Y(Nlist_b3,Plistb3)

	#Updating X and Y values
	for i in range(len(b3)):
		X[b3[i]-1][0]=X_b3[b3[i]-1][0]
		Y[b3[i]-1][0]=Y_b3[b3[i]-1][0]

	#	!!!!!!!!!!!!!!!         Starting B4         !!!!!!!!!!!!!!!

	#Initialising Pads for B4 Containment
	Plistb4=np.zeros(1,dtype=int)
	Plistb4=Plistb4.tolist()
	Plistb4=init_pads(Plistb4,Plist,(xlb+xub)/2,xub,ylb,(ylb+yub)/2)
	'''print("Pads for B4 Containment")
	print(Plistb4)'''

	#Creating Uncontained Gate List for B4
	#def uncont_list(contgates,Glist,gates):
	uncontb4=uncont_list(b4,Glist,gates)

	#Calling Containment Function for B4
	#def Containment(Plistn,Nlist,Glist,Plist,contgates,uncontgates,xmin,xmax,ymin,ymax):
	Plistb4=Containment(Plistb4,Nlist,Glist,Plist,b4,uncontb4,(xlb+xub)/2,xub,ylb,(ylb+yub)/2,4)
	'''print("!!	Pads for B4 Containment  :")
	print(Plistb4)
	print()'''

	#Forming Net List for b4 Gates
	Nlist_b4=np.zeros([nets+1,1],dtype=int)
	Nlist_b4=Nlist_b4.tolist()
	for i in range(len(b4)):
		for j in range(Glist[b4[i]-1][2]):
			Nlist_b4[Glist[b4[i]-1][j+3]].append(Glist[b4[i]-1][1])

	#Forming Modified Net List for b4 Gates
	for i in range(1,len(Plistb4)):
		Nlist_b4[Plistb4[i][1]][0]=1
	'''print("Modified Net-List for b4 Gates:")
	print(Nlist_b4)
	print()'''

	#Calling X_Y Function
	X_b4,Y_b4=X_Y(Nlist_b4,Plistb4)

	#Updating X and Y values
	for i in range(len(b4)):
		X[b4[i]-1][0]=X_b4[b4[i]-1][0]
		Y[b4[i]-1][0]=Y_b4[b4[i]-1][0]
	
	#def recursive(Nlist,Glist,Plist,gates,nets,xlb,xub,ylb,yub,sort,X,Y):
	if(rec<=2):	
		X,Y=recursive(Nlist,Glist,Plist,gates,nets,xlb,(xlb+xub)/2,(ylb+yub)/2,yub,b1,X,Y)
		X,Y=recursive(Nlist,Glist,Plist,gates,nets,xlb,(xlb+xub)/2,ylb,(ylb+yub)/2,b2,X,Y)
		X,Y=recursive(Nlist,Glist,Plist,gates,nets,(xlb+xub)/2,xub,(ylb+yub)/2,yub,b3,X,Y)
		X,Y=recursive(Nlist,Glist,Plist,gates,nets,(xlb+xub)/2,xub,ylb,(ylb+yub)/2,b4,X,Y)
	else:
		return(X,Y)
	rec=1
	return(X,Y)
	

##########################################################################################################################################

# Main Function
f=open("data1.txt","r")
a=f.readline()
l=a.split()

gates=int(l[0])
nets=int(l[1])

#Forming Net and Gate List
Nlist=np.zeros([nets+1,1],dtype=int)
Glist=np.zeros([gates,1],dtype=int)
Nlist=Nlist.tolist()
Glist=Glist.tolist()
k=0;
for i in range(gates):
	a=f.readline()
	l=a.split()
	for k in range(len(l)):
		Glist[i].append(int(l[k]))
	for j in range(int(l[1])):
		if(Nlist[int(l[j+2])]!=None):
			Nlist[int(l[j+2])].append(int(l[0]))
		else:
			Nlist[int(l[j+2])]=int(l[0])

#Forming Pad List
a=f.readline()
Plist=np.zeros([int(a)+1,1],dtype=int)
Plist=Plist.tolist()
for i in range(int(a)):
	a=f.readline()
	l=a.split()
	Plist[int(l[0])].append(int(l[1]))
	Plist[int(l[0])].append(int(l[2]))
	Plist[int(l[0])].append(int(l[3]))
f.close()

#Forming Modified Net List
for i in range(1,len(Plist)):
	Nlist[Plist[i][1]][0]=1

print("Net List:")
print(Nlist)
print()
print("Pad List:")
print(Plist)
print()
print("Gate List: ")
print(Glist)
print()
	
#Calling X_Y Function
X,Y=X_Y(Nlist,Plist)

#Sorting 1(Entire Set)
sort=np.zeros(gates,dtype=int)
for i in range(gates):
	sort[i]=i+1

rec=0
#Call Recursion
#def recursive(Nlist,Glist,Plist,gates,nets,xlb,xub,ylb,yub,sort,X,Y):
X,Y=recursive(Nlist,Glist,Plist,gates,nets,0,100,0,100,sort,X,Y)

print("Final X:")
print(X)
print()
print("Final Y:")
print(Y)
print()

#Writing into file output.txt
o=open("output.txt","w")
for i in range(len(X)):
	o.write(str(i+1))
	o.write(" ")
	o.write(str(X[i][0]))
	o.write(" ")
	o.write(str(Y[i][0]))
	o.write("\n")

print("-------------------Done-------------------")
