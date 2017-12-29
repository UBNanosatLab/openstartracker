import sys
stars=open("hip_main.dat","r").readlines()
result=open(sys.argv[1],"r").readlines()
result_real=open(sys.argv[2],"r").readlines()
assert len(result)==len(result_real)

stardict={}
for i in range(0,len(stars)):
	stardict[int(stars[i].split("|")[1])]=i

scores=[]
for i in range(0,len(result)):
	thisscore=[]
	t=0
	c=0
	w=0
	a=result[i].split(",")
	b=result_real[i].split(",")
	assert len(a)==len(b)
	for j in range(0,len(a)):
		if int(a[j])!=-1:
			t+=1
			thisscore+=[stardict[int(a[j])]]
		if int(b[j])!=-1:
			if int(a[j])==int(b[j]):
				c+=1
			else:
				w+=1
	thisscore=[max((c-2*w)/float(max(t,1)),-1)]+thisscore
	if t>0:
		scores+=[[1.0]+thisscore]
	else:
		scores+=[[0.0]+thisscore]

a=0
b=0		
for i in range(0,len(scores)):
	if scores[i][1]<scores[i][0]:
		print i," ".join([str(j) for j in scores[i]])
	a+=scores[i][0]
	b+=scores[i][1]
	
print a,b
