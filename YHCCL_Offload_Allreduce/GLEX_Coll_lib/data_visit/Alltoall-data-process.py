name="2-512.out"
f=open(name)
MPI_Vec = []
GLEX_Vec = []
for noden in range(0,9):
	mpitmpv=[]
	glextmpv=[]
	TEXT=f.readline()
	for size in range(0,20):
		tmpv=f.readline().split()
		mpitmpv.append(float(tmpv[1]))
		glextmpv.append(float(tmpv[3]))
	MPI_Vec.append(mpitmpv)
	GLEX_Vec.append(glextmpv)
f.close()

for a in MPI_Vec:
	for b in a:
		print(b,end=' ')
	print("")
for a in GLEX_Vec:
	for b in a:
		print(b,end=' ')
	print("")
