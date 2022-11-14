import sys
#IN1: R1 read file; IN2: R2 read file
IN1=open(sys.argv[1])
IN2=open(sys.argv[2])
#index file: each line like "CGAGTAAC	soybean1"
index_file=open(sys.argv[3])
index={}
for line in index_file:
	tmp=line.strip().split('\t')
	index[tmp[0]]=open("./"+tmp[1]+'.fastq','w')
	
for R1 in IN1:
	R1=R1+IN1.next()+IN1.next()+IN1.next()
	R2=[IN2.readline(),IN2.readline(),IN2.readline(),IN2.readline()]
	if R2[1][:8] in index:
		index[R2[1][:8]].write(R1)

