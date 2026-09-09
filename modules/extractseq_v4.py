import sys, os

if len(sys.argv) < 4:
	print("usage: phyton3 extractseq_genomev4.py <genome.tab> <query.tab> <genes.gff>")
	exit(1)

genome_path = sys.argv[1]
query_path = sys.argv[2]
gff_path = sys.argv[3]

os.mkdir(os.path.splitext(sys.argv[3])[0]+"_query")
os.mkdir(os.path.splitext(sys.argv[3])[0]+"_genome")

print("Building HashMap for query dataset...")
query_dict = {}
with open(query_path, "r") as query_file:
	for i in query_file.readlines():
	  	query_dict.update({i.split('\t')[0] : i.split('\t')[1]})

print("Building HashMap for genome...")
genome_dict = {}
with open(genome_path, "r") as genome_file:
	for i in genome_file.readlines():
	  	genome_dict.update({i.split('\t')[0] : i.split('\t')[1]})

with open(gff_path, "r") as gff_file:
    lines=list(gff_file.readlines())
    n=0
    print("Extracting sequences... (total: "+str(len(lines))+")")
    for i in lines:
        genomekey = i.split('\t')[0]
        startpos = int(i.split('\t')[3])
        endpos = int(i.split('\t')[4])
        fastakey = i.split('\t')[8].split('=')[1].rstrip('\n')

        f = open(os.path.splitext(sys.argv[3])[0]+"_query/sequence_"+str(n)+".faa", "x")
        f.write(">"+fastakey+"\n")
        f.write(query_dict.get(fastakey)+"\n")

        f = open(os.path.splitext(sys.argv[3])[0]+"_genome/sequence_"+str(n)+".fa", "x")
        f.write(">"+genomekey+":"+str(startpos)+"-"+str(endpos)+"\n")
        f.write(genome_dict.get(genomekey)[startpos-1:endpos]+"\n")

        n+=1

        if n<10000 and n%1000==0:
            print(str(n)+" sequences done")
        if n<100000 and n%10000==0:
            print(str(n)+" sequences done")
        elif n%100000==0:
            print(str(n)+" sequences done")

print("Analysis complete")
