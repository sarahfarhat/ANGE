import sys, os

if len(sys.argv) < 3:
	print("usage: phyton3 extractseq_genomev4.py <genome.tab> <genes.gff>")
	exit(1)

genome_path = sys.argv[1]
gff_path = sys.argv[2]

print("Building HashMap for genome...")
genome_dict = {}
with open(genome_path, "r") as genome_file:
	for i in genome_file.readlines():
	  	genome_dict.update({i.split('\t')[0] : i.split('\t')[1]})

with open(gff_path, "r") as gff_file, open(os.path.splitext(sys.argv[2])[0]+".fasta", "x") as output:
	lines=list(gff_file.readlines())
	print("Extracting sequences... (total: "+str(len(lines))+")")
	n=0
	for i in lines:
		genomekey = i.split('\t')[0]
		gene = i.split('\t')[2]
		startpos = int(i.split('\t')[3])
		endpos = int(i.split('\t')[4])
		match = i.split('\t')[8].split('=')[1].rstrip('\n')

		if gene in ("HSP", "gene", "mRNA", "exon"):
			output.write(">"+match+"_"+str(n)+"\n")
			output.write(genome_dict.get(genomekey)[startpos-1:endpos]+"\n")
			n+=1

		if n<10000 and n%1000==0:
			print(str(n)+" sequences done")
		if n<100000 and n%10000==0:
			print(str(n)+" sequences done")
		elif n%100000==0:
			print(str(n)+" sequences done")

print("Analysis complete")
