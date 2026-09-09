import sys
from os.path import splitext, split, join


class Entry:
	def __init__(self, line):
		line = line.strip()
		seq_name, source, type_, start, end, score, strand, phase, attributes = line.split('\t')
		self.seq_name = seq_name
		self.source = source
		self.type = type_
		self.start = int(start)
		self.end = int(end)
		try:
			self.score = float(score)
		except ValueError:
			self.score = 0.0
		self.strand = strand
		self.phase = phase
		self.attributes = attributes
		self.seq_id = self.attributes[self.attributes.find("=") + 1:self.attributes.find(";")]

	def __str__(self):
		return "\t".join([
			self.seq_name,
			self.source,
			self.type,
			str(self.start),
			str(self.end),
			f"{self.score:g}" if self.score != 0.0 else ".",
			self.strand,
			self.phase,
			self.attributes
		])

	def extend_with(self, other):
		self.end = max(self.end, other.end)
		self.score = other.score
		self.attributes = other.attributes
		self.seq_id = other.seq_id

	def sort_key(self):
		return self.strand, self.seq_name, self.start, self.seq_id, 0 if self.type == "gene" else 1

	def __lt__(self, other):
		return self.sort_key() < other.sort_key()


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Input file is missing")
		exit(1)
	data_path = sys.argv[1]
	to_remove = set()
	if len(sys.argv) >= 3:
		remove_path = sys.argv[2]
		with open(remove_path, "r") as remove_file:
			to_remove = set(map(str.strip, remove_file.readlines()))
	base_path, file_path = split(data_path)
	file_name, file_ext = splitext(file_path)
	with open(data_path, "r") as data_file:
		entries = list(filter(lambda e: e.seq_id not in to_remove, map(Entry, data_file.readlines())))
	id_count = {}
	for entry in entries:
		if entry.seq_id in id_count:
			id_count[entry.seq_id] += 1
		else:
			id_count[entry.seq_id] = 1
	monoexonic_entries = []
	non_monoexonic_entries = []
	for entry in entries:
		# NB: plain if/elif rather than a `match` statement (Python 3.10+ syntax) so this
		# runs on whatever python3 happens to be on the execution host, not just newer ones.
		count = id_count[entry.seq_id]
		if count == 1:
			print(f"WARNING: id {entry.seq_id} appears only once in the file")
		elif count == 2:
			monoexonic_entries.append(entry)
		else:
			non_monoexonic_entries.append(entry)

	# scrittura del file con righe monoexonic
	with open(join(base_path, f"{file_name}_monoexonic.gff3"), "w") as monoexonic_file:
		monoexonic_file.writelines(map(lambda e: f"{e}\n", monoexonic_entries))
	# scrittura del file con righe non monoexonic
	with open(join(base_path, f"{file_name}_nonmonoexonic.gff3"), "w") as non_monoexonic_file:
		non_monoexonic_file.writelines(map(lambda e: f"{e}\n", non_monoexonic_entries))
	monoexonic_entries.sort()
	# scrittura del file con righe monoexonic ordinate
	with open(join(base_path, f"{file_name}_monoexonic_sorted.gff3"), "w") as monoexonic_file:
		monoexonic_file.writelines(map(lambda e: f"{e}\n", monoexonic_entries))

	non_redundant_monoexonic_entries = []
	if len(monoexonic_entries) > 0:
		non_redundant_monoexonic_entries.extend([monoexonic_entries[0], monoexonic_entries[1]])
		i = 2
		while i < len(monoexonic_entries):
			entry = monoexonic_entries[i]
			last_entry = non_redundant_monoexonic_entries[-2]
			if entry.seq_name == last_entry.seq_name and last_entry.start <= entry.start <= last_entry.end:
				last_entry.extend_with(entry)
				non_redundant_monoexonic_entries[-1].extend_with(monoexonic_entries[i + 1])
			else:
				non_redundant_monoexonic_entries.append(entry)
				non_redundant_monoexonic_entries.append(monoexonic_entries[i + 1])
			i += 2
	else:
		print("WARNING: found 0 monoexonic entries")
	non_redundant_monoexonic_entries.sort()
	with open(join(base_path, f"{file_name}_monoexonic_NR.gff3"), "w") as non_redundant_monoexonic_file:
		non_redundant_monoexonic_file.writelines(map(lambda e: f"{e}\n", non_redundant_monoexonic_entries))

	result = non_redundant_monoexonic_entries + non_monoexonic_entries
	result.sort(key=lambda e: (e.seq_id, 0 if e.type == "gene" else 1, e.start))
	with open(join(base_path, f"{file_name}_NR.gff3"), "w") as result_file:
		result_file.writelines(map(lambda e: f"{e}\n", result))
