import sys
import re

# Combines the per-step BUSCO short_summary.txt files (each named "<label>.short_summary.txt"
# by buscoassess.nf) into one TSV table, one row per pipeline step.
#
# NB: I couldn't run busco here to confirm the exact wording/spacing of its short-summary text
# for v6.1.0 - this is written against the standard, long-stable BUSCO one-line summary format:
#   C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:978
# and the "Lineage dataset is: <name> (...)" line. If parsing comes back empty for a real run,
# send me one of the short_summary.txt files and I'll fix the regex.

if len(sys.argv) < 2:
    print("usage: python3 busco_summary.py <label1.short_summary.txt> [label2.short_summary.txt ...]")
    exit(1)

SCORE_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)

rows = []
for path in sys.argv[1:]:
    label = path[:-len(".short_summary.txt")] if path.endswith(".short_summary.txt") else path
    lineage = ""
    match = None
    with open(path) as summary_file:
        for line in summary_file:
            line = line.strip()
            if line.startswith("C:"):
                match = SCORE_RE.search(line)
            if "Lineage dataset is:" in line:
                lineage = line.split("Lineage dataset is:", 1)[1].split("(")[0].strip()
    if not match:
        print(f"WARNING: could not parse a BUSCO one-line summary out of {path}")
        continue
    complete, single, dup, frag, missing, n_total = match.groups()
    rows.append((label, lineage, complete, single, dup, frag, missing, n_total))

with open("busco_summary.tsv", "w") as out:
    out.write("\t".join([
        "step", "lineage", "complete_pct", "single_copy_pct",
        "duplicated_pct", "fragmented_pct", "missing_pct", "n_total",
    ]) + "\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

print(f"Wrote busco_summary.tsv with {len(rows)} row(s)")
