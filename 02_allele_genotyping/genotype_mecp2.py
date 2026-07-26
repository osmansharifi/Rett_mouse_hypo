#!/usr/bin/env python3
"""Track B: per-cell allele call at the Mecp2 e1 start codon.

Runs `samtools mpileup --output-extra CB,UB` at a single position, walks the
pileup string into one symbol per read, zips it against the per-read CB/UB tag
lists, collapses to one call per (barcode, UMI), and tallies WT vs MUT per
barcode. Emits every barcode seen at the locus (incl. ambient / no-CB reads as
an AMBIENT aggregate row), so filtered-cell calls and the soup baseline both
come out of one pass.

WT (plus strand) = T, MUT (ATG->TTG start loss) = A, at chrX:74085586.
"""
import argparse, csv, subprocess, sys
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--chrom", default="X")
ap.add_argument("--pos", type=int, default=74085586)
ap.add_argument("--wt", default="T")
ap.add_argument("--mut", default="A")
ap.add_argument("--out", required=True)
a = ap.parse_args()

region = f"{a.chrom}:{a.pos}-{a.pos}"
cmd = ["samtools", "mpileup", "-B", "-Q", "0", "-q", "0",
       "-r", region, "--output-extra", "CB,UB", a.bam]
res = subprocess.run(cmd, capture_output=True, text=True)
if res.returncode != 0:
    sys.stderr.write(res.stderr)
    sys.exit(res.returncode)

header = ["sample", "barcode", "wt_umi", "mut_umi",
          "other_umi", "total_umi", "vaf_mut", "call"]
lines = [l for l in res.stdout.splitlines() if l.strip()]
if not lines:
    with open(a.out, "w", newline="") as f:
        csv.writer(f).writerow(header)
    print(f"{a.sample}: no reads at {region}")
    sys.exit(0)

col   = lines[0].split("\t")
bases = col[4]
cb_list = col[6].split(",") if len(col) > 6 and col[6] != "" else []
ub_list = col[7].split(",") if len(col) > 7 and col[7] != "" else []

# --- walk the pileup string: one symbol per read, in read order ---
syms = []; i = 0; n = len(bases)
while i < n:
    c = bases[i]
    if c == "^":                       # start marker + mapq char; base follows
        i += 2; continue
    if c == "$":                       # end marker on the previous base
        i += 1; continue
    if c in "+-":                      # indel on previous read: skip len+seq
        j = i + 1; num = ""
        while j < n and bases[j].isdigit():
            num += bases[j]; j += 1
        i = j + (int(num) if num else 0); continue
    syms.append(c); i += 1             # base letter, *, <, or >

if not (len(syms) == len(cb_list) == len(ub_list)):
    sys.stderr.write(f"WARN {a.sample}: length mismatch "
                     f"syms={len(syms)} cb={len(cb_list)} ub={len(ub_list)}\n")
m = min(len(syms), len(cb_list), len(ub_list))

# --- bin reads: barcode -> umi -> Counter(base) ; ambient kept separately ---
bc_umi  = defaultdict(lambda: defaultdict(Counter))
amb_umi = defaultdict(Counter)
anon = 0
for k in range(m):
    s = syms[k].upper()
    if s not in "ACGTN":               # skip ref-skips (<>) and deletions (*)
        continue
    cb = cb_list[k]; ub = ub_list[k]
    cb = None if cb in ("", "*") else cb
    ub = None if ub in ("", "*") else ub
    if cb is None:
        key = ub if ub is not None else f"__anon{anon}"
        if ub is None: anon += 1
        amb_umi[key][s] += 1
    else:
        u = ub if ub is not None else f"__anon{anon}"
        if ub is None: anon += 1
        bc_umi[cb][u][s] += 1

def collapse(umis):
    """One consensus base per UMI, then a per-base UMI tally."""
    t = Counter()
    for _, cnt in umis.items():
        t[cnt.most_common(1)[0][0]] += 1
    return t

def row(name, tally):
    wt = tally.get(a.wt, 0); mut = tally.get(a.mut, 0)
    other = sum(v for b, v in tally.items() if b not in (a.wt, a.mut))
    total = wt + mut + other
    vaf = (mut / (wt + mut)) if (wt + mut) > 0 else ""
    if total == 0:            call = "NA"
    elif mut == 0 and wt > 0: call = "WT"
    elif wt == 0 and mut > 0: call = "MUT"
    else:                     call = "MIXED"
    if name == "AMBIENT":     call = "AMBIENT"
    return [a.sample, name, wt, mut, other, total, vaf, call]

rows = [row(cb, collapse(umis)) for cb, umis in bc_umi.items()]
amb  = row("AMBIENT", collapse(amb_umi))
rows.append(amb)

with open(a.out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(header); w.writerows(rows)

print(f"{a.sample}: {len(bc_umi)} barcodes genotyped at {region}; "
      f"ambient WT={amb[2]} MUT={amb[3]}")
