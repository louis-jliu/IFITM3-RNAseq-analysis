"""00_build_gene_program_table.py

Author : Jinrui (Louis) Liu (Computational and Systems Biology Program, MSKCC).
         Contact: liuj11@mskcc.org

Build data/gene_programs/microglia_programs_Gan2024.tsv, the microglial gene programs used
for Fig. 5e, from Table S1 ("Reference gene sets" sheet) of Gan et al., Cancer Cell 42,
1693-1712 (2024), doi:10.1016/j.ccell.2024.08.015 (supplementary file mmc2.xlsx).

The TSV is already included in the repository; this script documents how it was made.

Programs used (component gene sets are merged, in this order):
  ATM               : ATM
  IRM               : IRM
  TGFbeta-dependent : TGF-beta-dependent (in vitro) down, (in vitro) up, (in vivo) down, (in vivo) up
  DAM-program       : DAM (core program), DAM
  LRRC33-dependent  : LRRC33-dependent upregulated, LRRC33-dependent downregulated
Mouse symbols are converted to human by upper-casing.

Usage (from the repository root):
  python scripts/00_build_gene_program_table.py data/gene_programs/Gan2024_TableS1.xlsx
"""

import csv
import re
import sys

import openpyxl

PROGRAMS = [
    ("ATM", [r"^ATM$"]),
    ("IRM", [r"^IRM$"]),
    ("TGFbeta-dependent", [r"^TGF.*-dependent \(in vitro\) downregulated$",
                           r"^TGF.*-dependent \(in vitro\) upregulated$",
                           r"^TGF.*-dependent \(in vivo\) downregulated$",
                           r"^TGF.*-dependent \(in vivo\) upregulated$"]),
    ("DAM-program", [r"^DAM \(core program\)$", r"^DAM$"]),
    ("LRRC33-dependent", [r"^LRRC33-dependent upregulated$", r"^LRRC33-dependent downregulated$"]),
]


def main(xlsx, out="data/gene_programs/microglia_programs_Gan2024.tsv"):
    ws = openpyxl.load_workbook(xlsx, read_only=True).worksheets[1]
    rows = list(ws.iter_rows(values_only=True))
    header_idx = next(i for i, r in enumerate(rows) if r and r[0] == "Gene set name abbreviated")
    sets = {}
    for r in rows[header_idx + 1:]:
        if r and r[0]:
            sets[r[0]] = [str(g) for g in r[3:] if g not in (None, "")]  # cols 1-3: name, full name, source

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["program", "source_gene_set", "gene_mouse", "gene"])
        for program, patterns in PROGRAMS:
            for pat in patterns:
                hits = [s for s in sets if re.match(pat, s)]
                if len(hits) != 1:
                    sys.exit(f"pattern {pat!r} matched {hits}")
                for g in sets[hits[0]]:
                    w.writerow([program, hits[0], g, g.upper()])
    print("wrote", out)


if __name__ == "__main__":
    main(*sys.argv[1:])
