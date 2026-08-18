#!/usr/bin/env python3
"""Convert one sample's pileupCaller PLINK output to a .calls.txt.

Output is one 0/2/9 per line, in the sample .bim order. pileupCaller --randomHaploid
produces haploid calls, so 0=hom ALT, 2=hom REF, 9=missing (any '1'/NA also -> 9).
"""
import argparse
import subprocess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bfile", required=True, help="sample PLINK prefix (.bed/.bim/.fam)")
    ap.add_argument("--out", required=True, help="output prefix (writes OUT.calls.txt)")
    args = ap.parse_args()

    # Read the .bim marker order and A2 (REF) allele first, so we are robust to
    # extra .raw columns and can orient the plink2 counted-allele suffix back to
    # REF/ALT. pileupCaller convention (see pileupcaller_shared_call.sh): A2=REF.
    bim_ids = []
    bim_ref = {}  # SNP ID -> A2 (REF) allele
    with open(args.bfile + ".bim") as fh:
        for line in fh:
            f = line.split()
            sid = f[1]
            bim_ids.append(sid)
            bim_ref[sid] = f[5]

    raw = args.out + ".raw"
    subprocess.run(["plink2", "--bfile", args.bfile, "--export", "A", "--out", args.out],
                   check=True)

    with open(raw) as fh:
        header = fh.readline().rstrip("\n").split()
        line = fh.readline().rstrip("\n").split()

    # plink2 --export A names variant columns "<ID>_<countedAllele>"
    # (e.g. "1np2833_T"). Map stripped ID -> (data index, counted allele).
    col = {}
    for i, name in enumerate(header):
        if "_" in name:
            sid, allele = name.rsplit("_", 1)
            col[sid] = (i, allele)

    calls = []
    for sid in bim_ids:
        if sid not in col:
            calls.append("9")
            continue
        idx, counted = col[sid]
        val = line[idx]
        ref = bim_ref[sid]
        if val == "2":
            # Two copies of the counted allele -> homozygous for that allele.
            calls.append("2" if counted == ref else "0")
        elif val == "0":
            # Zero copies of the counted allele -> homozygous for the other allele.
            calls.append("0" if counted == ref else "2")
        else:
            calls.append("9")  # '1' (unexpected for haploid) or NA -> missing

    out_calls = args.out + ".calls.txt"
    with open(out_calls, "w") as fh:
        fh.write("\n".join(calls) + "\n")
    print(f"wrote {out_calls}: {len(calls)} markers")


if __name__ == "__main__":
    main()
