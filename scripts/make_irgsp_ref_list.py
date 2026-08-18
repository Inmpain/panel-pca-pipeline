#!/usr/bin/env python3
"""Generate a 'SNP_ID<TAB>irgsp_base' reference-allele list for
`plink --a2-allele ... 2 1 --keep-allele-order --make-bed`.

For each biallelic panel site, the reference allele is the base present in
irgsp.fa (so for the clean-swap 3K/Civan panels this is the panel's ALT; for
the mixed 720 panel it is REF at ~91% of sites and ALT at ~9%). Mismatch sites
(FASTA base is a third allele or N) are skipped.
"""
import gzip
import sys


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_fasta(path):
    seqs = {}
    chrom = None
    chunks = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if chrom is not None:
            seqs[chrom] = "".join(chunks).upper()
    return seqs


def main(snp_path, fasta_path, out_path):
    seqs = load_fasta(fasta_path)
    n_ref = n_alt = n_skip = 0
    with open_maybe_gz(snp_path) as fh, open(out_path, "w") as out:
        for line in fh:
            f = line.split()
            if len(f) != 6:
                continue
            sid, chrom, _gp, pos_text, ref, alt = f
            ref = ref.upper()
            alt = alt.upper()
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            try:
                contig = "chr%02d" % int(chrom.lstrip("0") or "0")
            except ValueError:
                contig = chrom
            seq = seqs.get(contig)
            if seq is None:
                continue
            idx = pos - 1
            if idx < 0 or idx >= len(seq):
                continue
            fa = seq[idx]
            if fa == ref:
                out.write(f"{sid}\t{ref}\n")
                n_ref += 1
            elif fa == alt:
                out.write(f"{sid}\t{alt}\n")
                n_alt += 1
            else:
                n_skip += 1
    print(f"wrote {out_path}: ref_as_is={n_ref} ref_flip_to_alt={n_alt} skipped_mismatch={n_skip}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python3 make_irgsp_ref_list.py SNP FASTA OUT", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
