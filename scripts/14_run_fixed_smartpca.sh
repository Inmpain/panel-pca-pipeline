#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --help || ${1:-} == -h ]]; then echo "usage: $0 --config CFG --geno GENO --snp SNP --ind IND --poplist POP --label NAME --out-dir DIR"; exit 0; fi
CFG= GENO= SNP= IND= POP= LABEL= OUT=
while [[ $# -gt 0 ]]; do case "$1" in --config) CFG=$2; shift 2;; --geno) GENO=$2; shift 2;; --snp) SNP=$2; shift 2;; --ind) IND=$2; shift 2;; --poplist) POP=$2; shift 2;; --label) LABEL=$2; shift 2;; --out-dir) OUT=$2; shift 2;; *) echo "unknown argument: $1" >&2; exit 2;; esac; done
for x in CFG GENO SNP IND POP LABEL OUT; do [[ -n ${!x} ]] || { echo "missing argument" >&2; exit 2; }; done
command -v smartpca >/dev/null || { echo "smartpca not found" >&2; exit 1; }; mkdir -p "$OUT"; P="$OUT/$LABEL"
python3 - "$CFG" "$P.par" "$GENO" "$SNP" "$IND" "$POP" "$P" <<'PY'
import sys,yaml
c=yaml.safe_load(open(sys.argv[1]))['pca']; out,geno,snp,ind,pop,prefix=sys.argv[2:]
with open(out,'x') as f:
 for k,v in [('genotypename',geno),('snpname',snp),('indivname',ind),('evecoutname',prefix+'.evec'),('evaloutname',prefix+'.eval'),('poplistname',pop),('lsqproject','YES'),('numoutevec',c['num_pcs']),('numoutlieriter',c['numoutlieriter']),('numchrom',c['numchrom']),('numthreads',c['numthreads'])]: f.write(f'{k}: {v}\n')
PY
smartpca -p "$P.par" > "$P.smartpca.log" 2>&1; test -s "$P.evec"; test -s "$P.eval"; echo "PASS: $P"
