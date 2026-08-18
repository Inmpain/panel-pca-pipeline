#!/bin/bash
# Union two newline-delimited SNP ID lists (e.g. backbone ∪ covered), dedup,
# stable via sort -u. Takes column 1 of each line so it tolerates extra columns.
# Panel-agnostic: works for 720 hybrid and later 3K.
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "usage: $0 LIST_A LIST_B OUT" >&2
  exit 2
fi
A="$1"; B="$2"; OUT="$3"
[[ -s "$A" ]] || { echo "FATAL: missing/empty $A" >&2; exit 1; }
[[ -s "$B" ]] || { echo "FATAL: missing/empty $B" >&2; exit 1; }

awk 'NF{print $1}' "$A" "$B" | sort -u > "$OUT"
echo "union -> $OUT: $(wc -l < "$OUT") markers (A=$(wc -l < "$A"), B=$(wc -l < "$B"))"
