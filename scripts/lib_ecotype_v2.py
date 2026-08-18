"""Shared helpers for scripts/ecotype_pca_v2/. Config-driven, read-only unless
a script explicitly writes to results_v2_root. No script in this package may
touch results/ecotype_pca/ (the v1 tree) or run anything against ancient
sample intersections -- see docs/ECOTYPE_PCA_V2_SPEC.md.
"""
import argparse
import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


def load_config(path):
    if yaml is None:
        sys.exit("FATAL: pyyaml not importable. Install it (pip install pyyaml) "
                 "before running any ecotype_pca_v2 script.")
    p = Path(path)
    if not p.is_file():
        sys.exit(f"FATAL: config not found: {p}")
    with open(p) as fh:
        cfg = yaml.safe_load(fh)
    if cfg.get("version") != "ecotype_pca_v2":
        sys.exit(f"FATAL: config version mismatch in {p}: "
                  f"expected 'ecotype_pca_v2', got {cfg.get('version')!r}")
    return cfg


def base_argparser(description):
    ap = argparse.ArgumentParser(description=description,
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", required=True,
                     help="path to ecotype_pca_v2.yaml (single source of truth for all parameters)")
    ap.add_argument("--out-dir", required=True,
                     help="output directory for this script's products; must not already contain "
                          "this script's output files unless --overwrite is given")
    ap.add_argument("--overwrite", action="store_true",
                     help="allow overwriting pre-existing output files from a prior run of this script")
    return ap


def setup_logger(name, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / f"{name}.log"
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(log_path, mode="a")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger, log_path


def check_output_not_present(paths, overwrite, logger):
    existing = [str(p) for p in paths if Path(p).exists()]
    if existing and not overwrite:
        for p in existing:
            logger.error(f"refusing to overwrite existing output (pass --overwrite to allow): {p}")
        sys.exit(3)
    if existing and overwrite:
        for p in existing:
            logger.warning(f"--overwrite given, will replace: {p}")


def run_cmd(cmd, logger, check=True):
    logger.info("RUN: " + " ".join(str(c) for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stdout.strip():
        logger.info("STDOUT:\n" + proc.stdout.strip())
    if proc.stderr.strip():
        logger.info("STDERR:\n" + proc.stderr.strip())
    if check and proc.returncode != 0:
        logger.error(f"command exited {proc.returncode}: {' '.join(str(c) for c in cmd)}")
        sys.exit(proc.returncode)
    return proc


def run_cmd_stream_to_file(cmd, out_path, logger, check=True):
    """Like run_cmd, but stdout is streamed directly to out_path instead of being
    captured into memory and dumped whole into the log. Only stderr (expected to
    be small -- progress/warnings, not the actual data payload) is logged. Use
    this for any command whose stdout IS the output data (bedtools intersect,
    bedtools sort, etc) rather than a status message -- capturing multi-GB
    stdout into a Python string and then into a log file is both a memory risk
    and produces an unusable log."""
    logger.info(f"RUN (stdout -> {out_path}): " + " ".join(str(c) for c in cmd))
    with open(out_path, "w") as out_fh:
        proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if proc.stderr.strip():
        logger.info("STDERR:\n" + proc.stderr.strip())
    if check and proc.returncode != 0:
        logger.error(f"command exited {proc.returncode}: {' '.join(str(c) for c in cmd)}")
        sys.exit(proc.returncode)
    return proc


def run_cmd_discard_stdout(cmd, logger, check=True):
    """For validation-only commands whose stdout is neither the result we want
    nor something worth logging (e.g. `bedtools sort -i X` used only to confirm
    X parses) -- stdout goes to /dev/null, only stderr/returncode matter."""
    logger.info("RUN (stdout discarded, validation only): " + " ".join(str(c) for c in cmd))
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if proc.stderr.strip():
        logger.info("STDERR:\n" + proc.stderr.strip())
    if check and proc.returncode != 0:
        logger.error(f"command exited {proc.returncode}: {' '.join(str(c) for c in cmd)}")
        sys.exit(proc.returncode)
    return proc


def tool_version(cmd, logger):
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        out = (proc.stdout + proc.stderr).strip().splitlines()
        return out[0] if out else "(no output)"
    except FileNotFoundError:
        return None


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_eigenstrat_ind(path):
    """EIGENSTRAT .ind: whitespace-separated SampleID Sex Label. Returns list of dicts."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"{path}:{lineno}: expected 3 columns (ID Sex Label), got {len(parts)}: {line!r}")
            sample_id, sex, label = parts[0], parts[1], parts[2]
            rows.append({"id": sample_id, "sex": sex, "label": label})
    return rows


def read_eigenstrat_snp_header_probe(path, n=5):
    """Sniff column count of an EIGENSTRAT .snp file (4-col legacy vs 6-col with REF/ALT)."""
    with open(path) as fh:
        seen = []
        for line in fh:
            parts = line.split()
            seen.append(len(parts))
            if len(seen) >= n:
                break
    return seen


def iter_eigenstrat_snp(path):
    """Yield dicts per SNP line. Handles both 4-col (id chr genpos pos) and
    6-col (id chr genpos pos ref alt) EIGENSTRAT .snp formats."""
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            if len(parts) == 6:
                snpid, chrom, genpos, pos, ref, alt = parts
            elif len(parts) == 4:
                snpid, chrom, genpos, pos = parts
                ref = alt = None
            else:
                raise ValueError(f"{path}:{lineno}: unexpected .snp column count {len(parts)}: {line!r}")
            yield {"snpid": snpid, "chrom": chrom, "genpos": genpos, "pos": int(pos),
                   "ref": ref, "alt": alt}


def find_duplicate_ids(ids):
    """Return {id: count} for every id appearing more than once, preserving
    nothing about order -- used to hard-fail on duplicate sample/SNP IDs
    rather than silently keeping the first or last occurrence."""
    from collections import Counter
    counts = Counter(ids)
    return {k: v for k, v in counts.items() if v > 1}


def genomic_window_index(pos_1based, window_bp):
    """1-based genomic coordinate -> 0-based non-overlapping window index of
    size window_bp. Position 1..window_bp -> window 0, window_bp+1..2*window_bp
    -> window 1, etc. Using (pos-1)//window_bp (NOT pos//window_bp) is required
    for 1-based coordinates: pos==window_bp must fall in window 0, not window 1."""
    return (pos_1based - 1) // window_bp


def validate_panel_a_bim(bim_path):
    """Panel A primary-parameter structural checks (spec section 5): biallelic
    only, chromosomes 1-12 only, unique SNP ID, unique (chrom,pos) position.
    Returns a list of human-readable problem strings; empty list means clean.
    Does not modify anything -- callers must sys.exit on a nonempty result,
    not attempt to auto-fix."""
    problems = []
    valid_chroms = {str(c) for c in range(1, 13)}
    ids = []
    seen_pos = {}
    bad_chrom = set()
    non_biallelic = []
    dup_pos_examples = []
    with open(bim_path) as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            if len(parts) < 6:
                problems.append(f"{bim_path}:{lineno}: expected >=6 columns, got {len(parts)}")
                continue
            chrom, snpid, _cm, pos, a1, a2 = parts[:6]
            ids.append(snpid)
            if chrom not in valid_chroms:
                bad_chrom.add(chrom)
            a1u, a2u = a1.upper(), a2.upper()
            if (len(a1u) != 1 or len(a2u) != 1 or a1u == a2u
                    or a1u not in "ACGT" or a2u not in "ACGT"):
                non_biallelic.append((lineno, snpid, a1, a2))
            key = (chrom, pos)
            if key in seen_pos:
                dup_pos_examples.append((key, seen_pos[key], snpid))
            else:
                seen_pos[key] = snpid

    if bad_chrom:
        problems.append(f"chromosomes outside 1-12 found: {sorted(bad_chrom)}")
    if non_biallelic:
        sample = non_biallelic[:10]
        problems.append(f"{len(non_biallelic)} non-biallelic/non-ACGT allele rows, "
                         f"e.g. {sample}")
    dup_ids = find_duplicate_ids(ids)
    if dup_ids:
        sample = list(dup_ids.items())[:10]
        problems.append(f"{len(dup_ids)} duplicate SNP IDs, e.g. {sample}")
    if dup_pos_examples:
        sample = dup_pos_examples[:10]
        problems.append(f"{len(dup_pos_examples)} duplicate (chrom,pos) positions, "
                         f"e.g. {sample}")
    return problems


def write_manifest_tsv(path, rows, fieldnames):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


TRANSITION_PAIRS = {frozenset(("A", "G")), frozenset(("C", "T"))}
TRANSVERSION_PAIRS = {frozenset(("A", "C")), frozenset(("A", "T")),
                       frozenset(("C", "G")), frozenset(("G", "T"))}


PANEL_KEY_BY_LETTER = {"A": "panel_A_3k", "B": "panel_B_720", "C": "panel_C_civan"}


def resolve_marker_params(cfg, panel_letter, sensitivity):
    """Single authoritative lookup from (panel, sensitivity label) -> geno/maf/LD
    params, sourced only from config. No script may compute or override these
    numbers itself. Raises SystemExit on an unknown/invalid combination rather
    than silently defaulting."""
    panel_key = PANEL_KEY_BY_LETTER.get(panel_letter)
    if panel_key is None:
        raise SystemExit(f"unknown panel letter {panel_letter!r}")
    base = cfg[panel_key]
    if panel_key == "panel_B_720":
        window_kb, r2 = base["ld_primary"]["window_kb"], base["ld_primary"]["r2"]
    else:
        window_kb, r2 = base["ld_window_kb"], base["ld_r2"]
    maf = base["maf"]
    geno = base["geno"]

    if sensitivity == "primary":
        pass
    elif sensitivity in ("S1", "S2", "S3"):
        match = [s for s in cfg["sensitivity"] if s["label"] == sensitivity]
        if not match:
            raise SystemExit(f"sensitivity {sensitivity} not found in config.sensitivity")
        window_kb, r2 = match[0]["window_kb"], match[0]["r2"]
    elif sensitivity == "S4":
        if panel_key != "panel_A_3k":
            raise SystemExit("S4 (MAF=0.05 sensitivity) is defined only for panel_A_3k")
        maf = base["sensitivity_S4_maf"]
    else:
        raise SystemExit(f"unknown sensitivity label {sensitivity!r}")

    return {"panel_key": panel_key, "geno": geno, "maf": maf,
            "ld_window_kb": window_kb, "ld_r2": r2}


def is_transversion(a1, a2):
    a1, a2 = a1.upper(), a2.upper()
    if len(a1) != 1 or len(a2) != 1 or a1 == a2:
        return None
    pair = frozenset((a1, a2))
    if pair in TRANSVERSION_PAIRS:
        return True
    if pair in TRANSITION_PAIRS:
        return False
    return None
