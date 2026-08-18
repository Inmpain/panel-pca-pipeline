#!/usr/bin/env python3
"""Shared, dependency-light helpers for ecotype PCA v2 scripts 09--18."""

from __future__ import annotations

import csv
import hashlib
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


VALID_GENOTYPES = frozenset("0129")
POOLED_LIBRARY_TYPE = "pooled_mixed"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_int_seed(*parts: object) -> int:
    payload = ":".join(str(part) for part in parts).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def atomic_text_writer(final_path: str | Path):
    """Context manager-like writer with explicit commit/abort methods."""
    final = Path(final_path)
    final.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=final.name + ".", suffix=".tmp", dir=final.parent)
    handle = os.fdopen(fd, "w")

    class Writer:
        def write(self, value: str) -> None:
            handle.write(value)

        def commit(self) -> None:
            handle.close()
            os.replace(temporary, final)

        def abort(self) -> None:
            if not handle.closed:
                handle.close()
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass

    return Writer()


def refuse_existing(paths: list[str | Path], overwrite: bool = False) -> None:
    existing = [str(path) for path in paths if Path(path).exists()]
    if existing and not overwrite:
        raise ValueError("refusing to overwrite existing output(s): " + ", ".join(existing))


def read_ind(path: str | Path) -> list[dict[str, str]]:
    rows = []
    seen = set()
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected ID SEX LABEL")
            sample_id, sex, label = fields[:3]
            if sample_id in seen:
                raise ValueError(f"{path}:{line_no}: duplicate sample ID {sample_id!r}")
            seen.add(sample_id)
            rows.append({"id": sample_id, "sex": sex, "label": label})
    if not rows:
        raise ValueError(f"{path}: empty .ind file")
    return rows


def read_keep_ids(path: str | Path) -> list[str]:
    result = []
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            result.append(fields[1] if len(fields) >= 2 else fields[0])
    duplicates = [key for key, n in Counter(result).items() if n > 1]
    if duplicates:
        raise ValueError(f"{path}: duplicate keep IDs: {duplicates[:10]}")
    if not result:
        raise ValueError(f"{path}: empty keep file")
    return result


def iter_snp(path: str | Path):
    seen_ids = set()
    seen_positions = set()
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) not in (4, 6):
                raise ValueError(f"{path}:{line_no}: expected 4 or 6 .snp columns")
            snp_id, chrom, genetic_pos, pos_text = fields[:4]
            try:
                pos = int(pos_text)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid position {pos_text!r}") from exc
            if snp_id in seen_ids:
                raise ValueError(f"{path}:{line_no}: duplicate SNP ID {snp_id!r}")
            if (chrom, pos) in seen_positions:
                raise ValueError(f"{path}:{line_no}: duplicate coordinate {(chrom, pos)!r}")
            seen_ids.add(snp_id)
            seen_positions.add((chrom, pos))
            ref = alt = None
            if len(fields) == 6:
                ref, alt = fields[4].upper(), fields[5].upper()
                if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                    raise ValueError(f"{path}:{line_no}: invalid REF/ALT {ref}/{alt}")
            yield {
                "id": snp_id, "chrom": chrom, "genetic_pos": genetic_pos,
                "pos": pos, "ref": ref, "alt": alt, "line": line.rstrip("\n"),
            }


def iter_panel_snp(path: str | Path):
    """Lenient reader for a raw upstream panel .snp file (e.g. civan_snp.snp,
    2.36M rows): silently skips any line that is not a clean 6-column (id
    chrom genpos pos ref alt) row, matching the tolerance the Stage 50
    prototype's own marker selection already used against this same file.
    Unlike iter_snp (strict; meant for already-curated fixed-marker subsets),
    this is for the raw upstream panel, which is not guaranteed clean."""
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                continue
            snp_id, chrom, _genpos, pos_text, ref, alt = fields
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            yield {"id": snp_id, "chrom": chrom, "pos": pos, "ref": ref.upper(), "alt": alt.upper()}


def read_calls(path: str | Path) -> str:
    calls = []
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            if len(value) != 1 or value not in VALID_GENOTYPES:
                raise ValueError(f"{path}:{line_no}: expected one of 0/1/2/9, got {value!r}")
            calls.append(value)
    if not calls:
        raise ValueError(f"{path}: empty call file")
    return "".join(calls)


def parse_sample_paths(values: list[str]) -> list[tuple[str, Path]]:
    result = []
    seen = set()
    for value in values:
        if "=" not in value:
            raise ValueError(f"expected SAMPLE=PATH, got {value!r}")
        sample, path_text = value.split("=", 1)
        if not sample or sample in seen:
            raise ValueError(f"empty or duplicate sample ID in {value!r}")
        path = Path(path_text)
        if not path.is_file():
            raise ValueError(f"call file not found: {path}")
        seen.add(sample)
        result.append((sample, path))
    if not result:
        raise ValueError("at least one SAMPLE=PATH entry is required")
    return result


def write_tsv(path: str | Path, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def broad_class(panel: str, label: str) -> str | None:
    if panel == "C":
        mapping = {
            "indica": "INDICA", "aus": "AUS", "aromatic": "AROMATIC",
            "japonica": "JAPONICA", "japonica_(temperate)": "JAPONICA",
            "japonica_(tropical)": "JAPONICA",
        }
    elif panel == "A":
        mapping = {key: key for key in ("IND", "AUS", "ARO", "TRJ", "TEJ")}
    else:
        return label
    return mapping.get(label)


def read_evec(path: str | Path, num_pcs: int = 10) -> list[dict]:
    rows = []
    with open(path) as handle:
        header = handle.readline()
        if not header.lstrip().startswith("#"):
            raise ValueError(f"{path}: smartpca .evec header must start with #")
        for line_no, line in enumerate(handle, 2):
            fields = line.split()
            if len(fields) < num_pcs + 2:
                raise ValueError(f"{path}:{line_no}: fewer than {num_pcs} PCs")
            try:
                pcs = [float(value) for value in fields[1:num_pcs + 1]]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid PC coordinate") from exc
            rows.append({"id": fields[0], "pcs": pcs, "label": fields[-1]})
    if not rows:
        raise ValueError(f"{path}: no sample rows")
    return rows


def euclidean(left: list[float], right: list[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def class_centroids(rows: list[dict], panel: str, dimensions: int) -> dict[str, list[float]]:
    grouped = defaultdict(list)
    for row in rows:
        group = broad_class(panel, row["label"])
        if group is not None:
            grouped[group].append(row["pcs"][:dimensions])
    return {
        group: [sum(vector[i] for vector in vectors) / len(vectors) for i in range(dimensions)]
        for group, vectors in grouped.items()
    }


def nearest_group(rows: list[dict], panel: str, query_pcs: list[float], dimensions: int, k: int) -> tuple[str, int, str]:
    candidates = []
    for row in rows:
        group = broad_class(panel, row["label"])
        if group is not None:
            candidates.append((euclidean(query_pcs[:dimensions], row["pcs"][:dimensions]), row["id"], group))
    if len(candidates) < k:
        raise ValueError(f"nearest-{k} requires at least {k} eligible modern references; found {len(candidates)}")
    nearest = sorted(candidates)[:k]
    counts = Counter(group for _, _, group in nearest)
    max_count = max(counts.values())
    winners = sorted(group for group, count in counts.items() if count == max_count)
    winner = winners[0] if len(winners) == 1 else "UNRESOLVED_TIE"
    detail = ";".join(f"{group}:{counts[group]}" for group in sorted(counts))
    return winner, max_count, detail


def format_contig(contig_format: str, chrom: str) -> str:
    """EIGENSTRAT .snp chrom column -> BAM-style contig name, e.g. '1' or '01'
    -> 'chr01' with the default contig_format. Shared by every script that
    must match panel chromosome labels against BAM reference names."""
    try:
        return contig_format % int(chrom.lstrip("0") or "0")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unsupported chromosome {chrom!r}") from exc


def tally_coverage(sample_covered: dict[str, set]) -> dict[object, set[str]]:
    """{sample: {site, ...}} -> {site: {sample, ...}}. Pure aggregation, no I/O."""
    tally: dict[object, set[str]] = defaultdict(set)
    for sample, sites in sample_covered.items():
        for site in sites:
            tally[site].add(sample)
    return dict(tally)
