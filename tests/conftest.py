"""Small synthetic rMATS fixtures; never read patient data."""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

EVENT_IDS = [
    "A5SS:CTDSP2:chr12:57827533:57827590:57827551:57827590:57826995:57827097",
    "A5SS:MAP4K4:chr2:101887087:101887261:101887087:101887237:101887777:101887937",
]
HEADERS = [
    "ID",
    "GeneID",
    "geneSymbol",
    "chr",
    "strand",
    "longExonStart_0base",
    "longExonEnd",
    "shortES",
    "shortEE",
    "flankingES",
    "flankingEE",
    "ID",
    "IJC_SAMPLE_1",
    "SJC_SAMPLE_1",
    "IJC_SAMPLE_2",
    "SJC_SAMPLE_2",
    "FDR",
]


@pytest.fixture
def model_file(tmp_path):
    path = tmp_path / "model.toml"
    path.write_text("""schema_version = 1
name = "synthetic"
reference_genome = "GRCh38"
intercept = -2.927847423
[[events]]
id = "A5SS:CTDSP2:chr12:57827533:57827590:57827551:57827590:57826995:57827097"
strand = "-"
direction = 1
weight = 2.0
[[events]]
id = "A5SS:MAP4K4:chr2:101887087:101887261:101887087:101887237:101887777:101887937"
strand = "+"
direction = -1
weight = 0.5
""")
    return path


@pytest.fixture
def rmats_rows():
    rows = []
    for i, (event, strand, counts) in enumerate(
        zip(EVENT_IDS, ["-", "+"], [("4", "0", "0,4", "4,4"), ("0", "4", "4,4", "0,4")])
    ):
        _, gene, chrom, *coords = event.split(":")
        rows.append(
            [
                str(i),
                "ENSG_SYNTHETIC",
                gene,
                chrom,
                strand,
                *coords,
                str(i),
                *counts,
                "1.0",
            ]
        )
    return rows


@pytest.fixture
def write_rmats(tmp_path, rmats_rows):
    def write(rows=None, headers=None, bom=False):
        path = tmp_path / "A5SS.MATS.JC.txt"
        with path.open("w", encoding="utf-8-sig" if bom else "utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(HEADERS if headers is None else headers)
            writer.writerows(rmats_rows if rows is None else rows)
        return path

    return write
