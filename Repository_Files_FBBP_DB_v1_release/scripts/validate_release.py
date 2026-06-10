from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "tables/proteins.csv": 1996,
    "tables/loop_annotations.csv": 3383,
    "tables/loop_flexibility_public_summary.csv": 3383,
    "tables/protein_flexibility_summary.csv": 1996,
}

def count_rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        return sum(1 for _ in reader)

failed = False
for rel, expected in EXPECTED.items():
    path = ROOT / rel
    observed = count_rows(path)
    status = "OK" if observed == expected else "FAIL"
    print(f"{status}\t{rel}\tobserved={observed}\texpected={expected}")
    failed = failed or observed != expected

raise SystemExit(1 if failed else 0)
