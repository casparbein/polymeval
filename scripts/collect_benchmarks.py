import csv
import json

CELLS      = snakemake.params.cells
BENCHMARKS = snakemake.params.benchmarks

AARDVARK_TYPE = {"Snv": "SNP", "Indel": "INDEL"}
FIELDS = ["sample", "engine", "caller", "benchmark", "build", "cls", "stratum", "filter",
          "truth_total", "tp", "fn", "query_total", "fp", "recall", "precision", "f1"]

def num(x):
    """hap.py leaves cells empty and writes literal nan; keep those as blanks."""
    x = str(x or "").strip()
    if x.lower() in ("", "nan", "na", "none", "."):
        return ""
    try:
        v = float(x)
    except ValueError:
        return ""
    ## keep counts as integers so the table is readable
    return int(v) if v.is_integer() and "." not in x and "e" not in x.lower() else v


def read_happy(path):
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            yield r["Type"], r["Filter"], dict(
                truth_total=num(r["TRUTH.TOTAL"]), tp=num(r["TRUTH.TP"]), fn=num(r["TRUTH.FN"]),
                query_total=num(r["QUERY.TOTAL"]), fp=num(r["QUERY.FP"]),
                recall=num(r["METRIC.Recall"]), precision=num(r["METRIC.Precision"]),
                f1=num(r["METRIC.F1_Score"]))


def read_truvari(path, cls):
    with open(path) as fh:
        d = json.load(fh)
    yield cls.upper(), "NA", dict(
        truth_total=num(d.get("base cnt")), tp=num(d.get("TP-base")), fn=num(d.get("FN")),
        query_total=num(d.get("comp cnt")), fp=num(d.get("FP")),
        recall=num(d.get("recall")), precision=num(d.get("precision")), f1=num(d.get("f1")))


def read_aardvark(path):
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ## BASEPAIR rows count bases, not variants, and would not be comparable
            ## with the other two engines. Stratification regions are kept separately.
            if r["comparison"] != "GT" or r["region_label"] != "ALL":
                continue
            yield AARDVARK_TYPE.get(r["variant_type"], r["variant_type"].upper()), r["filter"], dict(
                truth_total=num(r["truth_total"]), tp=num(r["truth_tp"]), fn=num(r["truth_fn"]),
                query_total=num(r["query_total"]), fp=num(r["query_fp"]),
                recall=num(r["metric_recall"]), precision=num(r["metric_precision"]),
                f1=num(r["metric_f1"]))


rows = []
for engine, caller, truth, sample, path in CELLS:
    b = BENCHMARKS[truth]
    if engine == "happy":
        parsed = read_happy(path)
    elif engine in ("truvari", "truvari_refine"):
        parsed = read_truvari(path, b["cls"])
    else:
        parsed = read_aardvark(path)
    for stratum, filt, metrics in parsed:
        rows.append(dict(sample=sample, engine=engine, caller=caller, benchmark=truth,
                         build=b["build"], cls=b["cls"], stratum=stratum, filter=filt, **metrics))

rows.sort(key=lambda r: tuple(str(r[k]) for k in
                              ("sample", "caller", "benchmark", "stratum", "filter", "engine")))

with open(snakemake.output.long, "w", newline="") as fh:
    w = csv.DictWriter(fh, FIELDS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

## Wide view: one row per cell, one column per engine, F1 as the value.
engines = sorted({r["engine"] for r in rows})
index = ["sample", "caller", "benchmark", "stratum", "filter"]
grid = {}
for r in rows:
    grid.setdefault(tuple(r[k] for k in index), {})[r["engine"]] = r["f1"]

with open(snakemake.output.matrix, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(index + engines)
    for key in sorted(grid):
        w.writerow(list(key) + [grid[key].get(e, "") for e in engines])