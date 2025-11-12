import json, sys, pathlib as p
j = json.load(open("raw/encsr.json")) if p.Path("raw/encsr.json").exists() else json.load(open("encsr.json"))
fx = [f for f in j.get("files", []) if f.get("file_format")=="fastq" and f.get("status")=="released"]
pairs = {}
for f in fx:
    pe = f.get("paired_end")
    if pe not in ("1","2"): continue
    rep = (f.get("biological_replicates") or ["rep?"])[0]
    size = int(f.get("file_size", 0) or 0)
    acc  = f["accession"]
    url  = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.fastq.gz"
    pairs.setdefault(rep, {}).setdefault(pe, []).append((size, acc, url))
for rep, d in pairs.items():
    if "1" in d and "2" in d and d["1"] and d["2"]:
        r1 = max(d["1"])[2]; r2 = max(d["2"])[2]
        print(r1); print(r2)
        sys.exit(0)
sys.exit("No paired-end FASTQs found")
