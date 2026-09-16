import os, glob, re, sys
cancer = sys.argv[1]
raw = f"/public/home/fjhui/ZW/{cancer}/WSI/raw"
out = f"/public/home/fjhui/ZW/scripts/slide_type_map_{cancer}.tsv"
rows = set()
for p in glob.glob(raw + "/**/*.svs", recursive=True):
    sid = os.path.basename(p)[:-4]
    case = re.sub(r"-\d+$", "", sid)
    rows.add((sid, "Primary Tumor", case))
with open(out, "w") as f:
    f.write("slide_submitter_id\tsample_type\tcase_submitter_id\n")
    for sid, st, case in sorted(rows): f.write(f"{sid}\t{st}\t{case}\n")
print(f"{cancer}: wrote {len(rows)} slides (all Primary Tumor) -> {out}")
