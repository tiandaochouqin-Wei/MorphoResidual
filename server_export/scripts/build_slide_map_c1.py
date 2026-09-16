#!/usr/bin/env python3
"""
GDC slide-entity sample_type map for a given list of CPTAC slide ids, written in the
same format and via the same lookup as build_slide_map_gdc.py (fetch_slide_entities is
imported, not re-implemented). Needed for C1-UCEC because the map has to exist before
extraction, so ids come from the DICOM ContainerIdentifier rather than .pt filenames.
Slides without a GDC slide entity get no row, exactly as in the discovery maps.

  python build_slide_map_c1.py ids.txt out.tsv        # ids.txt: one slide id per line
"""
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_slide_map_gdc as B  # noqa: E402


def main():
    ids_path, out = sys.argv[1], sys.argv[2]
    ids = sorted({line.strip() for line in open(ids_path, encoding="utf-8") if line.strip()})
    cases = set()
    for sid in ids:
        m = B.SLIDE_RE.match(sid)
        if not m:
            sys.exit(f"not a CPTAC slide id: {sid!r}")
        cases.add(m.group(1))
    entities, _ = B.fetch_slide_entities(cases)
    mapped = {sid: entities[sid] for sid in ids if sid in entities}
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("slide_submitter_id\tsample_type\tcase_submitter_id\tsection_location\n")
        for sid in sorted(mapped):
            stype, case, loc = mapped[sid]
            f.write(f"{sid}\t{stype}\t{case}\t{loc if loc else ''}\n")
    tally = Counter(v[0] for v in mapped.values())
    print(f"{len(ids)} slide ids / {len(cases)} cases: {len(mapped)} with a GDC slide entity, "
          f"{len(ids) - len(mapped)} unmatched; {dict(tally)}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
