"""Build SuppTable_provenance.tex: one row per registered rule, pre-specification, amendment and
results record cited in the Methods (tab:provenance).

Nothing in the printed hash and time columns is typed by hand:
  * SHA-256: the first hash on the first non-comment line of the record's sidecar
    review/<document>.sha256 (the hash of the document itself, as printed there);
  * "Sidecar written (UTC)": the "# written <ISO time>" line of that sidecar, where present;
  * tag creation time: read from the local git repository (git for-each-ref, read only) and
    converted to UTC; the script fails if a listed tag does not exist.
The record list below (document, tag, what it governs) is the only hand-written part; every path
in it is checked to exist.

Usage:  python make_supptable_provenance.py      (writes SuppTable_provenance.tex next to it)
"""
import datetime as dt
import hashlib
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "SuppTable_provenance.tex")

# kind, document (relative to HERE), sidecar (relative to HERE) or None, public tag or None,
# what the record governs (LaTeX).  Results records carry their own per-file sums file, given
# as the document, and have no review/*.sha256 sidecar.
RECORDS = [
    ("Dated rule (not registered or publicly released)", "review/C1_UCEC_FROZEN_RULE_2026-09-15.md", None, None,
     r"UCEC set-level confirmatory test, with its two amendments (Methods~\S\ref{sec:c1ucecmeth})"),
    (r"Registered rule (signed; OSF \texttt{2d6mb})", "review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md",
     "review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md.sha256", "c1-luad-signed-rule-2026-10-01",
     r"LUAD set-level confirmatory test (Methods~\S\ref{sec:c1luadmeth})"),
    ("Post-signature addenda", "review/C1_LUAD_ADDENDA.md", "review/C1_LUAD_ADDENDA.md.sha256",
     "c1-luad-addenda-2026-10-02",
     r"LUAD run log, unblinding bound and the out-of-order check (Supplementary Note~\ref{snote:c1luad})"),
    ("Addenda, Entry~11", "review/C1_LUAD_ADDENDA.md", "review/C1_LUAD_ADDENDA.md.sha256",
     "c1-luad-addenda-entry11-2026-10-10", r"Rewritten Abstract and its fixed wording"),
    ("Addenda, Entry~12", "review/C1_LUAD_ADDENDA.md", "review/C1_LUAD_ADDENDA.md.sha256",
     "c1-luad-addenda-entry12-2026-11-20", r"Second rewrite of the Abstract; where the fixed wording stands"),
    ("Pre-specification", "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md",
     "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md.sha256", "c1-luad-posthoc-prespec-2026-10-02",
     r"Post hoc plex and medium axes, FFPE-only and Disqualified arms (Methods~\S\ref{sec:batchmeth})"),
    ("Amendment~1", "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md",
     "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md.sha256", "c1-luad-posthoc-amend1-2026-10-03",
     r"Arm builder's analysed population; Disqualified-aliquot findings"),
    ("Amendment~2", "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md",
     "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md.sha256", "c1-luad-posthoc-amend2-2026-10-03",
     r"Python-version defect in the disclosure writer; a stale comment"),
    ("Amendment~2, correction", "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md",
     "review/POSTHOC_BATCH_AXES_PRESPEC_2026-10-02.md.sha256", "c1-luad-posthoc-amend2-corr-2026-10-03",
     r"One-sentence factual correction to Amendment~2"),
    ("Results record", "server_export/pinned/posthoc", None, "c1-luad-posthoc-results-2026-10-04",
     r"Batch-axis results (addenda entries 9--10; Supplementary Table~\ref{tab:posthoc})"),
    ("Pre-specification", "review/POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md",
     "review/POSTHOC_UCEC_C1_SETS_PRESPEC_2026-10-08.md.sha256", "c1-ucec-posthoc-sets-prespec-2026-10-08",
     r"UCEC confirmatory re-read with alternative selection sets (Supplementary Note~\ref{snote:c1posthoc})"),
    ("Results record", "server_export/results/c1_ucec_posthoc_sets", None,
     "c1-ucec-posthoc-sets-results-2026-10-08", r"Supplementary Tables~\ref{tab:c1ucecsets} and~\ref{tab:c1ucecsetsref}"),
    ("Pre-specification", "review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md",
     "review/POSTHOC_UCEC_C1_PLEX_PRESPEC_2026-10-08.md.sha256", "c1-ucec-posthoc-plex-prespec-2026-10-08",
     r"UCEC confirmatory re-read within TMT plex (Supplementary Note~\ref{snote:c1posthoc})"),
    ("Results record", "server_export/results/c1_ucec_plex/SHA256SUMS.txt", None,
     "c1-ucec-posthoc-plex-results-2026-10-09", r"Supplementary Table~\ref{tab:c1ucecplex}"),
    ("Pre-specification", "review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md",
     "review/POSTHOC_MEAS_ERROR_PRESPEC_2026-10-08.md.sha256", "meas-error-posthoc-prespec-2026-10-08",
     r"Own-mRNA measurement-error test (Methods~\S\ref{sec:txomemeth})"),
    ("Amendment~1", "review/POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md",
     "review/POSTHOC_MEAS_ERROR_AMENDMENT1_2026-10-10.md.sha256", "meas-error-posthoc-amendment1-2026-10-10",
     r"UCEC rerun and the LUAD, GBM and PDAC runs (Supplementary Note~\ref{snote:measerr})"),
    ("Results record", "review/MEAS_ERROR_2026-10-08/RESULTS_SHA256SUMS.txt", None,
     "meas-error-posthoc-results-2026-10-10", r"Supplementary Table~\ref{tab:measerr}"),
    ("Pre-specification", "review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md",
     "review/POSTHOC_BRCA_TS_PRESPEC_2026-10-08.md.sha256", "brca-ts-posthoc-prespec-2026-10-08",
     r"TCGA-BRCA frozen-section arm (Methods~\S\ref{sec:kircmeth})"),
    ("Amendment~1", "review/POSTHOC_BRCA_TS_AMENDMENT1_2026-10-08.md",
     "review/POSTHOC_BRCA_TS_AMENDMENT1_2026-10-08.md.sha256", "brca-ts-posthoc-amendment1-2026-10-08",
     r"Diagnostic-slide embedding folder of the re-run (Supplementary Note~\ref{snote:brcats})"),
    ("Results record", "review/BRCA_TS_ARM_2026-10-08/RESULTS_SHA256SUMS.txt", None,
     "brca-ts-posthoc-results-2026-10-10", r"Supplementary Table~\ref{tab:brcats}"),
    ("Pre-specification", "review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md",
     "review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md.sha256", "c1-covariates-posthoc-prespec-2026-11-21",
     r"Confirmatory re-reads with transcriptome-derived covariates (Methods~\S\ref{sec:c1luadmeth}; Supplementary Note~\ref{snote:c1cov})"),
    ("Stage~H record", "review/C1_COVARIATES_2026-11-21/stage_h_record/c1_cov_stage_h_manifest.sha256", None,
     "c1-covariates-posthoc-stageh-2026-11-22", r"Manifest of the confirmatory arms (310 files; file hashes and logs only)"),
    ("Results record", "review/C1_COVARIATES_2026-11-21/results_record/RESULTS_SHA256.txt", None,
     "c1-covariates-posthoc-results-2026-11-23", r"Supplementary Tables~\ref{tab:c1cov} and~\ref{tab:c1covsets}"),
    ("Addenda, Entry~13", "review/C1_LUAD_ADDENDA.md", "review/C1_LUAD_ADDENDA.md.sha256",
     "c1-luad-addenda-entry13-2026-11-23", r"Abstract after the covariate re-reads (Case 1; which reading decided it)"),
]

CAPTION = (
    r"\textbf{Registered rules, pre-specifications, amendments and results records.} One row per "
    r"record cited in the Methods. SHA-256, the hash of the document as printed in its sidecar "
    r"(\texttt{review/<document>.sha256}) or, for an earlier version of a file that later changed, the hash "
    r"of the file as committed at its tag (read with \texttt{git show}); written, the sidecar's own UTC time stamp (the three "
    r"oldest sidecars carry a bare hash and no time). Where one file carries several tags (the "
    r"addenda; the batch-axis specification and its amendments), only the current version has a "
    r"sidecar; rows marked earlier version give the hash of the file at that tag and no written time. Tag "
    r"times are the time stamps git records for each public tag (the commit time for a lightweight "
    r"tag, the tagger time for an annotated one), in UTC, not GitHub server times; tag names carry the local "
    r"(UTC+8) date, so a tag can be named one day later than its UTC time. Results records have no "
    r"sidecar; their per-file hashes are in the sums file named, or, for the two results folders "
    r"without one, in the tagged commit. The UCEC rule has no sidecar, tag or third-party "
    r"timestamp. The signed LUAD rule keeps its pre-signature file name and header because it "
    r"could not be edited after hashing; a receipt in the repository links it to the OSF "
    r"registration and the release; every script and input was hash-asserted at run time, and "
    r"smoke runs printed counts only. Generated by \texttt{make\_supptable\_provenance.py} from "
    r"the sidecars and the repository's tags."
)


def tt(s):
    """Typewriter text with break points after / _ - ."""
    t = s.replace("\\", "/").replace("_", r"\_")
    for ch in ("/", r"\_", "-", "."):
        t = t.replace(ch, ch + r"\allowbreak{}")
    return r"\texttt{" + t + "}"


def read_sidecar(path):
    h, written = None, None
    with open(os.path.join(HERE, path), encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if s.startswith("# written"):
                    written = s.split("written", 1)[1].strip()
                continue
            if h is None:
                h = s.split()[0]
    assert h is not None and len(h) == 64 and all(c in "0123456789abcdef" for c in h), (path, h)
    return h, written


def tag_time_utc(tag):
    out = subprocess.run(["git", "-C", HERE, "for-each-ref", "--format=%(creatordate:iso-strict)",
                          "refs/tags/" + tag], capture_output=True, text=True, check=True).stdout.strip()
    assert out, "tag not found: " + tag
    t = dt.datetime.fromisoformat(out).astimezone(dt.timezone.utc)
    return t.strftime("%Y-%m-%d %H:%M:%S")


def sha_at_tag(tag, doc):
    """SHA-256 of the document's bytes as committed at the tag (git show, read only)."""
    blob = subprocess.run(["git", "-C", HERE, "show", tag + ":" + doc],
                          capture_output=True, check=True).stdout
    return hashlib.sha256(blob).hexdigest()


def fmt_written(w):
    if not w:
        return "---"
    t = dt.datetime.fromisoformat(w.replace("Z", "+00:00")).astimezone(dt.timezone.utc)
    return t.strftime("%Y-%m-%d") + r"\newline " + t.strftime("%H:%M:%S")


def main():
    rows = []
    for kind, doc, side, tag, gov in RECORDS:
        assert os.path.exists(os.path.join(HERE, doc)), doc
        if side:
            h, w = read_sidecar(side)
            if tag:
                ht = sha_at_tag(tag, doc)
                if ht != h:  # an earlier version of a file that later changed
                    h, w = ht, None
                    kind = kind + r" (earlier version)"
            hcell = r"\texttt{" + h[:32] + r"}\newline\texttt{" + h[32:] + "}"
        else:
            h, w = None, None
            hcell = "no sidecar"
        if tag:
            tcell = tt(tag) + r"\newline " + tag_time_utc(tag)
        else:
            tcell = "none"
        cells = [kind, tt(doc), hcell, fmt_written(w), tcell, gov]
        rows.append(" & ".join(r"{\raggedright " + c + r"\par}" for c in cells) + r" \\")

    head = r"Record & Document & SHA-256 & Written\newline (UTC) & Public tag, time (UTC) & Governs \\"
    out = [
        "% Supplementary Table: provenance of registered rules, pre-specifications, amendments and",
        "% results records (tab:provenance). Generated by make_supptable_provenance.py; do not edit by hand.",
        r"\begin{landscape}",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\renewcommand{\arraystretch}{1.15}",
        r"\setlength{\LTcapwidth}{\linewidth}",
        r"\begin{longtable}{p{2.5cm}p{4.0cm}p{4.9cm}p{1.5cm}p{4.0cm}p{4.1cm}}",
        r"\captionsetup{font={scriptsize,stretch=1.0}}",
        r"\caption{" + CAPTION + r"}\label{tab:provenance}\\",
        r"\toprule",
        head,
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{6}{l}{\emph{Supplementary Table~\ref{tab:provenance}, continued}}\\",
        r"\toprule",
        head,
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endlastfoot",
    ] + rows + [
        r"\end{longtable}",
        r"\endgroup",
        r"\end{landscape}",
    ]
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    print("wrote", OUT, "rows", len(rows))


if __name__ == "__main__":
    main()
