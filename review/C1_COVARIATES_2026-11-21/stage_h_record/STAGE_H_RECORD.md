# Stage H record: C1 covariates post hoc re-reads (2026-11-22)

Record required by section 9 of `review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md` before any local stage.
It holds file hashes and logs only. **No statistic of any Stage H output has been read by anyone**: the
first reading is Stage L1, which runs only after this record is released.

## What was run

All 16 production jobs of the pre-specification were submitted once, together, at 2026-11-21 10:01 CST
(02:01 UTC) with `submit_stage_h.sh` (the commands of `SERVER_RUNBOOK.md` sections 4b-4e; it first checks
the sha256 of the four frozen files on the cluster). Sidecar sha256
`e5391d15d70781b80d0e9aa1257cfa04a9bb725500c2924ae60e474a2ea9cf07`. Pinned hosts: UCEC s001, LUAD s002
(`PINNED_HOSTS.txt`). Nothing is closed (`CLOSED_PARTS.txt` has no entry). Every run used numpy 1.26.4 (private
archive sha256 `92f2aeb5603a18ed35b2242db5fa07f6ac43a96104e48a2b0ffa76f3663ef8ba`).

| Cohort | Job (id) | Run time |
|---|---|---|
| UCEC | A_C3 (76713879) | 34.2 min |
| UCEC | A_R20 (76713880) | 40.7 min |
| UCEC | V_C3 (76713881) | 27.1 min |
| UCEC | V_R20 (76713882) | 31.9 min |
| UCEC | placebo draws 1-30, 31-60, 61-90, 91-119 (76713883 to 76713886) | 35.1, 35.1, 33.7, 31.0 min |
| LUAD | A_C3 (76713887) | 71.0 min |
| LUAD | A_R20 (76713888) | 83.5 min |
| LUAD | V_C3 (76713889) | 58.4 min |
| LUAD | V_R20 (76713890) | 69.3 min |
| LUAD | placebo draws 1-30, 31-60, 61-90, 91-119 (76713891 to 76713894) | 68.0, 66.7, 66.0, 64.9 min |

Every job completed at attempt 1 (`DONE:` in its log); no job was interrupted, refused, or ended in an
exception. Each cohort holds 4 arm matrices, 119 placebo matrices and one covariates record.

## The manifest

`c1_cov_stage_h_manifest.sha256` (this folder) lists 310 files of the two output folders (sha256 of each, bare
names). **Manifest sha256 `ca39784874776c7fe0e8faf3c315b1ff8fcb642a9ea06992b9f53813db8c3e13`**; this is the
value Stage L1 and Stage L2 are given on their command lines (section 9). It was written once, by job 76720641
(2026-11-21 18:29 UTC). A second, concurrent submission of the same command (job 76720640) was refused
(`a leftover ... .partial exists`) because the first job's temporary file existed; it wrote nothing. Both logs are
in `joblogs/`.

## Transfer and verification

The 310 files were packed on the cluster into `c1_cov_stage_h.tgz` (1,119,232,294 bytes, sha256
`b4b360543f0bd087411432eb8e7a4b0c734bef9e804037ea75885e8313fd1c20`) and the job logs into
`c1_cov_stage_h_joblogs.tgz` (14,983 bytes, sha256
`0b823ad08d6700149ccadbf7a16e5e3dc69627a34bf9ee8e33ccac49490ffdf8`), downloaded, unpacked, and each of the 310 files
was compared with its manifest entry by sha256: 310 identical, 0 missing, 0 different. The folder also holds the
two discovery score tables (not in the manifest; hashed in the sidecar) and 16 empty run-lock files.

## What this folder contains and what is withheld

Released here: the manifest, the 18 job logs (`joblogs/`), and the 48 small run files (attempt and run records
and logs; `run_files/`) that are not already in `gate_logs/`. They contain hashes, hosts, versions and run times; the
logs print no result statistic. Withheld, fixed by the manifest: the 248 matrix files (`.npz`), the 8
arm matrices and 238 placebo matrices and the 2 covariates records.

## Next

Stage L1 (`--stage L1 --manifest-sha ca39784874776c7fe0e8faf3c315b1ff8fcb642a9ea06992b9f53813db8c3e13`), then Stage L2
with the same value, then the results record.
