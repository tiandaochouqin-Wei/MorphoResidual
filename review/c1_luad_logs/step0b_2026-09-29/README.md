# C1-LUAD Step 0b logs (metadata reads and DICOM header range reads, 2026-09-28/29 UTC)

Copies of the logs written by the pre-signature Step 0b metadata checks of the signed
C1-LUAD rule (`review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md`, preamble and D3).
`query_log_idc.tsv` records every HTTP range read of the 135 confirmatory DICOM series.
The raw API responses they point to (`raw/...`) are not published (they contain
case-level biospecimen records).

One change was made to the copies: the Windows user name in the local idc-index cache
path was replaced by `<user>`. Nothing else differs. sha256 of each original and copy:

| File | original sha256 | published copy sha256 | prefixes replaced |
|---|---|---|---|
| `query_log_gdc.tsv` | `cf9d7676e3cac1528ef02d0bfeb2fa0420af2c0617e18b7b92bd0e1ccdaa3ce9` | `cf9d7676e3cac1528ef02d0bfeb2fa0420af2c0617e18b7b92bd0e1ccdaa3ce9` | 0 |
| `query_log_idc.tsv` | `e723b3b37200f4c1dfd902e17a396bc8ff648ffeb6c1ec72580fda35e3b8b767` | `4fd4f19a1886f0750b6a84980e08adca49cc398db0138d82b039245f0d87ab0f` | 2 |
| `query_log_pdc.tsv` | `fa1e66361e639e885e8e953eafa5141c2123736817cc451b2daf20cfeb92b18b` | `fa1e66361e639e885e8e953eafa5141c2123736817cc451b2daf20cfeb92b18b` | 0 |
| `query_log_protein.tsv` | `e538f2f96e4c935b52068c09412bba7c35b9c4675ee64e2215e245d59a15abac` | `e538f2f96e4c935b52068c09412bba7c35b9c4675ee64e2215e245d59a15abac` | 0 |
| `idc_step3_run.log` | `b61edde108abb02b216ef38dc34da4c41f956336d90952b618c217261198e385` | `85a87f7f38e6835ba35f5839a74d34b579d3285316c7ce5462c9bb7685b3b104` | 2 |
| `run_log_build_slide_map_c1.txt` | `c24d9319357cb9ef5da24ced5750fcc37c79be28698716deaeed1a428d083f40` | `c24d9319357cb9ef5da24ced5750fcc37c79be28698716deaeed1a428d083f40` | 0 |
| `run_log_gdc_step1.txt` | `6f549d3a61f2ba5c1a45e3c28f4c2871f204e967a4b195262839139a36564a52` | `6f549d3a61f2ba5c1a45e3c28f4c2871f204e967a4b195262839139a36564a52` | 0 |
| `run_log_gdc_step2.txt` | `2e26d4ff92f159ee1b4237501969caa1e0d49859759ae140977a021ab3299719` | `2e26d4ff92f159ee1b4237501969caa1e0d49859759ae140977a021ab3299719` | 0 |
| `run_log_gdc_step3.txt` | `4146154a85d433a20a426b7639e04566e16f9991af86c59d956821926ea9f1f3` | `4146154a85d433a20a426b7639e04566e16f9991af86c59d956821926ea9f1f3` | 0 |
| `run_log_protein_selftest.txt` | `4b163144b261b275e42fa45b672c6c31cc4021be6b0e1e6f8988799bfae7d3e2` | `4b163144b261b275e42fa45b672c6c31cc4021be6b0e1e6f8988799bfae7d3e2` | 0 |
| `run_log_protein_step1.txt` | `384bf933e87bca7b655f82ab83dde816ffffb4c6473fe0b79b5e6eb0c16a2af9` | `384bf933e87bca7b655f82ab83dde816ffffb4c6473fe0b79b5e6eb0c16a2af9` | 0 |
