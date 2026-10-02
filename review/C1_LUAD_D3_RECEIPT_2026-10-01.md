# D3 third-party timestamp: receipt

Recorded outside `review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md` per that file's own
§7.3 ("sha256 of this file (recorded outside the file)") and §0 decision D3. This
receipt does not modify, and is not part of, the signed file; the signed file's sha256
(`a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022`) is unchanged by
anything below.

## D3 step 1-2: signature

- Signed by: Wei Zhang
- Date and time, with time zone: 2026-10-01 04:16 CST (UTC+8)
- Signed file: `review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md`
- sha256 of the signed file: `a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022`
  (recorded in the sidecar `review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md.sha256`,
  computed immediately after signature was saved and before any further edit)

## D3 step 3: OSF registration

- Registry: OSF Registries, "Secondary Data Preregistration" template
- Registration: https://osf.io/2d6mb/
- Registration state: embargoed (`public: false`, `embargo_end_date:
  2028-10-01T00:00:00Z`), not archiving, not pending approval, not withdrawn,
  confirmed via the OSF v2 API 2026-10-01
- Date registered (OSF, per the API `date_registered`): 2026-09-30T22:00:24Z UTC =
  2026-10-01 06:00 CST
- Registration DOI: not yet minted by DataCite as of this receipt (OSF creates the
  DOI record asynchronously after approval; the identifier slot exists,
  `id=6abd872aacef6206e3a2780a`, `category=doi`, value pending). **To be added to
  this receipt once it appears** -- check https://osf.io/2d6mb/ (Metadata panel,
  "Registration DOI") or `GET https://api.osf.io/v2/registrations/2d6mb/identifiers/`.
- Attachments on the registration (all verified present with the correct, current
  OSF Storage file IDs before submission): the signed rule, its sha256 sidecar, the
  Step 0 note, the Step 0b workflow summary, the dated errata, and the frozen
  113-case confirmatory ID list (`c1_case_ids_luad.txt`, separately, as the
  Codebook Documentation attachment).
- Associated OSF project (created automatically to host the registration's files):
  https://osf.io/apmv6/

### Correction, recorded for the audit trail

The first submission attempt (registration node `2uv8m`, created 2026-09-30T21:52:59Z)
selected "Make registration public immediately" rather than embargo. This was caught
before approval (the pending-approval email states plainly: "the registration will be
made public"), and that pending registration was cancelled before any admin approved
it -- cancelling reverted the draft to editable state with no loss of its attached
files. The draft was then resubmitted with "Enter registration into embargo" selected,
producing the `2d6mb` registration recorded above. `2uv8m` itself now returns HTTP 410
("no longer available") and carries no content. No OSF DOI or public registration ever
existed under `2uv8m`, so this correction changes nothing about the D3 timestamp basis
-- the timestamp that counts is the one above, from the registration that was actually
approved.

A second, unrelated defect was found and fixed before either submission's final
approval: two of the five files attached to the "Data Collection Procedures
Documentation" question (the signed rule and the Step 0 note) had stale OSF Storage
file-object references left over from an earlier delete-and-reupload cycle during
drafting -- the question's stored reference pointed at a file ID that no longer
matched the current object in storage, which is what OSF's own pre-submission check
flagged as "file(s) ... not part of a component being registered." Confirmed via the
OSF v2 API (`GET /v2/draft_registrations/6abd71d668cf0ed12f54b258/`,
`registration_responses['72-22']`) and fixed by deselecting and reselecting each
affected file in the Data Description step, which updated the stored reference to the
correct, current file ID. Verified after the fix: every `file_id` referenced by
question `72-22` and `72-26` matched the corresponding file's current ID in the
node's OSF Storage tree (`GET /v2/nodes/apmv6/files/osfstorage/`), and the rule
file's hash in that storage matched the signed sha256 above. This defect existed only
in the draft-registration form's internal bookkeeping; no confirmatory LUAD value was
touched by it, and it does not affect the timestamp basis.

## D3 step 4: GitHub Release

- Repository: https://github.com/tiandaochouqin-Wei/MorphoResidual (public)
- Commit: `18c6bfb1ed33b59c9d9b959aac780d1b6af875be`, "Ship the signed C1-LUAD
  confirmatory-test rule (D3 third-party timestamp)"
- Files added at a tracked path (a `.gitignore` exception carved out of the
  otherwise-private `review/*`, matching the existing pattern used for
  `review/recalc/`): `review/C1_LUAD_FROZEN_RULE_DRAFT_2026-09-27.md` and its
  `.sha256` sidecar. Confirmed byte-identical to the signed copy: the sha256 of the
  git blob at commit time was verified to equal
  `a68dfbec65e0186bb629eaae1f81c8e117549e224b39483d67513c9d1270d022` before the
  commit was made.
- Pushed to `origin/master`: `0588b5f..18c6bfb`
- Release tag: `c1-luad-signed-rule-2026-10-01`
- Release URL:
  https://github.com/tiandaochouqin-Wei/MorphoResidual/releases/tag/c1-luad-signed-rule-2026-10-01
- Published at (GitHub, UTC): 2026-09-30T22:12:54Z = 2026-10-01 06:12:54 CST
- Tag target, verified locally after `git fetch origin tag
  c1-luad-signed-rule-2026-10-01`: commit `18c6bfb1...`, matching the commit above
  exactly.

## D3 step 6: confirmatory data may now be downloaded

Both the OSF registration (step 3) and the GitHub Release (step 4) exist, each with
a durable, externally-verifiable timestamp, before any confirmatory LUAD slide, RNA,
or protein value has been opened. Per the signed rule's own D3 clause, confirmatory
data may be downloaded from this point on; PDC/GDC/IDC download logs are to be kept
as the rule requires.

**Outstanding:** the OSF Registration DOI has not yet been minted; when it appears,
add it to this receipt (it does not require re-signing or re-hashing anything, since
neither the signed rule nor this receipt's existing content changes).
