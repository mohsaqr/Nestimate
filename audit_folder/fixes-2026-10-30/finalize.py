from pathlib import Path
import csv,json,runpy,html
from html.parser import HTMLParser
p=Path(__file__).resolve().parent
root=p.parent.parent
rows=list(csv.DictReader((p/'final-full-tests.csv').open()))
number=lambda x: 1 if x=='TRUE' else 0 if x=='FALSE' else int(x)
totals={k:sum(number(x[k]) for x in rows) for k in ['passed','failed','error','warning','skipped']}
assert totals['failed']==totals['error']==0
assert 'Status: OK' in (p/'final-check.log').read_text()
assert 'htna R CMD check is clean' in (p/'revdep-htna-final.log').read_text()
new=[r for r in rows if r['file'] in ['test-fix-review-sequences.R','test-fix-review-models.R','test-fix-review-higherorder.R','test-fix-review-landscapes.R']]
summary=dict(totals,files=len({r['file'] for r in rows}),test_blocks=len(rows),new_passed=sum(number(r['passed']) for r in new),package_check='Status: OK; 0 errors, warnings or notes',htna_check='Status: OK; 0 errors, warnings or notes',visual_review='No visual approval is claimed.')
(p/'validation-summary.json').write_text(json.dumps(summary,indent=2))
res=json.loads((p/'resolutions.json').read_text())
res['R08']['resolution']+=' Surface stability measures retain every prepared position, including start markers.'
for r in res.values():r['status']='fixed; final full-suite regression passed'
(p/'resolutions.json').write_text(json.dumps(res,indent=2))
runpy.run_path(str(p/'update-report.py'),run_name='__main__')
report=(p.parent/'review-2026-10-30'/'report.html').read_text()
findings=json.loads((p.parent/'review-2026-10-30'/'findings.json').read_text())
assert len(findings)==12 and all(f['resolution_status']=='fixed' for f in findings)
assert 'VALIDATION_PLACEHOLDER' not in report
assert all('id="'+f['id']+'"' in report for f in findings)
parser=HTMLParser();parser.feed(report);parser.close()
entry=f'''\n### 2026-10-30 (review fixes completed)
- [review closure]: all 12 findings and the case-drop target concern are resolved. Original review is preserved as audit_folder/review-2026-10-30/report.before-fixes.html; report.html and findings.json now record resolutions. New evidence is in audit_folder/fixes-2026-10-30/.
- [estimator replay]: landscapes prepare concatenation/markers and retain full-sequence weighting before counts and bootstrap units; surface metrics must read all prepared positions (unname prevents start markers being inferred as metadata). Covariates must be constant within merged units. Unreconstructable transformed networks error explicitly.
- [downstream shape]: htna requires the original case-drop result field names. Store resolved target in attr(result, "target") rather than a new list field. Default target includes all cells for association/transformed networks; explicit existing/all remain available.
- [inference contracts]: composite bootstraps retain/replay custom weighting; older saved custom fits need rebuilding. Post-hoc clustering covariates are excluded before distances. Gaussian ordinary-fit intervals use residual-df t limits. Named graph comparisons align axes.
- [final validation]: {summary['passed']} passed, 0 failures/errors, {summary['warning']} warnings, {summary['skipped']} skips across {summary['files']} files; {summary['new_passed']} new regression assertions. Final --as-cran --no-manual and isolated htna 0.3.1 checks: Status OK. Archive matches all 28 changed packaged files byte for byte. CI is configured for Linux/macOS/Windows but not run remotely here.
'''
l=root/'LEARNINGS.md';s=l.read_text()
if '### 2026-10-30 (review fixes completed)' not in s:l.write_text(s.rstrip()+'\n'+entry)
c=root/'CHANGES.md';s=c.read_text()
c.write_text(f'''### 2026-10-30 — fixed all package-review findings
- R/build_network.R, pattern_outcome.R, sequence_patterns.R: numeric missing positions, metadata-safe generated outcomes, documented classed rejection of wildcard-only state labels.
- R/cluster_data.R, mcml_pc.R, compare_network.R, outcome_model.R, centrality_measures.R: prevent covariate leakage, replay composite weighting, align named matrices, correct Gaussian intervals, normalize small strengths correctly.
- R/hypergraph_centrality.R, mogen.R, casedrop_reliability.R: log-space tensor directions, numeric DOF summaries, explicit case-drop targets with backward-compatible result fields.
- R/stability_landscape.R: replay weighting/concat/markers, resample merged units, validate merged covariates and reconstructed weights, keep boundary positions in metrics, fall back to full-range grids for tied trimmed endpoints.
- Four test-fix-review regression files: {summary['new_passed']} assertions; relevant Rd help and NEWS synchronized. _pkgdown.yml includes all five missing topics; .github/workflows/R-CMD-check.yaml configures three-platform full testing/checks.
- Updated review HTML/JSON with all resolutions; retained before-fix report, code patch, source-archive hashes, test/check logs and rendered plot preview tmp/review-fixes-landscapes.html.
- Tests: final full suite {summary['passed']} passed / 0 failures/errors / {summary['warning']} warnings / {summary['skipped']} skips. Final CRAN-style --no-manual and isolated htna checks: Status OK. No commits or pushes.

'''+s)
(root/'HANDOFF.md').write_text(f'''# Session Handoff — 2026-10-30

## Completed
- User requested fixing every issue in the extensive package review. All 12 defects and the separate case-drop target concern are fixed, documented and regression-tested.
- Current report: audit_folder/review-2026-10-30/report.html; original: report.before-fixes.html. findings.json carries fixed status and regression references; resolutions.json and full evidence are under audit_folder/fixes-2026-10-30/.
- Production fixes span 12 R files; synchronized help/NEWS; four regression files add {summary['new_passed']} passing assertions. Added missing pkgdown topics and three-platform CI.
- No commits, pushes, dependency installs/upgrades or sibling-repository modifications.

## Current State
- Version remains 1.0.1 above dcdd7135, CRAN baseline 0.9.24 / 74617a5e.
- Final complete development run outside sandbox: {summary['passed']} passed, 0 failures, 0 errors, {summary['warning']} warnings, {summary['skipped']} skips, {summary['files']} files / {summary['test_blocks']} blocks. See final-full-tests.csv/.rds/.log.
- Final source build with vignettes succeeded; final R CMD check --as-cran --no-manual Status OK. Isolated htna 0.3.1 R CMD check also Status OK. Final archive byte-matches all 28 changed packaged files.
- All 12 structured findings marked fixed. Report HTML, JSON, anchors, docs and CI syntax verified.
- Rendered plot preview: tmp/review-fixes-landscapes.html, presented through an optional visual-review question. No visual approval claimed.

## Key Decisions
- Case-drop target appended to signature; auto selects existing edges only for unscaled/unthresholded transition networks, all cells otherwise. Resolved target is an ATTRIBUTE to preserve htna's exact result fields. First reverse-dependency run caught a new-field break; corrected and rechecked clean.
- Patterns reject stars-only literal labels explicitly (wildcards reserved). Generated outcomes never overwrite metadata.
- Custom composite weighting specification is retained for new fits; old custom fits lacking it require rebuild.
- Landscapes replay estimator settings and use merged bootstrap units; conflicting per-unit covariates and transformed/pruned weights that cannot be reconstructed are refused. Surface measures use every prepared position, including boundary markers.
- Report keeps original defect evidence separate from resolution status.

## Open Issues
- No outstanding defects from the requested review.
- Existing diagnostic/deprecation warnings remain in development tests; see logs.
- Windows/Linux workflow is configured but not executed remotely in this session. Minimum dependency versions, full pkgdown build and PDF manual not validated. Local runtime has development idiographic and newer cograph.

## Next Steps
1. Review changes and optional plot preview; commit only if explicitly authorized.
2. Before release, inspect remote CI and validate minimum dependency versions and PDF manual as needed.

## Context
- R 4.5.2/macOS Tahoe; psychnets 0.5.2, idiographic 0.4.0.9004, cograph 2.7.11, tna 1.3.1.
- processx/lavaan need approved unsandboxed execution here. Full final run/checks used escalation; no sandbox errors remain in final totals.
- Parallel agents used under inherited AGENTS.md. Independent landscape recheck exposed an additional boundary-position case, now fixed with a numeric oracle.
- Pre-existing untracked user files preserved. audit_folder artifacts and session notes are not committed. Final changed-file hash manifest and implementation.patch retained under fixes evidence.
''')
print(json.dumps(summary,indent=2))
print('Report, manifest and session records validated and finalized.')
