# sentinel-findings.md

*Written by Sentinel — WARN-tier findings.*

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1090`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:18.317497+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1091`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:18.317522+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1093`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:18.317531+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1096`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:18.317537+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:18.317542+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1090`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:05.914242+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1091`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:05.914276+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1093`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:05.914287+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1096`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:05.914294+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:05.914301+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1090`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:18.186939+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1091`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:18.186966+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1093`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:18.186975+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1096`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:18.186981+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:18.186987+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1090`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:25:59.559825+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1091`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:25:59.559853+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1093`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:25:59.559864+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1096`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:25:59.559871+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:25:59.559877+00:00

## [WARN] P1-license-noncompliance
- **Location:** `(repo-wide)`
- **Detail:** package 'pyreadr' has license 'GNU Affero General Public License v3 or later (AGPLv3+)' which conflicts with MIT-licensed publication. If this repo is intentionally GPL/AGPL, suppress with a sentinel:skip-line marker; otherwise swap to a permissive-licensed alternative.
- **Fix hint:** Options: (1) replace the dep with an MIT/BSD/Apache-2.0 alternative; (2) change the repo's LICENSE to GPL/AGPL to match the dep; (3) if dual-licensed in your favor, add a comment near the dep declaration with the chosen license and add 'sentinel:skip-line P1-license-noncompliance' nearby.
- **Source:** lessons.md#license-compliance-2026-04-29
- **When:** 2026-05-15T15:39:08.267157+00:00
