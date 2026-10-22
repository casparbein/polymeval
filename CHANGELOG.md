# Changelog

## [0.2.0] - 2026-10-22

### Added
- `--samples` added to reference mode: a run can be restricted to named samples. Useful if 
  a standard run was performed but only some of the reads and assemblies should be used
  for reference
- Combine mode's 'pairwise' sampling extended to six input read sets (15 pairs)
- `seqkit` now removes reads with duplicate names, since flye chokes on this
- bumped compleasm version to 0.2.10, which fixes db download speed

### Changed
- Allocated memory raised for meryl, compleasm, hifiasm and flye (needed for big assemblies)

### Fixed
- Bugfix: sniffles2: `--report_seq` replaced by `--max-del-seq-len` (broke sniffles2 runs before)

## [0.1.0] - 2026-10-02

### Added
- First release.
  - pinned all pinnable packages in src/polymeval/workflow/envs and environment.yaml
  - snakemake wrappers updated to available versions from 2026-10-02 (cli.py)
  - wd-related bindings for deepvariant container

## [0.1.0-pre1] - 2026-10-01

### Added
- First prerelease.
- Extensions from pre-release:
  - standard/downsample:
    - assemblies now possible with LJA, Flye, Verkko
  - HG002 vcf benchmarking:
    - auto-fetch of truth sets implemented (--fetch_benchmarks)
    - new caller: longcallD, new evaluator: aardvark
    - new analysis: allelic imbalance on truth het sites