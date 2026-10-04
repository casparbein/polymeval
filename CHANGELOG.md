# Changelog

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