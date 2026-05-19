# Changelog

All notable changes to Venus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-05-19

### Added
- Initial release of Venus as standalone repository
- Three analysis modes: experiment_only, control_only, experiment_control
- Support for WGS, WES, and Panel sequencing
- Modular workflow files using `[[include]]` directive
- Conda environment definitions
- CLI commands: generate, validate, list-steps
- VEP annotation integration
- Clinical report generation

### Changed
- Extracted from oxo-flow workspace into independent repository
- Configuration format updated with `[[samples]]` syntax

### Removed
- Removed from oxo-flow workspace Cargo.toml
