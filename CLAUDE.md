# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Reusable **components** for the SMaHT-DAC variant calling pipelines. This repo contains only building blocks — the individual tool steps. The MetaWorkflows (MWF/MWFR) that chain these steps into full pipelines live in a **separate `main-pipelines` repo**, and shared file-format definitions live in a **`shared_pipelines` repo**. Do not add MetaWorkflow definitions here; they were deliberately moved out (see git history: "Move MWF to main-pipelines", "Moved file formats to shared_pipelines").

These components are consumed by the SMaHT portal via the `portal-pipeline-utils` / `smaht_pipeline_utils` tooling, which reads the `portal_objects/` YAMLs and `descriptions/` CWL to register pipelines.

There is no build system, linter, or test suite in this repo — no Makefile, CI config, or package manifest. Changes are validated by building the Docker images and by the downstream portal/pipeline deployment tooling.

## Architecture: three synchronized layers per tool

Each variant caller / utility is defined across three parallel locations that **must be kept consistent**. When adding or changing a step, update all three:

1. **`dockerfiles/<tool>/`** — a `Dockerfile` plus the bash entrypoint script(s) it installs (e.g. [dockerfiles/delly/run_delly_sr.sh](dockerfiles/delly/run_delly_sr.sh)). The Dockerfile `COPY`s each script into the image's PATH and `chmod +x`s it. Shared/utility scripts (bcftools normalization, minipileup, filtering, phasing, caller-merging) all live under [dockerfiles/calling_utils/](dockerfiles/calling_utils/), with Python helpers in [dockerfiles/calling_utils/scripts/](dockerfiles/calling_utils/scripts/). Base images come from `public.ecr.aws/smaht-dac/base-ubuntu2204-py3x`.

2. **`descriptions/<step>.cwl`** — a CWL `CommandLineTool` (cwlVersion v1.0) wrapping one bash entrypoint. The `baseCommand` is the bare script name (e.g. `[run_delly_sr.sh]`), which resolves because the Dockerfile puts it on PATH. The Docker image is referenced as a placeholder `dockerPull: ACCOUNT/<tool>:VERSION` — `ACCOUNT` and `VERSION` are substituted by the deployment tooling, not hardcoded. `inputBinding.prefix` values (e.g. `-s`, `-r`, `-x`) map directly to the `getopts` flags in the corresponding bash script.

3. **`portal_objects/workflows/<step>.yaml`** — portal Workflow metadata. Declares `runner.main` (the CWL file), `software` (as `Name@version`, matching `software.yaml`), and typed `input`/`output` using `argument_type: file.<format>` with `secondary_files`. The input/output arguments here must line up with the CWL `inputs`/`outputs` ids.

Not every CWL has a portal workflow YAML (there are more `.cwl` descriptions than workflow YAMLs); some CWL steps are internal sub-steps composed inside MetaWorkflows in the `main-pipelines` repo.

### Portal object registries (`portal_objects/`)

- **`software.yaml`** — one YAML document per tool version (`---` separated). Fields: `name`, `version`, `source_url`, `category`, and `code` (required for the version to appear in output file names).
- **`file_reference.yaml`** — reference/resource files (VEP cache, gnomAD, region BEDs, PON, tandem-repeat regions, etc.), one document each. Each carries a `uuid` and `accession` that enable sync with the reference S3 bucket — treat these as stable identifiers; do not regenerate or change them for an existing file.

## Conventions

- **Bash entrypoints** use `set -euo pipefail`, parse args with `getopts`, validate required args, and wrap each pipeline command with `|| { echo "Error: ..."; exit 1; }`.
- **Adding a tool version:** bump it in the Dockerfile install step, add/update the `software.yaml` document, and update the `@version` reference in any workflow YAML. Keep the `code:` field when the tool name must appear in output filenames.
- **File formats** (`file.cram`, `file.vcf_gz`, etc.) and their secondary-file suffixes follow the shared conventions defined in the `shared_pipelines` repo.
