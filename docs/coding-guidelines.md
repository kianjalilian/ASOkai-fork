# ASOkai Coding & Design Guidelines

## 1. Purpose and Scope

- **Goal**: Define how to extend ASOkai (new analyses, targets, data types) in a **modular**, **reusable**, and **serializable** way.
- **Scope**: This document covers:
  - New analysis **stages**, **groups**, and **workflows** in the ASOkai pipeline (see [§5 Pipeline Contracts](#5-pipeline-contracts-stages-groups-workflows)).
  - New **domain models** (targets, sites, pangenomes, etc.) (see [§4 Domain Contracts](#4-domain-contracts-targets-sites-sequences)).
  - **Serialization** rules (what must be serializable and how) (see [§6 Serialization](#6-serialization)).
- **Audience**: Developers adding or modifying analysis logic, domain types, or pipeline structure in ASOkai.

---

## 2. High-Level Architecture

- **Module layout (domain & antisense)**:
  - Targets and related domain types live in `ASOkai.targets` (e.g. `src/ASOkai/targets/target.py`, `target_gene.py`, `target_gene_creator.py`).
  - Sites and their variants live in `ASOkai.sites` (e.g. `src/ASOkai/sites/genomic_site.py`, `transcript_site.py`).
  - Antisense constructs and ASOs live in `ASOkai.antisense` (e.g. `src/ASOkai/antisense/aso.py`, `antisense_construct.py`).
  - Analysis code, biochemistry, and other domain logic live in their respective packages under `src/ASOkai/` (e.g. `analysis/`, `biochemistry/`).
- **Configuration model**:
  - A single **configuration file** (Default at `config.yaml`) contains all configurable parameters for workflows, groups, and stages.
  - Some configuration keys/arguments are **mandatory** (the pipeline will fail fast if they are missing), others are **optional** with documented defaults.
  - All configuration values from the file can be **overridden or added** via CLI arguments; the CLI layer is always the final source of truth.
- **CLI structure** (see `README.md` for details):
  - `asokai run <workflow_name>`:
    - Positional `<workflow_name>` is interpreted as the workflow to run, e.g. `asokai run standard`.
    - Alternatively, you can specify the workflow explicitly:
      - `--workflow <workflow_name>`: run a predefined workflow (a complete pipeline).
    - Advanced usage:
      - `--groups <group1,group2,...>`: run one or more analysis groups (in order).
      - `--stages <stage1,stage2,...>`: run one or more individual stages (in order).
    - If both a positional workflow and `--groups` / `--stages` are provided, the explicit flags take precedence.
  - Utility/introspection commands:
    - `asokai list workflows|groups|stages`: list available workflows, groups, or stages.
    - `asokai describe workflow|group <name>`: show details about a workflow or group.
    - `asokai config validate -c config.yaml`: validate a configuration file (mandatory/optional keys).
    - `asokai config show -c config.yaml`: show the fully resolved configuration (including defaults and overrides).
- **Conceptual layers**:
  - **Domain layer**: targets, sites, genome/pangenome representations (see [§4](#4-domain-contracts-targets-sites-sequences)).
  - **Analysis layer**: intrinsic features, off-target analyses, kinetic models, etc.
  - **Pipeline layer**: stages, groups, workflows, CLI, configuration, orchestration (see [§5](#5-pipeline-contracts-stages-groups-workflows)).
  - **Serialization layer**: reading/writing objects and results (files, tables, etc.) (see [§6](#6-serialization)).

Describe briefly here how these layers map to the actual packages/modules in `src/ASOkai/` (e.g. `targets/`, `sites/`, `analysis/`, `utils/`, …).

---

## 3. Extension Points

The main extension points in ASOkai are:

- **New workflows**: Predefined combinations of groups/stages for common analyses (see [§5.4](#54-workflows)).
- **New groups**: High-level analyses that orchestrate multiple stages (see [§5.2](#52-groups)).
- **New stages**: Low-level building blocks, each doing one focused task (see [§5.1](#51-stages)).
- **New domain models**: New types for targets, sites, pangenomes, or other biological entities (see [§4](#4-domain-contracts-targets-sites-sequences)).
- **New serializable result types**: New metrics or result tables written to disk (see [§6](#6-serialization)).

For each change, decide which of these you are extending and follow the relevant sections below.

---

## 4. Domain Contracts (Targets, Sites, Sequences)

### 4.1 Target Objects

- **Responsibilities**:
  - Represent a **biological target** for ASO design (e.g. a gene, transcript, locus) and its relevant context.
  - Serve as the primary input for most analysis groups and stages.
- **Required properties/methods** (adapt this list to your actual classes):
  - A stable identifier (e.g. `id`, `gene_id`, `transcript_id`, or similar).
  - Access to underlying genome / transcriptome coordinates where applicable.
  - Access to **sites** associated with the target (e.g. repeated sites, accessible sites, etc.) (see [§4.2 Site Objects](#42-site-objects)).
  - `to_dict()` / `from_dict()` or equivalent serialization hooks (see [§6.3](#63-rules-for-new-serializable-types)).
- **Design rules**:
  - Avoid mixing I/O logic into target classes (loading/saving is handled elsewhere).
  - Keep invariants clear (e.g. coordinates always in genomic 5'→3' orientation, consistent reference genome).

### 4.2 Site Objects

- **Examples**: genomic sites, transcript sites, or other binding site concepts.
- **Required fields** (adapt to your real types):
  - Coordinate identifiers and bounds appropriate to the site’s coordinate system (e.g. chromosome + start/end/strand for genomic coordinates; transcript_id + offsets for transcript coordinates).
  - `sequence: str`.
- **Behavior**:
  - Provide a consistent `.sequence` interface (see [§4.3 Sequences](#43-sequences)).
  - Provide serialization into primitive types (`to_dict()`) (see [§6.3](#63-rules-for-new-serializable-types)).
- **Design rules**:
  - New site types should **reuse** common logic where possible (inherit from a base class or use mixins).
  - Do not assume a particular genome implementation internally; rely on a clear interface.

### 4.3 Sequences

- **Canonical type**: `str` is used to represent sequences. Sites expose this via their `.sequence` interface (see [§4.2](#42-site-objects)).
- **Encoding rules**:
  - Uppercase characters.
  - Allowed alphabet: specify here (e.g. `A`, `C`, `G`, `T`, `U`, `N`).
- **Helpers**:
  - Use a **central helper** for reverse complement and validation (e.g. `utils.sequence.*`).
- **Design rules**:
  - Do not reimplement reverse-complement or validation logic in analysis code; call the shared utilities.

---

## 5. Pipeline Contracts (Stages, Groups, Workflows)

### 5.1 Stages

Stages are **low-level units** in the pipeline, usually corresponding to the `stages` section in the CLI.

- **Responsibilities**:
  - Perform one well-defined operation (e.g. download genome, build genome, compute one type of metric).
  - Read inputs from known locations (config, serialized files, etc.).
  - Write outputs in a consistent, serializable format (see [§6](#6-serialization)).
- **Typical layout**:
  - **Core logic** in a pure/mostly-pure function (e.g. `analysis.intrinsic_features.compute(...)`).
  - **Stage wrapper** that:
    - Reads configuration.
    - Loads required domain objects.
    - Calls the core logic.
    - Serializes results.

- **Design rules**:
  - Stages should be **idempotent** where practical:
    - Re-running should either overwrite outputs or detect existing outputs and skip with a clear rule.
  - Clearly document:
    - Inputs (files, objects, config keys).
    - Outputs (files and their schema).
  - Avoid side effects (like temporary global state) that other stages rely on implicitly.

### 5.2 Groups

Groups are **logical collections of stages** that perform a complete analysis.

- **Responsibilities**:
  - Orchestrate multiple [stages](#51-stages) in the correct order.
  - Handle dependencies (target gene object, genome, etc.) (see [§4.1](#41-target-objects)).
  - Provide a higher-level CLI entry point.
- **Dependency handling**:
  - Each group declares:
    - Required **stages** it depends on.
    - Required **groups** it depends on (if any).
  - On execution:
    - If dependencies are already satisfied (files/objects exist), re-use them.
    - If not, run the necessary preceding stages/groups.

- **Design rules**:
  - Groups mostly coordinate; **heavy logic belongs in stages/analysis modules**.
  - The same analysis functionality should be usable from Python API, not only via CLI.
  - Clearly document the **expected preconditions** (e.g. “requires instantiated target gene”).

### 5.3 Dependency Contracts (Draft)

- Every [stage](#51-stages) must declare:
  - Its required **inputs** (artifacts, config keys, domain objects) in a way that can be checked.
  - Its produced **outputs** (artifacts, updated domain objects, result tables).
- Every [group](#52-groups) must ensure that, within the group:
  - For each stage in its internal sequence, all declared inputs are satisfied by:
    - Outputs of earlier stages in the group, and/or
    - Artifacts guaranteed by the group’s own declared dependencies.
- Every **workflow** must be a **clean pipeline**:
  - For each group or stage in the workflow order, all of its declared inputs must be satisfied by:
    - Outputs of earlier steps in the same workflow, and/or
    - Global prerequisites of the workflow (e.g. “requires target gene instantiated”).
- Hidden or implicit dependencies should be avoided; we will later make this contract explicit and (ideally) machine-checkable.

### 5.4 Workflows

Workflows are **named, high-level pipelines** that combine groups and/or stages into commonly used sequences.

- **Responsibilities**:
  - Provide a single entry point for complex multi-step analyses.
  - Encapsulate ordering and configuration of groups/stages.
- **Design rules**:
  - A workflow should be a **thin composition layer** on top of groups and stages.
  - Prefer reusing existing groups rather than duplicating logic.
  - Clearly document:
    - Which groups/stages it runs.
    - What inputs it requires (config keys, existing artifacts).

---

## 6. Serialization

### 6.1 What Must Be Serializable

At minimum, the following should be serializable:

- **Domain objects**:
  - Target gene objects (see [§4.1](#41-target-objects)).
  - Site collections (see [§4.2](#42-site-objects)).
- **Analysis results**:
  - Per-site metrics (intrinsic features, accessibility scores, off-target counts, etc.).
  - Per-target metrics (summary statistics, kinetic model outputs).
- **Pipeline artifacts**:
  - Any intermediate result you expect to re-use across stages/groups.

### 6.2 How We Serialize

- **Central module**: All serialization should go through `ASOkai.utils.serializer` (adapt name if different).
- **Allowed formats** (define for this project, e.g.):
  - JSON for small structured configs/metadata.
  - CSV/Parquet for tabular per-site or per-target results.
  - Other formats only via clearly documented helpers.
- **Schema and versioning**:
  - Each serialized artifact should have a **schema** (columns, types, units) and a **version**.
  - When changing schema:
    - Bump the version.
    - Keep backwards-compatible loaders when feasible or provide migration utilities.

### 6.3 Rules for New Serializable Types

- Implement `to_dict()` / `from_dict()` (or an equivalent protocol) on domain classes that need serialization (see [§4.1](#41-target-objects), [§4.2](#42-site-objects)).
- Use only **primitive types** in public serialized form:
  - Strings, numbers, arrays, and simple dicts.
  - Map enums and custom Python types to strings.
- Do **not** rely on `pickle` for long-term storage or cross-version compatibility.
- Any analysis that writes files must:
  - Use shared serializer utilities (e.g. `serializer.save_*`).
  - Document the output schema and version in the function/stage docstring.

---

## 7. Coding & Design Patterns for Analyses

### 7.1 General Principles

- **Separation of concerns**:
  - Domain logic inside domain/analysis modules.
  - I/O and CLI glue code in pipeline/CLI modules.
- **Modularity**:
  - Prefer small, composable functions and classes over monolithic ones.
  - Reuse existing building blocks instead of duplicating logic.
- **Configurability**:
  - All user-facing parameters must be:
    - Representable in `config.yaml`.
    - Passed explicitly into analysis code (no hidden global config).

### 7.2 Functions vs Classes

- Prefer **pure functions** operating on `TargetGene` / `Site` collections for analysis calculations (see [§4.1](#41-target-objects), [§4.2](#42-site-objects)).
- Use small **config/data classes** for grouping parameters.
- Avoid large stateful classes that mix computation, I/O, and configuration.

### 7.3 I/O Boundaries

- Core analysis functions:
  - Accept in-memory data structures (targets, sites, tables).
  - Return in-memory results (tables, lists of records, etc.).
- Stages and CLI:
  - Are responsible for loading inputs (via serializer) and writing outputs (see [§6](#6-serialization)).
  - Pass only domain objects and primitive data into analysis functions.

### 7.4 Configuration Patterns

- All tunable parameters:
  - Must have a documented key in `config.yaml`.
  - Should be mapped into typed config objects where practical (e.g. `IntrinsicFeaturesConfig`).
- Distinguish **mandatory** vs **optional** configuration:
  - Mandatory keys must be validated early (ideally at CLI/config loading time) with clear error messages.
  - Optional keys must have documented defaults and behavior when omitted.
- CLI arguments:
  - May override existing keys from `config.yaml`.
  - May introduce additional keys, which should still be validated and documented.
- Avoid magic constants in analysis logic:
  - Put them in configuration with clear default values and descriptions.

### 7.5 Reusability and API Design

- Any analysis that is exposed via CLI should also be callable from Python:
  - e.g. `ASOkai.analysis.intrinsic_features.compute(target_gene, config)`.
- Document public Python entry points so they can be reused in external workflows or notebooks.

---

## 8. Testing Guidelines

### 8.1 Test Layout

- **Unit tests**:
  - Location: `tests/unit/`.
  - One test module per main module or closely related set of functions.
- **Integration tests**:
  - Location: `tests/integration/`.
  - Cover realistic workflows (e.g. running a group end-to-end with a small example dataset).

### 8.2 Expectations for New Code

- **New stages** (see [§5.1](#51-stages)):
  - At least one **unit test** for the core analysis logic.
  - If they change the pipeline behavior, add or update **integration tests**.
- **New groups/workflows** (see [§5.2](#52-groups), [§5.4](#54-workflows)):
  - At least one **integration test** covering the full group/workflow on a small test input.

### 8.3 Test Patterns

- Use fixtures for:
  - Small genomes or genome fragments.
  - Example target genes (see [§4.1](#41-target-objects)).
  - Synthetic sites/metrics (see [§4.2](#42-site-objects)).
- Prefer tests that are:
  - Deterministic.
  - Fast enough to run frequently in CI and locally.

---

## 9. How to Add a New Analysis (Worked Template)

### 9.1 Example: New Intrinsic Feature (Per-Site Metric)

1. **Define the metric**:
   - Input: site (sequence, coordinates, etc.) (see [§4.2](#42-site-objects), [§4.3](#43-sequences)).
   - Output: numeric or categorical value (e.g. GC skew, motif score).
2. **Implement core logic** in the relevant `analysis/` module:
   - A pure function (or small set of functions) that consumes:
     - A `TargetGene` or site collection (see [§4.1](#41-target-objects), [§4.2](#42-site-objects)).
     - A config object or simple parameters.
   - Returns a table or list of records containing:
     - Site identifier.
     - Metric value(s).
3. **Integrate with serialization** (see [§6](#6-serialization)):
   - Convert results to a serializable structure (e.g. DataFrame, list of dicts).
   - Save via `utils.serializer` using a new or existing schema.
4. **Create/extend a stage** (see [§5.1](#51-stages)):
   - Define a new stage (e.g. `stages.compute-new-feature`) or extend an existing feature stage.
   - Stage should:
     - Load the target gene object.
     - Read configuration.
     - Call your analysis function.
     - Serialize results.
5. **Optionally add a group or workflow** (see [§5.2](#52-groups), [§5.4](#54-workflows)):
   - If this feature is part of a broader analysis, create or extend a group/workflow.
   - Declare dependencies and ensure it plays nicely with existing groups/stages.
6. **Add tests** (see [§8 Testing Guidelines](#8-testing-guidelines)):
   - Unit tests on the core analysis logic.
   - If necessary, an integration test on the group/stage/workflow.

---

## 10. Contributor Checklists

### 10.1 Adding a New Workflow

- [ ] Define workflow name, description, and the groups/stages it runs (see [§5.4](#54-workflows)).
- [ ] Ensure it only orchestrates existing groups/stages (no hidden logic).
- [ ] Document required inputs and configuration.
- [ ] Add or update integration tests that cover the workflow (see [§8](#8-testing-guidelines)).

### 10.2 Adding a New Group

- [ ] Define CLI name and description.
- [ ] Declare group dependencies (other groups/stages) (see [§5.2](#52-groups)).
- [ ] Implement group orchestration (calls to stages).
- [ ] Ensure results are serialized using `utils.serializer` (see [§6](#6-serialization)).
- [ ] Add/update integration tests (see [§8](#8-testing-guidelines)).
- [ ] Document group behavior and outputs (README or docs).

### 10.3 Adding a New Stage

- [ ] Implement core logic in an appropriate analysis/domain module (see [§5.1](#51-stages)).
- [ ] Implement a stage wrapper for config + I/O.
- [ ] Register the stage with the CLI.
- [ ] Add unit tests for the core logic (see [§8](#8-testing-guidelines)).
- [ ] If needed, add integration tests.

### 10.4 Adding a New Domain/Data Type (e.g. Pangenome)

- [ ] Choose module/package location.
- [ ] Define its relationship to existing types (e.g. `GenomeUtils` objects) (see [§4](#4-domain-contracts-targets-sites-sequences)).
- [ ] Specify required methods/properties (contracts).
- [ ] Define serialization format and version (see [§6](#6-serialization)).
- [ ] Add tests.
- [ ] Update relevant sections of this document (architecture and domain contracts) (see [§2](#2-high-level-architecture), [§4](#4-domain-contracts-targets-sites-sequences)).

---

## 11. Maintaining This Document

- **Owner**: _<name/role>_.
- **Update policy**:
  - Any new extension point (new kind of workflow/group/stage/domain type) should be reflected here.
  - Changes to serialization formats or major APIs must be documented.
- **Process**:
  - Update this Markdown file via Pull Requests.
  - Keep examples in sync with actual code (prefer links to specific modules/tests where possible).