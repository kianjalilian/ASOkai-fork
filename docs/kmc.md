# `kmc` and `kmc_tools` in ASOkai

Small Python helpers around the [KMC](https://github.com/refresh-bio/KMC) counter (`kmc`) and the `kmc_tools` utility. Requires the upstream binaries on `PATH`, or an explicit path passed to the constructor.

KMC is third-party GPLv3-only software. See the project [README](../README.md) (third-party software section).

## Imports

```python
from ASOkai.Utils import KMC
from ASOkai.Utils.KMCTools import Simple, Transform, Filter, Complex
```

## `KMC` — running `kmc`

- **Constructor:** `KMC(kmc_executable: str = "kmc")` — default executable name on `PATH`; a filesystem path selects a specific binary.
- **Main API:** `run(input_path, output_db_prefix, working_directory, **options)` — returns a `subprocess` completed-process object (`returncode`; `stdout` / `stderr` captured when `debug=False`).

Arguments follow the `kmc` CLI: input file or `@list_file`, output database prefix (paths to `{prefix}.kmc_pre` / `{prefix}.kmc_suf`), working directory for temporary files.

Parameters correspond to [KMC 3.2.x](https://github.com/refresh-bio/KMC) `kmc` options. Full keyword list: `KMC.run` docstring.

| Parameter | Description |
|-----------|-------------|
| `k` | K-mer length (1–256); default **25** |
| `m` | Max RAM in GB (1–1024); default **12** |
| `f` | Format: `fa`, `fq`, `fm`, `fbam`, `fkmc`. Always emitted as `-f…`; default **`fm`**. Standalone `kmc` uses FASTQ when `-f` is omitted. |
| `sm` | Strict memory mode (`-sm`) |
| `t` | Threads (`-t`); omitted → KMC default (all CPU cores) |
| `p` | Signature length 5–11; omitted → KMC default **9** |
| `b` | **`True`:** pass `-b` (non-canonical k-mers). **`False` (default):** omit `-b` (canonical k-mers). |
| `hc` | Homopolymer-compressed k-mers (`-hc`) |
| `r` | RAM-only mode (`-r`) |
| `ci`, `cs`, `cx` | `-ci`, `-cs`, `-cx`; omitted → **2**, **255**, **1e9** |
| `o` | `kmc` or `kff`; omitted → **kmc** |
| `w` | No database output (`-w`) |
| `e` | Histogram estimate only (`-e`) |
| `hp` | Hide progress (`-hp`) |
| `v` | Verbose (`-v`) |
| `j` | JSON summary path (`-j`) |
| `sf`, `sp`, `sr` | FASTQ read / split / stage-2 thread counts |
| `n` | Bins (`-n`) |
| `opt_out_size` | `--opt-out-size` |
| `debug` | `True`: DEBUG logging, stream subprocess output; `False`: capture output |
| `check` | `True`: non-zero exit → `KMCExecutionError` |
| `additional_args` | Extra command-line tokens appended after built-in options |

Input paths resolve to absolute paths. List inputs use the `@` prefix; the path after `@` is resolved.



---

## `kmc_tools` wrappers

`Simple`, `Transform`, `Filter`, `Complex`: executable resolution matches `KMC` (`kmc_tools` on `PATH` or constructor path). Failed runs with `check=True` expose `returncode`, `cmd`, `stdout`, `stderr` on the raised exception.

Shared options:

| Parameter | Description |
|-----------|-------------|
| `t` | Threads (`-t`) |
| `v` | Verbose (`-v`) |
| `hp` | Hide progress (`-hp`) |
| `debug` | Log command; stream vs capture |
| `check` | Raise on non-zero exit |
| `cwd` | Subprocess working directory |

### `Transform`

```python
Transform().reduce("err_kmers", cx=10).reduce("valid_kmers", ci=11).histogram("histo.txt").dump("dump.txt").run("db", t=8)
```

`sort`, `reduce`, `compact`, `histogram`, `dump`, `set_counts`; then `run(input_prefix, …)` with shared options from the table above.

### `Simple`

```python
Simple().union("kmers1_kmers2_union", cs=65536, oc="left").intersect("intersect_kmers1_kmers2").intersect("intersect_max_kmers1_kmers2", oc="max").run("kmers1", "kmers2", input1_ci=3, input1_cx=70000, t=8)
```

`intersect`, `union`, `kmers_subtract`, `counters_subtract`, `reverse_kmers_subtract`, `reverse_counters_subtract`; each takes `output` and optional `ci`, `cx`, `cs`, `o`, `oc`. `run(input1, input2, …)` adds `input1_ci` / `input1_cx` / `input2_ci` / `input2_cx` and shared options.

### `Filter`

```python
Filter().run("kmc_db", "input.fastq", "filtered.fastq", db_ci=3, read_ci=0.5, read_cx=1.0, t=8)
```

`trim` / `hm` (filter-level `-t` / `-hm`); `db_ci` / `db_cx` on the KMC DB; `read_ci` / `read_cx` / `read_f` on input reads (`read_ci` / `read_cx` may be int or fraction in `[0, 1]`); `output_f` on the output (omit to match input). Thread count is the global `-t` from the shared options table, not `trim`.

### `Complex`

```python
Complex().run("/path/to/ops.txt", t=8)
```

Passes one **operations definition file** path: `INPUT:` (`name = db_path` with optional `-ci` / `-cx`), then `OUTPUT:` (expression with `* - ~ +` for intersect / kmers_subtract / counters_subtract / union, optional `[c_mode]`, optional output `-ci` / `-cx` / `-cs`). Full syntax is in `kmc_tools complex` help.
