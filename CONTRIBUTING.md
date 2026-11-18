# Contributing to cfizz

> Chromosome Flow Integration with Zone-Terminal Visualization
> End-to-end Hi-C differential analysis and zone-terminal visualization.

Thank you for your interest in cfizz. This guide covers the development setup, code and documentation conventions, examples, and the pull-request workflow.

## Quick start

```bash
git clone https://github.com/Cuixiaojian21/cfizz.git
cd cfizz

micromamba env create -f environment.yml
micromamba activate cfizz

pip install -r requirements.txt
pip install -e ".[all]"
```

Scripts under `examples/` add `src/` to `sys.path`, so an editable install is not strictly required for running them. For regular development, `pip install -e ".[all]"` is recommended.

There is no separate `pytest` suite. The examples are the smoke tests: if the relevant example runs and produces the expected figure or table, the feature works.

## Design principles

cfizz is built around four commitments.

### 1. End-to-end flow

Every analysis step should be reachable without hand-wiring intermediate files. The intended chain is:

```text
FASTQ (optional) → .mcool → Compartment / TAD / Loop → altered features → publication-ready figure
```

The FASTQ preprocessing scripts are an optional reference implementation. Users who already have `.mcool` files can start there.

### 2. Integration

cfizz treats Hi-C, BigWig, GTF and BED as first-class inputs. Tracks are composed through dictionaries:

```python
tracks = [
    {"file": "genes.gtf", "gtf_style": "flybase", "labels": True},
    {"file": "CTCF.bw", "color": "#264653", "name": "CTCF"},
    {"file": "enhancers.bed", "color": "purple", "name": "Enhancers"},
]
```

Hi-C panels use the same idea:

```python
hics = [
    {
        "file": "sample.mcool",
        "name": "sample",
        "balance": True,
        "resolution": 10_000,
        "triangle_ratio": 1.0,
        "insulation_path": "tad/1_0.sample.10000.insulation.tsv",
        "window_size": 100_000,
        "loops_path": "loops/sample.10k.loops.txt",
    }
]
```

### 3. Zone-Terminal visualization

The Zone-Terminal module is cfizz's visualization layer:

- **Zonal**: overlays adapt to Compartment, TAD and Loop regions.
- **Terminal**: sits at the end of the pipeline, turning results into publication-grade figures.

Figures expose physical sizing (`width_cm` / `height_cm` or `plot_size`) and high-dpi output. Heatmaps embedded in SVG use `svg_dpi` (default 1000) so they remain sharp when zoomed.

### 4. Example = test = documentation

Every public feature should have at least one runnable example in `examples/`. An example is:

- a **smoke test** — running it verifies the feature;
- **documentation** — it is the canonical usage pattern;
- a **demo** — the generated figure is the visual proof.

When adding or changing a function, add or update its example.

## Project layout

```text
cfizz/
├── src/cfizz/              # installable package (src layout)
│   ├── analyze/            # Compartment / TAD / Loop / O/E / distance
│   ├── api/integrated/     # quick_plot / saddle / APA / tracks
│   ├── io/                 # cooler / BED / GTF / BigWig / path helpers
│   ├── preprocessing/      # optional FASTQ → mcool wrappers
│   ├── utils/              # coordinates / intervals / BEDPE
│   ├── viz/                # heatmap / compartment / TAD / loop / pileup
│   └── config.py
├── examples/               # runnable examples (smoke tests + docs)
├── demo/                   # 15 MB demo dataset and outputs
├── assets/figures/         # images referenced by the README
├── README.md               # English main document
├── README.zh-CN.md         # Simplified Chinese version
├── CONTRIBUTING.md         # this file
├── CONTRIBUTING.zh-CN.md   # Chinese version
└── LICENSE
```

Do not put the package directly at `cfizz/`, `io/`, `analyze/` and so on. The package lives under `src/cfizz/`; the repository root only contains project files.

## Coding conventions

### Naming

- Functions: `verb_noun`, e.g. `compute_insulation`, `plot_multi_heatmap`.
- Classes: `PascalCase`.
- Constants: `UPPER_SNAKE_CASE`.
- Config-dict keys: `snake_case`.
- Internal helpers: leading underscore, e.g. `_helper_function`.

### Public API

- Public functions live at module level.
- Every public function has a docstring with a one-line summary, parameters, returns and a callable example where practical.
- Internal helpers are prefixed with `_`.
- Keep backwards compatibility when adding optional parameters; append new keyword arguments rather than changing positional order.

### Imports and typing

- Import order: standard library → third-party → cfizz.
- Avoid `import *`.
- Public APIs use type hints.
- Use `Optional[X]`, not `X | None`, for Python 3.9 compatibility.

### Figures and output

- Use `width_cm` / `plot_size` for physical sizing; do not hard-code pixel dimensions for publication figures.
- Save both PNG and SVG for user-facing figures.
- Preserve the distinction between vector text/lines and rasterized heatmaps.
- Prefer `svg_dpi=1000` or higher for heatmaps embedded in SVG; document any exception.

## Example conventions

Examples should:

- use relative paths and the bundled demo data where possible;
- be runnable from the repository root;
- print the output location;
- write to `demo/output/`;
- fail loudly on missing required inputs;
- avoid machine-specific absolute paths in committed examples (document overrides such as `CFIZZ_FASTA` instead).

When adding an example:

1. Add the script under `examples/` (or `examples/diff/`, `examples/integrated/`).
2. Use `demo/data/` and `demo/output/` unless the example explicitly documents external data.
3. Run it.
4. Inspect the generated figure or table.
5. Add a short section to `README.md` and `README.zh-CN.md`.
6. If the output is shown in the README, add the image to `assets/figures/` or reference a committed demo output.

## Documentation conventions

- `README.md` is the English main document; `README.zh-CN.md` is the Simplified Chinese version.
- `README.md` uses absolute GitHub raw image URLs so the same content renders on GitHub and PyPI.
- `README.zh-CN.md` uses relative paths for GitHub and local preview.
- Keep the two language versions structurally aligned.
- Use sentence-case headings and one H1 per document.
- Use fenced code blocks with a language tag (`bash`, `python`, `text`, `mermaid`).
- Use the project name `cfizz` in lowercase, including in commands and code.
- In prose, use `altered feature` for the 6_x concept; code and file names keep `differential` / `diff`.
- Update the table of contents when adding or renaming headings.

## Testing and verification

There is no `pytest` suite. Verify changes by running the relevant example(s):

```bash
python examples/5_1_primary_analysis_template.py
python examples/5_2_primary_analysis_visualization.py
python examples/diff/6_1_differential_compute.py
```

For a quick syntax check:

```bash
python -m py_compile path/to/changed_file.py
```

For documentation changes, check:

- image paths exist;
- links and anchors resolve;
- code blocks are balanced;
- numbers and file names match the current outputs.

## Commit and pull-request workflow

1. Fork and clone the repository.
2. Create a branch: `git checkout -b feat/your-feature`.
3. Make the change.
4. Run and inspect the relevant example(s).
5. Commit with a clear message.
6. Open a pull request describing what changed, why, and which example demonstrates it.

Commit messages: a short first line (≤72 chars), followed by a body explaining **why**.

Example:

```text
Add SVG rasterization dpi control for quick_plot_integrated

- Pass svg_dpi through to HeatmapTracks.plot
- Keep PNG dpi independent
- Verified with examples/5_2_primary_analysis_visualization.py
```

## Bug reports

Please include:

- cfizz version (`python -c "import cfizz; print(cfizz.__version__)"`);
- Python version and operating system;
- minimal reproduction code or the example you ran;
- expected vs actual output;
- traceback or relevant log excerpt.

Open issues at <https://github.com/Cuixiaojian21/cfizz/issues>.

## Community

- Issues: <https://github.com/Cuixiaojian21/cfizz/issues>
- Maintainer: cfizz developers

Be respectful and constructive in discussions.

## License

By contributing, you agree that your contributions are licensed under the MIT License. See [LICENSE](LICENSE) for details.
