# Contributing to cfizz

> Chromosome Flow Integration with Zone-Terminal Visualization
> End-to-end Hi-C differential analysis and zone-terminal visualization.

Thank you for your interest in cfizz! This guide covers the essentials for contributing code, bug fixes, documentation, and examples to the project.

## Quick Start

```bash
# Clone and install for development
git clone https://github.com/<username>/cfizz.git
cd cfizz
pip install -e ".[all]"

# Run the test suite (examples = tests, see test philosophy below)
python -m pytest cfizz/ -x  # optional; we primarily test via examples
```

## Design Philosophy

cfizz is built around four first-class commitments:

### 1. End-to-end pipeline (Flow)

Every analysis step should be reachable from a one-liner. From `.mcool` → TAD/Loop/A-B compartment call → differential comparison → publication-ready figure, the chain should not require hand-wiring.

### 2. Multi-omics integration (Integration)

cfizz treats Hi-C + BigWig + GTF + BED as first-class inputs. Tracks compose via a unified config dictionary:

```python
tracks = [
    {"type": "gtf",       "path": "genes.gtf"},
    {"type": "bigwig",    "path": "CTCF.bw"},
    {"type": "loop",      "path": "loops.bedpe"},
]
```

### 3. Zone-Terminal visualization (Zonal Visualization)

The Zone-Terminal module is cfizz's signature: a unified visualization layer that takes any genomic region × any combination of tracks × any number of samples, and produces publication-ready figures.

**Two semantic dimensions**:
- **Zonal (分区式适配)**: Adapts to A/B compartments, TADs, loops as zonal overlays.
- **Terminal (终端模块)**: Sits at the terminal step of the analysis pipeline, as the last-mile companion that converts analytical results into publication figures.

### 4. Example = Test = Documentation

Every public function ships with at least one runnable example in `examples/`. The example serves as:

- **Test**: If `python3 examples/...py` runs and produces a figure, the function works.
- **Documentation**: The example is the canonical usage pattern.
- **Demo**: The resulting figure is the visual proof.

When you add or change a function, **add or update its example**.

## Coding Conventions

### Naming
- Functions: `verb_noun` (e.g. `compute_insulation`, `plot_multi_heatmap`)
- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Config dict keys: `snake_case`
- Module-private: leading underscore (`_helper_function`)

### Public API surface
- Public functions live at the module level (no deep nesting in public APIs).
- Add a docstring with: one-line summary, Parameters, Returns, Example (callable).
- Internal helpers (used only inside cfizz) are prefixed with `_`.

### Imports
- Import order: stdlib → third-party → cfizz.
- Each import on its own line in production code.
- Avoid `import *`.
- cfizz must not depend on out-of-tree code; all functionality is implemented locally.

### Type hints
- All public APIs use type hints.
- Internal helpers may omit.
- Use `Optional[X]` not `X | None` for compatibility with Python 3.9.

## Submitting Changes

### Pull request workflow

1. Fork and clone the repository.
2. Create a feature branch: `git checkout -b feat/your-feature`.
3. Make your changes.
4. Run and verify the relevant example(s):
   ```bash
   python3 examples/<your-area>/<your-example>.py
   ```
   Verify the output figure looks correct.
5. Commit your work with a descriptive message.
6. Open a Pull Request describing: what changed, why, and which example(s) demonstrate the change.

### Commit messages

Use a clear, descriptive summary on the first line (≤72 chars). Add a body explaining **why** rather than **what**.

Good:
```
Add quick_plot_integrated support for insulation overlay

- Resolves #123
- Enables TAD boundary display next to the heatmap
- Tested via examples/integrated/7_3_multi_gene.py
```

### Bug reports

When filing a bug, please include:
- cfizz version (`python3 -c "import cfizz; print(cfizz.__version__)"`)
- Python version
- Operating system
- Minimal reproducing code (or the example you're running)
- Expected vs. actual output
- Traceback (if applicable)

## Project Structure

```
cfizz/
├── __init__.py
├── io/             # I/O: coolers, BED/GTF/BigWig readers
├── analyze/        # Algorithms: insulation, compartments, TADs, loops, O/E
├── viz/            # Visualization: heatmaps, tracks, pileups, layouts
├── api/            # High-level user-facing API
│   └── integrated/ # Zone-Terminal visualization
├── utils/          # Cross-cutting utilities (coordinates, range)
├── examples/       # Example = Test = Documentation
│   ├── heatmap/
│   ├── compartment/
│   ├── tad/
│   ├── loop/
│   ├── distance/
│   ├── diff/
│   ├── integrated/
│   └── io/
└── ...
```

## Community

- Issues: <github-issues-url>
- Discussions: <github-discussions-url>
- Maintainer: cfizz developers

## License

By contributing, you agree that your contributions will be licensed under the MIT License. See [LICENSE](LICENSE) for details.
