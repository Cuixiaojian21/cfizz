"""
cfizz.preprocessing - upstream Hi-C preprocessing reference scripts.

This subpackage bundles three reference scripts (two bash, one Python)
that wrap the standard open2C Hi-C preprocessing toolchain (fastp,
bwa-mem2, pairtools, cooler) into a parameterized, skip-if-exists,
well-logged command-line interface.

It is **not a Python API**: the scripts are meant to be invoked from
the shell, not imported. They are packaged with the cfizz distribution
so that the upstream pipeline is discoverable alongside the downstream
analysis.

See ``cfizz/preprocessing/README.md`` for the full usage guide.
"""
