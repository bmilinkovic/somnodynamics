# Emergence Engine multiscale snapshot

The base snapshot of the standalone `multiscale/` package comes from
[Emergence Engine](https://github.com/bmilinkovic/emergence-engine), commit
`c9e9c7e16286f4b238a5dd8dc6c25703a895b216`. The `objective.py`,
`optimise.py`, and `parallel.py` files were updated from commit `adc6b05`
to distribute optimisation starts across workers and speed up the exact
proxy matrix calculation.

The copy includes the package licence and source attribution. SomnoDynamics
uses this local copy so a transfer of this repository to the cluster carries
the numerical engine with it. Its import name is `multiscale`.

To refresh it from a local Emergence Engine checkout, copy the contents of
`emergence-engine/multiscale/` into this directory, excluding `build/`,
`*.egg-info/`, caches, and other generated files, then update the commit
above. Re-run the DK68 synthetic smoke test and
inspect the numerical changes before using a new snapshot for an analysis.
