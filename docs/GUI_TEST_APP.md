# GUI Test App

## Purpose

This is a minimal internal Swing GUI for validating scaffold discovery on real datasets.

It is not intended to be the final user-facing chemistry GUI. Its purpose is:

- load a simple SMILES file
- run scaffold discovery
- inspect candidate scaffolds and their score components

## Entry point

Main class:

- `tech.molecules.structurized.gui.ScaffoldDiscoverySwingApp`

## Supported input format

The GUI expects a simple text file with one compound per line.

Accepted conventions:

- first token on each line is treated as the SMILES string
- additional columns are ignored
- empty lines are ignored
- lines starting with `#` are ignored
- a first header token `smiles` is ignored

Examples:

```text
c1ccccc1C toluene
c1ccccc1F fluorobenzene
COc1ccccc1 anisole
```

## What the GUI shows

For every discovered scaffold candidate, the table shows:

- rank
- combined score
- support count
- average explained fraction
- scaffold heavy atom count
- observed exit vector count
- scaffold SMILES

Selecting a row shows:

- scaffold depiction
- scaffold SMILES
- scaffold IDCode
- detailed score values
- discovery/support compound indices
- observed exit-vector atoms and symmetry classes

## Current limitations

1. This is a local validation tool.
   It is functional, but intentionally minimal.

2. There is no progress granularity yet.
   The UI only shows a busy indicator while discovery runs.

3. The table is candidate-centric only.
   There is no separate support-compound browser yet.

4. Discovery parameters are still minimal.
   The GUI currently exposes neighbor count, minimum scaffold size, and minimum support.
