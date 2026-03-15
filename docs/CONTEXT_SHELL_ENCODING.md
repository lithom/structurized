# Context Shell Encoding

## Purpose

`structurized` stores the preserved-core attachment context of a transformation as
`TransformationSignature.contextShellIdcode`.

This object is meant to be the canonical raw context representation for downstream algorithms.
It must therefore preserve:

- the copied shell graph around the attachment atoms,
- the identity and ordering of the attachment atoms,
- chemically relevant parent-core semantics that may be lost when the shell truncates a larger ring system.

## What is stored

The stored `contextShellIdcode` is an OpenChemLib IDCode of a copied core-only shell fragment.
In addition, every signature now stores `expandedRawContextIdcode`, which uses the same encoding but a
larger fixed radius (`3` by default).

The shell is built by:

1. taking the strict preserved core on molecule `A`,
2. selecting all core atoms within radius `r` of the ordered attachment atoms,
3. copying that induced subgraph with `copyMoleculeByAtoms(...)`,
4. assigning custom atom labels to every copied shell atom,
5. canonizing with `Canonizer.ENCODE_ATOM_CUSTOM_LABELS`.

The shell graph itself is therefore a normal OpenChemLib fragment, but the labels make it
attachment-aware and parent-aware.

The two stored context objects have different roles:

- `contextShellIdcode`: compact canonical context used in the transformation signature hash
- `expandedRawContextIdcode`: larger auxiliary raw context for downstream feature derivation and
  exploratory analysis

## Why custom labels are needed

OpenChemLib correctly copies the selected subgraph, but helper-array properties on the copied shell
are evaluated on that copied shell graph only.

That means:

- if the shell still contains a complete ring, shell-local ring/aromatic information remains valid,
- if the shell cuts through a larger ring system, shell-local `isRingAtom()` / `isAromaticAtom()`
  may no longer match the original preserved core.

To avoid losing this information, the splitter writes parent-core annotations into the shell atom
labels before canonization.

## Per-atom annotation payload

Every shell atom receives a custom label encoded by `ContextLabelCodec`.

The payload currently contains:

- `ATT`: ordered attachment index, if the atom is an attachment atom
- `PR`: whether the atom was a ring atom in the full copied core
- `PA`: whether the atom was aromatic in the full copied core
- `PP`: parent-core `atomPi`
- `PS`: smallest ring size in the full copied core
- `PAIR`: pairwise relation data from this attachment to later attachment indices

`PAIR` entries contain:

- other attachment index
- same ring
- same aromatic system
- same fused ring system
- topological distance on the full copied core

Only the lower-index attachment atom stores a given pair relation. This keeps labels shorter while
remaining sufficient for downstream reconstruction.

## Source of truth for parent annotations

The parent annotations do not come from the truncated shell.

Instead, `TransformationSplitter` first copies the entire preserved core into an internal
`CoreGraphData` object and computes:

- ring membership,
- fused ring-system membership,
- aromatic ring-system membership,
- parent-core shortest-path distances.

The later shell annotations are derived from this full copied core graph.

## Downstream feature extraction

`ContextFeatureExtractor` decodes the `contextShellIdcode` and reconstructs:

- shell summary features,
- per-attachment atom features,
- pairwise attachment features,
- attachment distance matrix,
- min/max attachment distance.

When parent-aware annotations are present, the extractor prefers them over shell-local helper-array
values for:

- aromaticity,
- ring membership,
- `atomPi`,
- smallest ring size,
- pairwise attachment relations.

If no annotations are present, the extractor falls back to shell-local calculations.

## Current design intent

The current design intentionally separates:

- raw canonical context storage: `contextShellIdcode`
- larger auxiliary raw context storage: `expandedRawContextIdcode`
- derived algorithm-specific interpretation: `ContextFeatureExtractor`

This keeps `TransformationSignature` stable and lets surrounding algorithms choose which context
features they want to use.

The expanded raw context exists for cases where a broader neighborhood is useful for downstream
analysis, but where including that broader neighborhood directly in the canonical signature would make
transformation identities too specific and less reusable.

## Current limitation

The shell is still structurally a truncated fragment.

The parent-aware annotations preserve important chemistry semantics, but they do not turn the shell
back into the full original core. Algorithms that need information beyond the encoded payload still
need access to the original pair or a richer scaffold model.
