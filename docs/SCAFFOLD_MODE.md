# Scaffold Mode

## Current scope

The first scaffold-mode implementation treats a predefined scaffold as the preserved core and
reuses the existing core-relative transformation splitter to decompose a compound into
scaffold-relative substitution events.

Implemented objects:

- `ScaffoldTemplate`
- `ExitVector`
- `ScaffoldMatch`
- `ScaffoldDecomposition`
- `SubstitutionEvent`
- `ScaffoldAnalyzer`

## Main idea

Instead of comparing compound `A` against compound `B`, scaffold mode compares:

- `scaffold -> compound`

The scaffold atoms are treated as the full preserved core. Everything present in the compound
outside the scaffold match is interpreted as an added region.

This means the existing `TransformationSplitter` can already do most of the heavy lifting.

## What a `ScaffoldTemplate` stores

The template is intentionally minimal:

- the unsubstituted scaffold graph
- a canonical scaffold IDCode
- one symmetry class per scaffold atom
- one candidate `ExitVector` per scaffold atom

No explicit manually defined exit vectors are required.
Observed occupied positions emerge from real scaffold-to-compound matches.

## Matching strategy

`ScaffoldAnalyzer` uses OpenChemLib `SSSearcher` with:

- scaffold copied as a fragment query
- unique-match mode

If multiple unique matches exist, the analyzer currently evaluates every match, runs the splitter
for each one, converts the result into substitution events, and then picks the lexicographically
smallest event signature set. This gives a deterministic first-pass symmetry handling strategy.

## Decomposition result

A successful `ScaffoldDecomposition` contains:

- the selected `ScaffoldMatch`
- the raw `TransformationGroup`s from scaffold-to-compound splitting
- scaffold-centric `SubstitutionEvent`s

Each `SubstitutionEvent` contains:

- the underlying `TransformationGroup`
- the concrete scaffold atom indices
- the corresponding scaffold symmetry classes
- the `ExitVector`s
- the added fragment IDCode
- a scaffold-centric event type

## Current substitution event typing

The current event typing is intentionally simple:

- `SINGLE_ATTACHMENT`
- `MULTI_ATTACHMENT`

This is enough to distinguish ordinary substituents from bridge / annulation / multi-anchor cases.
Finer distinctions such as `BRIDGE` versus `ANNULATION` are not implemented yet.

## Current limitations

1. This is scaffold matching, not scaffold discovery.
   The scaffold must already be provided.

2. Exit vectors are implicit.
   All scaffold atoms are candidate exit vectors; no chemical filtering of “reasonable”
   substitution sites is applied yet.

3. Match selection is deterministic but still heuristic.
   The analyzer chooses the smallest decomposition key among unique substructure matches.

4. Multi-attachment events are grouped, but not yet subclassified into bridge / cyclization /
   fused-ring extension categories.

5. No dataset-level scaffold aggregation exists yet.
   The current implementation is one scaffold against one compound at a time.

## Recommended next steps

1. Add dataset-level grouping:
   scaffold template + many compounds + per-exit-vector statistics.

2. Add more specific multi-anchor event classes:
   bridge, annulation, cyclization, multi-anchor extension.

3. Add richer symmetry-aware aggregation:
   concrete scaffold atom versus symmetry class versus matched-orientation information.

4. Decide whether to introduce optional chemical filtering of candidate exit vectors.

5. Later, add scaffold discovery separately from scaffold decomposition.
