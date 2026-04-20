# Structurized Concepts

## Purpose

`structurized-core` is a cheminformatics library for analyzing chemical structures relative to a
preserved core.

The central idea is simple:

- identify what is structurally preserved,
- identify what changed around that preserved part,
- represent those changes in a canonical and reusable way.

This supports two closely related analysis modes:

- pairwise transformation analysis between two compounds,
- scaffold-centric analysis of compounds relative to a predefined scaffold.

The same conceptual machinery underlies both.

## The core-relative view of molecular change

The project treats structural analysis as a comparison between:

- a preserved core,
- and one or more non-core regions attached to it.

In pairwise mode, the preserved core is the strict mapped common substructure between two
compounds.

In scaffold mode, the preserved core is a predefined scaffold template.

In both cases, everything outside the preserved core is treated as structural variation.

This gives a medchem-native decomposition:

- the core explains what remains the same,
- the non-core regions explain what was added, removed, or substituted,
- the attachment pattern explains where those changes act on the core.

## Pairwise transformation analysis

Pairwise analysis starts from two compounds `A` and `B`.

The implemented workflow is:

1. compute or provide a strict MCS mapping,
2. treat the mapped atoms as the preserved core,
3. identify all atoms outside the core in both molecules,
4. split those non-core atoms into connected components,
5. determine which core atoms serve as attachment points,
6. connect removed components, attachment atoms, and added components in one tripartite graph,
7. split that graph into connected components,
8. turn each connected component into one `TransformationGroup`,
9. canonicalize each group into a `TransformationSignature`.

This means that one compound pair is not forced into one single change object.
If two compounds differ at two independent positions, the code produces two independent local
transformation groups inside one overall `PairTransformation`.

The public pairwise model is:

- `PairTransformation`: the full `A -> B` relation,
- `TransformationGroup`: one independent local edit inside that relation,
- `TransformationType`: the coarse class of that local edit,
- `TransformationSignature`: the canonical representation of that local edit.

## Transformation groups

A `TransformationGroup` is the main local edit object.

It represents one connected structural change relative to the preserved core.
Its classification is intentionally high-level:

- `REPLACEMENT`
- `INSERTION`
- `DELETION`
- `MERGE`
- `SPLIT`
- `MULTI_CENTER`

These classes do not attempt to encode full reaction chemistry.
They summarize the topological shape of the change:

- how many removed regions exist,
- how many added regions exist,
- and whether they are connected through the same attachment network.

So a `TransformationGroup` is not “everything that changed from `A` to `B`”.
It is one independent local change extracted from the full pair relation.

## Canonical transformation signatures

A `TransformationSignature` is the canonical identity of one local transformation group.

Its job is to answer the question:

“When should two local edits be considered the same transformation?”

The current signature is built from:

- the removed fragment,
- the added fragment,
- the attachment pattern,
- a local context shell around the attachment points,
- signature metadata such as radius and reaction class,
- a stable hash-like identifier `sigId`.

Canonicalization matters because downstream analysis should not depend on:

- atom ordering in the input molecule,
- which compound pair first exposed a transformation,
- or incidental graph traversal details.

The signature is therefore intended as a reproducible transformation key for medchem series
analysis, transformation frequency counting, and later transformation-based analytics.

## Attachment context

The transformation itself is only part of the story.
In medchem, the local core environment around an attachment point also matters.

Structurized therefore stores a canonical attachment context in
`TransformationSignature.contextShellIdcode`.

Conceptually, this is:

- a copied shell of core atoms around the attachment atoms,
- canonically encoded,
- with explicit attachment-aware labels.

The stored shell is also parent-aware.
Before canonization, the code encodes parent-core properties such as:

- ring membership,
- aromaticity,
- local pi character,
- smallest parent ring size.

This is important because a copied shell fragment may truncate larger ring systems.
Without those parent-derived annotations, some chemically important context would be lost.

## Compact context versus expanded raw context

The library stores two related context representations:

- `contextShellIdcode`
- `expandedRawContextIdcode`

The compact context shell is the canonical signature context.
It is intentionally local and is meant to stay reusable across related transformations.

The expanded raw context is a larger neighborhood around the same attachment region.
It is not part of the canonical signature identity.
Its purpose is downstream feature derivation and exploratory analysis.

This separation reflects an important design choice:

- compact context for canonical grouping,
- larger context for richer later interpretation.

## Context features

The project separates raw stored context from derived context features.

Instead of baking every context abstraction into the signature object, the signature stores the
raw canonical shell and the context package derives features from it later.

Current feature families include:

- shell summary features,
- attachment atom features,
- attachment pair features,
- shell-topology-style summaries.

This keeps the signature itself stable while allowing downstream code to decide how much context
detail it wants to use.

## Scaffold mode

Scaffold mode generalizes the same core-relative logic from:

- compound versus compound

to:

- scaffold versus compound

Here the scaffold is the full preserved core.
Everything present in the compound outside the scaffold match is treated as an added region.

This produces a scaffold-centric decomposition of a compound into:

- the matched scaffold,
- occupied exit vectors,
- attached fragments,
- scaffold-relative substitution events.

The key public objects are:

- `ScaffoldTemplate`
- `ScaffoldMatch`
- `ScaffoldDecomposition`
- `ExitVector`
- `SubstitutionEvent`

The practical motivation is to represent compounds the way medicinal chemists usually think about
them:

- one core scaffold,
- one or more decorated positions,
- specific substituents attached at those positions.

## Exit vectors and symmetry

`ScaffoldTemplate` stores one candidate exit vector per scaffold atom and one symmetry class per
scaffold atom.

The template itself stays intentionally minimal.
It does not require manually curated exit-vector definitions.
Instead, observed decorated positions emerge from real scaffold-to-compound matches.

Symmetry matters because two scaffold atoms may be equivalent in the scaffold graph even if a
specific match maps a substituent to one concrete atom index.

The library therefore distinguishes between:

- the concrete matched scaffold atom,
- and the corresponding scaffold symmetry class.

This supports later aggregation across symmetry-equivalent positions while still retaining the
specific matched orientation when needed.

## Substitution events

A `SubstitutionEvent` is the scaffold-centric counterpart of a pairwise transformation group.

It represents one added region attached to the scaffold match.
The current classification is intentionally simple:

- `SINGLE_ATTACHMENT`
- `MULTI_ATTACHMENT`

This already separates ordinary substituents from multi-anchor events such as bridges,
annulations, or other scaffold-spanning additions.

Finer subclassification can be added later, but the current model already preserves the crucial
topological distinction.

## Scaffold discovery

Scaffold discovery is the dataset-level extension of scaffold mode.

The implemented V1 workflow is local and pragmatic:

1. compute OpenChemLib FFP512 fingerprints,
2. build small nearest-neighbor neighborhoods,
3. generate pairwise strict-MCS scaffold candidates inside those neighborhoods,
4. canonicalize and merge identical scaffold candidates,
5. expand local support,
6. compute simple scaffold statistics and scores.

The result is not a forced single scaffold per compound.
Instead, the library mines a useful library of scaffold candidates with associated support and
metadata.

This is closer to real medchem practice, where compounds often belong meaningfully to multiple
related scaffold views.

## What the library assumes

The current implementation assumes:

- normalized input molecules,
- chemically meaningful strict MCS mappings,
- OpenChemLib as the low-level graph and search engine,
- canonicalization as the basis for stable downstream identities.

Normalization is intentionally treated as an upstream responsibility for now.

## What the library does not try to be

Structurized is not currently:

- a full reaction informatics engine,
- a full scaffold ontology,
- a full normalization/standardization pipeline,
- or a complete medchem analytics platform.

Its current role is narrower and more valuable:

- extract core-relative structural changes,
- represent them canonically,
- analyze scaffold-relative substitution structure,
- provide a strong foundation for later medchem analytics.

## How the pieces fit together

The library now has one consistent conceptual arc:

1. pairwise mode extracts local transformations from compound pairs,
2. context encoding stabilizes local attachment environments,
3. scaffold mode reuses the same machinery for scaffold-relative decomposition,
4. scaffold discovery lifts this toward dataset-level scaffold mining.

That makes `structurized-core` best understood as a core-relative structural analysis engine.
The exact definition of the “core” can change, but the underlying decomposition logic stays the
same.

## Recommended reading order

For a new developer or user, the most useful reading order is:

1. this document,
2. `STRUCTURIZED_REVIEW.md`,
3. `CONTEXT_SHELL_ENCODING.md`,
4. `SCAFFOLD_MODE.md`,
5. `SCAFFOLD_DISCOVERY.md`,
6. `OPENCHEMLIB_METHODS.md`.

That path moves from overall concepts to implementation specifics.
