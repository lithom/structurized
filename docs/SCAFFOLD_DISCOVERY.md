# Scaffold Discovery

## Scope of V1

The current discovery engine mines a useful library of scaffold candidates from a set of normalized
compounds. It does not try to assign one unique scaffold per compound. Compounds may support many
scaffold candidates.

The implementation is intentionally conservative:

- fragment fingerprints are used only for neighborhood generation
- scaffold candidates are generated only from strict pairwise common cores
- support expansion is local to the neighborhoods that discovered a scaffold
- all expensive chemistry operations are cached

## Main classes

- `ScaffoldDiscoveryEngine`
- `ScaffoldDiscoveryConfig`
- `ScaffoldDiscoveryResult`
- `ScaffoldCandidate`
- `CompoundRecord`

## Discovery algorithm

### 1. Fingerprints and neighborhoods

For every input compound, the engine computes an OpenChemLib `FFP512` fingerprint using
`DescriptorHandlerLongFFP512`.

Then it computes a top-k nearest-neighbor list by Tanimoto similarity.

These neighborhoods are used only to decide which local pairs are worth exploring for scaffold
generation. They are not used as final scaffold evidence.

### 2. Seed neighborhoods

The engine selects seed compounds, currently either:

- all compounds
- or the first / shuffled subset if `maxSeeds` is configured

For every seed `s`, the local compound set is:

- `{s} U neighbors(s)`

### 3. Pairwise scaffold generation

For every seed-neighbor pair:

- compute strict pairwise MCS using the provided strict MCS provider
- extract the mapped core from one side
- require a minimum scaffold heavy-atom count
- optionally require the candidate scaffold to be connected
- canonicalize the extracted core into a `ScaffoldTemplate`

All identical scaffold templates are merged by canonical scaffold IDCode.

### 4. Local support expansion

Every merged scaffold candidate stores the union of local neighborhood compounds from the discovery
events that produced it.

The engine then matches each scaffold candidate back only against this local support set. It does not
search the entire dataset in V1.

For every successful scaffold match:

- the compound is registered as a supporting compound
- the scaffold-relative decomposition is computed
- observed exit vectors are collected
- explained fraction is accumulated

### 5. Score components

For each scaffold candidate, the engine currently stores:

- `supportCount`
- `averageExplainedFraction`
- `scaffoldHeavyAtomCount`
- `observedExitVectorCount`
- `combinedScore`

The default `combinedScore` is intentionally simple and support-dominated:

- `support * 1000 + explainedFraction * 100 + scaffoldSize * 10 + exitVectorCount`

The raw score components remain available for later alternative ranking strategies.

## Caching

The V1 engine caches:

- pairwise scaffold candidates generated from strict MCS results
- scaffold-to-compound decomposition results during support expansion

This is necessary because overlapping seed neighborhoods would otherwise repeat the same expensive
operations many times.

## Current limitations

1. Discovery uses only pairwise common cores.
   Multi-compound common scaffold mining is not implemented yet.

2. Support expansion is local.
   Promising scaffolds are not yet expanded against the full dataset.

3. Candidate scoring is heuristic.
   The current default combined score is intentionally simple.

4. No scaffold hierarchy is built yet.
   Nested scaffolds are all kept as independent candidates.

5. No scaffold discovery GUI metadata is stored yet beyond the current result objects.

## Recommended next steps

1. Add approximate 3-compound scaffold generation from promising local pairwise candidates.
2. Add optional global support expansion for high-scoring local candidates.
3. Add scaffold containment / hierarchy relationships.
4. Add dataset-level exit-vector occupancy and substitution summaries.
