# OpenChemLib methods used in structurized

This note documents the OpenChemLib APIs currently used by `structurized`, based on the local source tree in `/home/lithom/dev/openchemlib`.

## `SmilesParser`

Used in:

- `OpenChemLibUtil.atomCountFromSmiles()`
- `TransformationBenchDemo.parseSmiles()`

Verified behavior:

- `new SmilesParser()` means SMILES mode, not SMARTS mode. The constructor comment states that SMARTS features are rejected and created molecules are not marked as fragments.
- `parse(StereoMolecule mol, String smiles)` clears and fills the passed molecule.

Relevant local source:

- `SmilesParser()` constructor: `openchemlib/.../SmilesParser.java:103-134`
- `parse(StereoMolecule, String)`: `openchemlib/.../SmilesParser.java:336-353`

Assessment:

- Current usage is correct.
- Calling `ensureHelperArrays(...)` after parsing is appropriate when downstream code depends on neighbor/ring/aromaticity information.

## `ensureHelperArrays(...)`

Used throughout the project before ring, aromaticity, bond, or stereo-sensitive operations.

Verified behavior:

- `ensureHelperArrays(int required)` computes neighbor, ring/aromaticity, parity, or CIP helper data depending on the requested level.
- `cHelperNeighbours` provides connectivity and `mPi`.
- `cHelperRings` adds ring and aromaticity information.

Relevant local source:

- `ExtendedMolecule.ensureHelperArrays(...)`: `openchemlib/.../ExtendedMolecule.java:3638-3665`

Assessment:

- Current usage is broadly correct.
- For substructure searching, OpenChemLib already calls `ensureHelperArrays(Molecule.cHelperNeighbours)` internally in `SSSearcher.setMolecule()` and `setFragment()`.
- For aromaticity-sensitive attachment and bond logic, using `cHelperRings` in `structurized` is appropriate.

## `MCSFast`

Used in:

- `TransformationBenchDemo.OCLMCSFastProvider`

Verified behavior:

- `new MCSFast()` defaults to `PAR_CLEAVE_RINGS`.
- `set(mol, frag)` stores the two molecules and initializes the search.
- `getMCS()` returns the largest common substructure or `null`.
- The returned substructure is explicitly marked as a fragment via `fragSubBonds.setFragment(true)`.

Relevant local source:

- constructor and ring-mode setup: `openchemlib/.../mcs/MCSFast.java:91-129`
- `getMCS()`: `openchemlib/.../mcs/MCSFast.java:392-410`
- returned MCS marked as fragment: `openchemlib/.../mcs/MCSFast.java:431-433`

Assessment:

- Current usage of `set(A, B)` and `getMCS()` is valid.
- The provider should still validate the returned mapping strictly, because `MCSFast` itself is a common-subgraph generator, not a full medchem-normalization or strict-identity policy layer.

## `SSSearcher`

Used in:

- `TransformationBenchDemo.OCLMCSFastProvider.firstMappingAtoB(...)`

Verified behavior:

- `setMol(fragment, molecule)` delegates to `setMolecule()` and `setFragment()`.
- `setFragment()` rejects inputs that are not marked as fragments: `!fragment.isFragment()` causes `mFragment = null`.
- `findFragmentInMolecule()` defaults to `cCountModeOverlapping`.
- `findFragmentInMolecule(int countMode, int matchMode)` lets the caller request first match only.
- `getMatchList()` returns fragment-to-molecule atom mappings.

Relevant local source:

- `setMol(...)`: `openchemlib/.../SSSearcher.java:155-164`
- `setFragment(...)` and fragment requirement: `openchemlib/.../SSSearcher.java:190-203`
- `findFragmentInMolecule(...)` overloads: `openchemlib/.../SSSearcher.java:459-497`
- `getMatchList()`: `openchemlib/.../SSSearcher.java:390-395`

Assessment:

- The important requirement is that the query must already be a fragment.
- Current usage is correct because the MCS returned by `MCSFast.getMCS()` is a fragment, and the code also calls `mcs.setFragment(true)` defensively.
- Using `cCountModeFirstMatch` is better than taking the first element after an unrestricted search when only one mapping is needed.

## `Canonizer`

Used in:

- `TransformationSplitter.extractFragmentWithAttachments(...)`
- `TransformationSplitter.canonizeInducedSubgraph(...)`

Verified behavior:

- `Canonizer.ENCODE_ATOM_CUSTOM_LABELS` is required if custom atom labels should affect ranking and be encoded into the IDCode.
- Without that flag, custom labels are not part of the canonical IDCode.

Relevant local source:

- `ENCODE_ATOM_CUSTOM_LABELS`: `openchemlib/.../Canonizer.java:64-65`
- custom-label encoding comment: `openchemlib/.../Canonizer.java:103-109`

Assessment:

- This was previously used incorrectly in `structurized` for attachment-labeled fragments.
- It is now fixed for fragment signatures.
- The context shell is now also canonized with attachment labels, so the decoded shell preserves ordered attachment atoms.

## `copyMoleculeByAtoms(...)`

Used in:

- `TransformationSplitter.extractFragmentWithAttachments(...)`

Verified behavior:

- `copyMoleculeByAtoms(destMol, includeAtom, recognizeDelocalizedBonds, atomMap)` copies selected atoms into a destination molecule.
- When the source molecule is a fragment, the destination molecule is also marked as a fragment.
- With `recognizeDelocalizedBonds=true`, OpenChemLib prepares ring helpers first.

Relevant local source:

- `ExtendedMolecule.copyMoleculeByAtoms(...)`: `openchemlib/.../ExtendedMolecule.java:100-119`

Assessment:

- Current usage is correct for extracting removed and added subgraphs before adding explicit attachment dummies.

## `IDCodeParser`

Used in:

- `ContextFeatureExtractor`

Verified behavior:

- `IDCodeParser` restores atom custom labels that were encoded with `Canonizer.ENCODE_ATOM_CUSTOM_LABELS`.
- This is what allows attachment-labeled context shells to be decoded later and analyzed externally.

Relevant local source:

- parser class: `openchemlib/.../IDCodeParser.java`
- custom-label decoding: `openchemlib/.../IDCodeParserWithoutCoordinateInvention.java:586-595`

Assessment:

- Current usage is correct for round-tripping attachment-aware context shells.

## Practical caveats for future scaffold work

1. First-match mappings are symmetry-arbitrary
   `SSSearcher` gives a valid match, but not a scaffold symmetry class.
   For scaffold-centric medchem analysis you will want explicit symmetry classes on scaffold atoms.

2. Strict identity is still your responsibility
   `MCSFast` plus `SSSearcher` gives you a workable graph mapping.
   Deciding whether that mapping is chemically strict enough for your project remains an application-level policy.

3. Fragment flags matter
   In OpenChemLib, fragment search is not just about graph content.
   The `isFragment` flag directly affects whether `SSSearcher` accepts the query.
