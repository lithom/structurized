structurized
=================

Minimal Java 23 Maven library for pairwise, core-relative cheminformatics analysis based on OpenChemLib.

Project coordinates: `tech.molecules:structurized:0.1.0-SNAPSHOT`

Build
-----
- Requires Java 23 and Maven 3.9+
- Commands:
  - `mvn package` – build the JAR
  - `mvn test` – run tests

Usage
-----
Example of counting heavy atoms from a SMILES string using OpenChemLib:

```java
import tech.molecules.structurized.OpenChemLibUtil;

public class Demo {
  public static void main(String[] args) {
    int atoms = OpenChemLibUtil.atomCountFromSmiles("c1ccccc1");
    System.out.println(atoms); // 6
  }
}
```

Notes
-----
- Base package: `tech.molecules`
- Java release: 23
- OpenChemLib dependency: `com.actelion.research:openchemlib:${openchemlib.version}` (see `pom.xml`)
- Main pairwise engine: `tech.molecules.structurized.transforms.TransformationSplitter`
- Public pairwise model types: `PairTransformation`, `TransformationGroup`, `TransformationSignature`, `TransformationType`
- Pairwise bench/demo: `tech.molecules.structurized.transforms.TransformationBench` and `TransformationBenchDemo`
- Scaffold-mode entry point: `tech.molecules.structurized.scaffolds.ScaffoldAnalyzer`
- Scaffold discovery entry point: `tech.molecules.structurized.scaffolds.ScaffoldDiscoveryEngine`
- Internal Swing validation GUI: `tech.molecules.structurized.gui.ScaffoldDiscoverySwingApp`
- Context feature extraction: `tech.molecules.structurized.context.ContextFeatureExtractor`
- Signatures store both a compact canonical context shell and a larger `expandedRawContextIdcode` for downstream analysis
- Parent-aware context shell spec: `docs/CONTEXT_SHELL_ENCODING.md`
- Scaffold-mode notes: `docs/SCAFFOLD_MODE.md`
- Scaffold discovery notes: `docs/SCAFFOLD_DISCOVERY.md`
- GUI test app notes: `docs/GUI_TEST_APP.md`
- Review and project notes: `docs/STRUCTURIZED_REVIEW.md`
- Verified OpenChemLib usage notes: `docs/OPENCHEMLIB_METHODS.md`
