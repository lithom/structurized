structurized
=================

Modular Java 23 Maven stack for cheminformatics, medchem endpoint standardization, and downstream analytics.

Project coordinates:
- Parent build: `tech.molecules:structurized:0.1.0-SNAPSHOT`
- Core module: `tech.molecules:structurized-core:0.1.0-SNAPSHOT`
- Endpoints module: `tech.molecules:structurized-endpoints:0.1.0-SNAPSHOT`
- Analytics module: `tech.molecules:structurized-analytics:0.1.0-SNAPSHOT`

Build
-----
- Requires Java 23 and Maven 3.9+
- Commands:
  - `mvn package` – build all modules
  - `mvn test` – run all module tests
  - `mvn -pl structurized-core test` – run core-module tests only

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
- Base package: `tech.molecules.structurized`
- Java release: 23
- Parent POM: [`pom.xml`](/home/lithom/dev_chem/structurized/pom.xml)
- Core module POM: [`structurized-core/pom.xml`](/home/lithom/dev_chem/structurized/structurized-core/pom.xml)
- Endpoints module POM: [`structurized-endpoints/pom.xml`](/home/lithom/dev_chem/structurized/structurized-endpoints/pom.xml)
- Analytics module POM: [`structurized-analytics/pom.xml`](/home/lithom/dev_chem/structurized/structurized-analytics/pom.xml)
- Existing cheminformatics implementation now lives in `structurized-core`
- Main pairwise engine: `tech.molecules.structurized.transforms.TransformationSplitter`
- Scaffold-mode entry point: `tech.molecules.structurized.scaffolds.ScaffoldAnalyzer`
- Scaffold discovery entry point: `tech.molecules.structurized.scaffolds.ScaffoldDiscoveryEngine`
- Internal Swing validation GUI: `tech.molecules.structurized.gui.ScaffoldDiscoverySwingApp`
- `structurized-endpoints` is reserved for the PRISM protocol and standardized endpoint abstractions
- `structurized-analytics` is reserved for analytics that combine structural methods with endpoints
- Parent-aware context shell spec: `docs/CONTEXT_SHELL_ENCODING.md`
- Scaffold-mode notes: `docs/SCAFFOLD_MODE.md`
- Scaffold discovery notes: `docs/SCAFFOLD_DISCOVERY.md`
- GUI test app notes: `docs/GUI_TEST_APP.md`
- Review and project notes: `docs/STRUCTURIZED_REVIEW.md`
- Verified OpenChemLib usage notes: `docs/OPENCHEMLIB_METHODS.md`
