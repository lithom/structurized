package tech.molecules.structurized.analytics.mmp;

import com.actelion.research.chem.StereoMolecule;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

/**
 * Resolves PRISM subject identifiers to chemical structures.
 */
public interface StructureProvider {
    Optional<StereoMolecule> findStructure(String subjectId);

    default Map<String, StereoMolecule> fetchStructures(Collection<String> subjectIds) {
        Objects.requireNonNull(subjectIds, "subjectIds");
        LinkedHashMap<String, StereoMolecule> structures = new LinkedHashMap<>();
        for (String subjectId : subjectIds) {
            findStructure(subjectId).ifPresent(molecule -> structures.put(subjectId, new StereoMolecule(molecule)));
        }
        return Map.copyOf(structures);
    }
}
