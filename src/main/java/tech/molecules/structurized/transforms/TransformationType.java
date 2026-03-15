package tech.molecules.structurized.transforms;

/**
 * High-level classification of one independent local transformation group.
 */
public enum TransformationType {
    REPLACEMENT,
    INSERTION,
    DELETION,
    MERGE,
    SPLIT,
    MULTI_CENTER
}
