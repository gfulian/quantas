# -*- coding: utf-8 -*-

"""Acoustic-wave propagation in anisotropic elastic solids."""

from .christoffel import ChristoffelSolver
from .isotropic import IsotropicSeismicVelocities, isotropic_seismic_velocities
from .medium import ElasticMedium
from .enhancement import (
    DirectionalEnhancementResult,
    EnhancementModeResult,
    cofactor_matrix,
)
from .field import (
    EnhancementFieldResult,
    EnhancementModeField,
    GroupFieldResult,
    GroupModeField,
    ModePairField,
    PhaseFieldResult,
    PhaseModeField,
    build_enhancement_field,
    build_group_field,
    build_phase_field,
)
from .group import DirectionalGroupResult, GroupModeResult
from .modes import (
    BRANCH_INDEX,
    BRANCH_ORDER,
    MODE_DESCRIPTIONS,
    MODE_INDEX,
    MODE_ORDER,
    MODE_PAIR_INDEX,
    MODE_PAIR_ORDER,
    MODE_PAIR_SYMBOLS,
    MODE_SYMBOLS,
    DirectionalPhaseResult,
    ModePair,
    PhaseModeResult,
    PolarizationBranch,
    WaveMode,
)
from .sampling import (
    SamplingLevel,
    SeismicFieldResult,
    sample_seismic_field,
)
from .tracking import (
    PolarizationBranchField,
    PolarizationTrackingResult,
    align_axial_vector,
    track_polarizations,
)

__all__ = [
    "BRANCH_INDEX",
    "BRANCH_ORDER",
    "MODE_DESCRIPTIONS",
    "MODE_INDEX",
    "MODE_ORDER",
    "MODE_PAIR_INDEX",
    "MODE_PAIR_ORDER",
    "MODE_PAIR_SYMBOLS",
    "MODE_SYMBOLS",
    "ChristoffelSolver",
    "ElasticMedium",
    "IsotropicSeismicVelocities",
    "DirectionalEnhancementResult",
    "DirectionalGroupResult",
    "EnhancementFieldResult",
    "EnhancementModeField",
    "EnhancementModeResult",
    "DirectionalPhaseResult",
    "GroupFieldResult",
    "GroupModeField",
    "GroupModeResult",
    "ModePair",
    "ModePairField",
    "PhaseFieldResult",
    "PhaseModeField",
    "PhaseModeResult",
    "PolarizationBranch",
    "PolarizationBranchField",
    "PolarizationTrackingResult",
    "WaveMode",
    "sample_seismic_field",
    "SeismicFieldResult",
    "SamplingLevel",
    "align_axial_vector",
    "build_enhancement_field",
    "build_group_field",
    "build_phase_field",
    "cofactor_matrix",
    "isotropic_seismic_velocities",
    "track_polarizations",
]
