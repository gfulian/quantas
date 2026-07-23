"""Scientific contracts for Mie--Grüneisen--Debye thermal pressure."""

from __future__ import annotations

import pytest

from quantas.core.physics.eos import (
    MGDNormalization,
    MGDParameters,
    MGDVariant,
    MGDVolumeBasis,
    ThermalPressureFamily,
    ThermalPressureModel,
    parse_thermal_pressure_model,
    thermal_pressure_model_contracts,
    thermal_pressure_parameter_names,
)


def test_thermal_pressure_registry_and_aliases() -> None:
    models = thermal_pressure_model_contracts()
    assert tuple(model.tag for model in models) == (
        "holland-powell-einstein",
        "mie-gruneisen-debye:full",
        "mie-gruneisen-debye:q-compromise",
    )
    assert parse_thermal_pressure_model("Einstein") == models[0]
    assert parse_thermal_pressure_model("Mie-Grüneisen-Debye") == models[1]
    assert parse_thermal_pressure_model("MGD:q-compromise") == models[2]
    assert parse_thermal_pressure_model("q compromise") == models[2]
    with pytest.raises(ValueError, match="unknown thermal-pressure model"):
        parse_thermal_pressure_model("unknown")


def test_family_variant_contract_is_explicit() -> None:
    einstein = ThermalPressureModel("holland-powell-einstein")
    assert einstein.family_name is ThermalPressureFamily.HOLLAND_POWELL_EINSTEIN
    assert einstein.mgd_variant is None
    assert not einstein.requires_mgd_normalization
    with pytest.raises(ValueError, match="has no variant"):
        ThermalPressureModel("einstein", "full")

    mgd = ThermalPressureModel("mgd")
    assert mgd.family_name is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
    assert mgd.mgd_variant is MGDVariant.FULL
    assert mgd.requires_mgd_normalization


def test_parameter_ownership_is_stable() -> None:
    assert thermal_pressure_parameter_names("einstein") == (
        "temperature_ref",
        "alpha_ref",
        "theta_e",
    )
    assert thermal_pressure_parameter_names("mgd") == (
        "temperature_ref",
        "theta_d0",
        "gamma0",
        "q",
    )
    assert thermal_pressure_parameter_names("mgd:q-compromise") == (
        "temperature_ref",
        "theta_d0",
        "gamma0",
    )


def test_cell_normalization_accepts_formula_and_z_or_atom_count() -> None:
    inferred = MGDNormalization.cell(formula="NaF", formula_units_per_cell=4)
    explicit = MGDNormalization.cell(atoms_per_cell=8)
    redundant = MGDNormalization.cell(
        atoms_per_cell=8,
        formula="NaF",
        formula_units_per_cell=4,
    )
    assert inferred == redundant
    assert inferred.basis is MGDVolumeBasis.CELL
    assert inferred.atoms_per_cell == pytest.approx(8.0)
    assert inferred.atoms_per_formula_unit == pytest.approx(2.0)
    assert inferred.canonical_volume_unit == "angstrom^3"
    assert explicit.atoms_per_cell == pytest.approx(8.0)

    with pytest.raises(ValueError, match="requires atoms_per_cell or formula"):
        MGDNormalization.cell(formula="NaF")
    with pytest.raises(ValueError, match="requires formula_units_per_cell"):
        MGDNormalization.cell(atoms_per_cell=8, formula="NaF")
    with pytest.raises(ValueError, match="requires a formula"):
        MGDNormalization("cell", 8, formula_units_per_cell=4)
    with pytest.raises(ValueError, match="inconsistent"):
        MGDNormalization.cell(
            atoms_per_cell=7,
            formula="NaF",
            formula_units_per_cell=4,
        )


def test_molar_normalization_uses_one_formula_unit() -> None:
    inferred = MGDNormalization.molar_formula_unit(formula="CaMg(CO3)2")
    explicit = MGDNormalization.molar_formula_unit(
        atoms_per_formula_unit=10,
        formula="CaMg(CO3)2",
    )
    assert inferred == explicit
    assert inferred.basis is MGDVolumeBasis.MOLAR_FORMULA_UNIT
    assert inferred.atoms_per_formula_unit == pytest.approx(10.0)
    assert inferred.atoms_per_cell is None
    assert inferred.canonical_volume_unit == "cm^3 mol^-1"
    with pytest.raises(ValueError, match="requires atoms_per_formula_unit"):
        MGDNormalization.molar_formula_unit()
    with pytest.raises(ValueError, match="only valid for cell-volume"):
        MGDNormalization("molar", 2, formula_units_per_cell=4)


def test_mgd_parameter_variants_do_not_hide_q_compromise() -> None:
    full = MGDParameters(
        temperature_ref=295.0,
        theta_d0=459.0,
        gamma0=1.547,
        q=0.94,
    )
    full.validate_for("mgd")
    with pytest.raises(ValueError, match="must omit q"):
        full.validate_for("mgd:q-compromise")

    compromise = MGDParameters(
        temperature_ref=295.0,
        theta_d0=459.0,
        gamma0=1.547,
    )
    compromise.validate_for("mgd:q-compromise")
    with pytest.raises(ValueError, match="require q"):
        compromise.validate_for("mgd")
    with pytest.raises(ValueError, match="require a Mie-Gruneisen-Debye"):
        compromise.validate_for("einstein")


def test_mgd_parameter_domains_are_validated() -> None:
    with pytest.raises(ValueError, match="temperature_ref cannot be negative"):
        MGDParameters(-1.0, 459.0, 1.5, 1.0)
    with pytest.raises(ValueError, match="theta_d0 must be positive"):
        MGDParameters(295.0, 0.0, 1.5, 1.0)
    negative_gamma = MGDParameters(295.0, 459.0, -0.5, 1.0)
    assert negative_gamma.gamma0 == pytest.approx(-0.5)
    with pytest.raises(ValueError, match="gamma0 must be finite"):
        MGDParameters(295.0, 459.0, float("nan"), 1.0)
    with pytest.raises(ValueError, match="q must be finite"):
        MGDParameters(295.0, 459.0, 1.5, float("nan"))
