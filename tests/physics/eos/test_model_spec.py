from __future__ import annotations

import pytest

from quantas.core.physics.eos import EOSFamily, EOSModel, parse_eos_model


@pytest.mark.parametrize(
    ("value", "tag", "family", "order"),
    [
        ("M", "M", EOSFamily.MURNAGHAN, None),
        ("BM", "BM3", EOSFamily.BIRCH_MURNAGHAN, 3),
        ("BM2", "BM2", EOSFamily.BIRCH_MURNAGHAN, 2),
        ("birch-murnaghan4", "BM4", EOSFamily.BIRCH_MURNAGHAN, 4),
        ("PT", "PT3", EOSFamily.NATURAL_STRAIN, 3),
        ("NS2", "PT2", EOSFamily.NATURAL_STRAIN, 2),
        ("V", "V3", EOSFamily.VINET, 3),
        ("V2", "V2", EOSFamily.VINET, 2),
        ("T4", "T4", EOSFamily.TAIT, 4),
    ],
)
def test_parse_eos_model_resolves_family_order_and_tag(value, tag, family, order):
    model = parse_eos_model(value)

    assert model.tag == tag
    assert model.family is family
    assert model.order == order


def test_model_reports_order_dependent_free_parameters():
    assert parse_eos_model("BM2").energy_parameter_names == ("E0", "K0", "V0")
    assert parse_eos_model("BM3").energy_parameter_names == ("E0", "K0", "KP", "V0")
    assert parse_eos_model("BM4").energy_parameter_names == (
        "E0",
        "K0",
        "KP",
        "KPP",
        "V0",
    )
    assert parse_eos_model("BM2").pressure_parameter_names == ("K0", "V0")
    assert parse_eos_model("BM4").pressure_parameter_names == (
        "K0",
        "KP",
        "KPP",
        "V0",
    )


def test_tait_is_pressure_only_and_murnaghan_has_no_order():
    assert not parse_eos_model("T3").supports_energy
    with pytest.raises(ValueError, match="does not define an EOS order"):
        EOSModel(EOSFamily.MURNAGHAN, 3)


def test_parser_rejects_unknown_unsupported_and_conflicting_orders():
    with pytest.raises(ValueError, match="unknown equation of state"):
        parse_eos_model("unknown")
    with pytest.raises(ValueError, match="unsupported order"):
        parse_eos_model("V4")
    with pytest.raises(ValueError, match="conflicting EOS orders"):
        parse_eos_model("BM2", order=3)


def test_available_model_registry_separates_pressure_and_energy_eos():
    from quantas.core.physics.eos import available_eos_models, available_eos_tags

    pressure_tags = tuple(model.tag for model in available_eos_models())
    energy_tags = tuple(
        model.tag for model in available_eos_models(require_energy=True)
    )

    assert pressure_tags == (
        "M",
        "BM2",
        "BM3",
        "BM4",
        "PT2",
        "PT3",
        "PT4",
        "V2",
        "V3",
        "T2",
        "T3",
        "T4",
    )
    assert energy_tags == pressure_tags[:-3]
    assert available_eos_tags(require_energy=True, include_default_aliases=True) == (
        "M",
        "BM",
        "BM2",
        "BM3",
        "BM4",
        "PT",
        "PT2",
        "PT3",
        "PT4",
        "V",
        "V2",
        "V3",
    )
