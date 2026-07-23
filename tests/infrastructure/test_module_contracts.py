"""Tests for uniform public contracts of scientific workflow modules."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.models import ModuleContract, ResultData, ResultMetadata
from quantas.modules.elasticity import MODULE_CONTRACT as ELASTICITY_CONTRACT
from quantas.modules.ha import MODULE_CONTRACT as HA_CONTRACT
from quantas.modules.qha import MODULE_CONTRACT as QHA_CONTRACT
from quantas.modules.thermoelasticity import (
    MODULE_CONTRACT as THERMOELASTICITY_CONTRACT,
)

CONTRACTS = (
    HA_CONTRACT,
    QHA_CONTRACT,
    ELASTICITY_CONTRACT,
    THERMOELASTICITY_CONTRACT,
)
PACKAGE_ROOT = Path(__file__).parents[2] / "src" / "quantas"


@pytest.mark.parametrize(
    ("contract", "name"),
    list(
        zip(
            CONTRACTS,
            ("ha", "qha", "elasticity", "thermoelasticity"),
            strict=True,
        )
    ),
)
def test_scientific_modules_expose_the_uniform_contract(
    contract: ModuleContract,
    name: str,
) -> None:
    """Each implemented module exposes the same frontend-neutral entry points."""
    assert contract.name == name
    assert contract.result_key == name
    assert callable(contract.read_input)
    assert callable(contract.run)
    assert callable(contract.read_hdf5)
    assert callable(contract.write_hdf5)
    assert callable(contract.build_report)
    assert callable(contract.build_plots)

    valid = ResultData(
        metadata=ResultMetadata(module=name, method="test"),
        results={name: object()},
    )
    contract.validate_result(valid)


@pytest.mark.parametrize("contract", CONTRACTS)
def test_module_contract_rejects_incompatible_results(
    contract: ModuleContract,
) -> None:
    """Contract validation catches wrong metadata and missing payloads."""
    with pytest.raises(ValueError, match="Expected"):
        contract.validate_result(
            ResultData(metadata=ResultMetadata(module="other", method="test"))
        )

    with pytest.raises(ValueError, match="does not contain"):
        contract.validate_result(
            ResultData(metadata=ResultMetadata(module=contract.name, method="test"))
        )


def test_scientific_module_packages_do_not_contain_cli_subpackages() -> None:
    """Click wiring is centralized under quantas.cli, outside workflows."""
    for module_name in ("ha", "qha", "elasticity", "thermoelasticity"):
        assert not (PACKAGE_ROOT / "modules" / module_name / "cli").exists()
