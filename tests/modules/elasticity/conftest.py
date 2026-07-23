"""Shared fixtures for elasticity-workflow tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.modules.elasticity.io.reader import ElasticityInputFileReader
from quantas.modules.elasticity.models import ElasticityInput


DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"


@pytest.fixture
def hydroxylapatite_input() -> ElasticityInput:
    """Return the reference hydroxylapatite elasticity input."""
    reader = ElasticityInputFileReader(DATA)
    return ElasticityInput(
        jobname=reader.jobname,
        stiffness=reader.stiffness,
        source=DATA,
    )


@pytest.fixture
def hydroxylapatite_tensor(
    hydroxylapatite_input: ElasticityInput,
) -> ElasticTensor:
    """Return the reference shared elastic tensor."""
    return ElasticTensor(hydroxylapatite_input.stiffness)
