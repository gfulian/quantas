from __future__ import annotations

import numpy as np

from quantas.core.events import Event, EventLevel
from quantas.cli.ha_observer import HATextObserver
from quantas.cli.ha_render import render_table
from quantas.modules.ha.models import HAInput, HAOptions, HAResult
from quantas.modules.ha.report import options_table


def test_render_table_formats_neutral_table():
    text = render_table(
        options_table(HAOptions(temperature_min=0.0, temperature_max=10.0))
    )

    assert "Selected options" in text
    assert "Temperature minimum" in text
    assert "Energy unit" in text


def test_text_observer_collects_structured_events(tmp_path):
    report = tmp_path / "ha.txt"
    observer = HATextObserver(report_file=report, silent=True, show_progress=True)

    input_data = HAInput(
        jobname="MgO",
        natoms=2,
        qpoints=1,
        volume=np.array([10.0]),
        energy=np.array([-1.0]),
        frequencies=np.ones((1, 6, 1)),
        weights=np.array([1.0]),
    )
    result = HAResult(
        jobname="MgO",
        volume=np.array([10.0]),
        temperature=np.array([300.0]),
        static_energy=np.array([-1.0]),
        free_energy=np.array([[-0.9]]),
        isochoric_heat_capacity=np.array([[1.0]]),
    )

    observer(
        Event(
            "Settings",
            level=EventLevel.RESULT,
            data={"kind": "settings", "options": HAOptions()},
        )
    )
    observer(
        Event(
            "Input",
            level=EventLevel.RESULT,
            data={"kind": "input", "input": input_data},
        )
    )
    observer(
        Event(
            "Progress",
            level=EventLevel.PROGRESS,
            data={"label": "F", "current": 1, "total": 2},
        )
    )
    observer(
        Event(
            "Thermo",
            level=EventLevel.RESULT,
            data={"kind": "thermodynamics", "result": result},
        )
    )
    observer.save()

    text = observer.text()
    assert "Selected options" in text
    assert "Input data" in text
    assert "Progress: F (1/2)" not in text
    assert "Thermodynamic properties" in text
    assert report.read_text(encoding="utf-8") == text


def test_text_observer_renders_multivolume_ha_result_without_shape_error():
    observer = HATextObserver(silent=True, show_progress=False)
    result = HAResult(
        jobname="MgO volume series",
        volume=np.array([18.0, 19.0]),
        temperature=np.array([0.0, 300.0]),
        static_energy=np.array([-274.0, -273.9]),
        zero_point_energy=np.array([[0.10, 0.09]]),
        thermal_energy=np.array([[0.0, 0.0], [0.02, 0.03]]),
        internal_energy=np.array([[-273.9, -273.81], [-273.88, -273.78]]),
        entropy=np.array([[0.0, 0.0], [0.001, 0.002]]),
        vibrational_free_energy=np.array([[0.10, 0.09], [0.08, 0.06]]),
        free_energy=np.array([[-273.9, -273.81], [-273.92, -273.84]]),
        isochoric_heat_capacity=np.array([[0.0, 0.0], [0.001, 0.002]]),
        metadata={
            "units": {
                "energy": "Ha",
                "entropy": "Ha cell^-1 K^-1",
                "heat_capacity": "Ha cell^-1 K^-1",
                "volume": "A^3",
                "temperature": "K",
            }
        },
    )

    observer(
        Event(
            "Thermodynamics",
            level=EventLevel.RESULT,
            data={"kind": "thermodynamics", "result": result},
        )
    )

    text = observer.text()
    assert "Static and zero-point energies" in text
    assert "18.00000000" in text
    assert "19.00000000" in text
    assert "Helmholtz free energy" in text
