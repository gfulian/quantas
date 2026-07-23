"""Fit a volume-only Mie--Gruneisen--Debye P-V-T model with Quantas."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

from quantas.api import eos

# Edit these values for the dataset under study.
DATA_FILE = Path(__file__).resolve().parents[1] / "PVT_NaF.dat"
PRESSURE_MODEL = "BM4"
MGD_MODEL = "mgd"  # Alternative: "mgd:q-compromise".
FORMULA = "NaF"
REFERENCE_TEMPERATURE = 295.0  # K
SOLVER = "effective-variance"  # Alternative: "ols".

# Starting values for the thermal-pressure parameters.
THETA_D0_INITIAL = 459.0  # K
GAMMA0_INITIAL = 1.5
Q_INITIAL = 1.0

PARAMETER_UNITS = {
    "K0": "GPa",
    "KP": "",
    "KPP": "GPa^-1",
    "V0": "cm^3/mol",
    "temperature_ref": "K",
    "theta_d0": "K",
    "gamma0": "",
    "q": "",
}


def _solver_options():
    """Return the selected typed solver options."""
    if SOLVER == "ols":
        return eos.OLSOptions(max_iterations=10000)
    if SOLVER == "effective-variance":
        return eos.EffectiveVarianceOptions(
            max_iterations=50,
            inner_max_iterations=10000,
        )
    raise ValueError(f"unknown solver: {SOLVER!r}")


def _constraints() -> tuple[eos.ParameterConstraint, ...]:
    """Return explicit MGD starting values and physical bounds."""
    constraints: list[eos.ParameterConstraint] = [
        eos.ParameterConstraint.fixed("temperature_ref", REFERENCE_TEMPERATURE),
        eos.ParameterConstraint.free(
            "theta_d0",
            THETA_D0_INITIAL,
            lower_bound=1.0,
            upper_bound=5000.0,
        ),
        eos.ParameterConstraint.free(
            "gamma0",
            GAMMA0_INITIAL,
            lower_bound=0.01,
            upper_bound=10.0,
        ),
    ]
    if MGD_MODEL == "mgd":
        constraints.append(
            eos.ParameterConstraint.free(
                "q",
                Q_INITIAL,
                lower_bound=-10.0,
                upper_bound=10.0,
            )
        )
    return tuple(constraints)


def _format_optional(value: float | None) -> str:
    """Format an optional floating-point value for terminal output."""
    return "—" if value is None else f"{value:.10g}"


def _print_result(result) -> None:
    """Print parameters, diagnostics, and the strongest correlations."""
    fit = result.fit
    print("\nMGD fit summary")
    print("===============")
    print(f"Model              : {result.request.model.tag}")
    print(f"Solver             : {fit.method.value if fit.method else 'unknown'}")
    print(f"Status             : {fit.status.value}")
    print(f"Success            : {fit.success}")
    print(f"Message            : {fit.message}")

    if not fit.success or fit.parameters is None:
        return

    print("\nFitted parameters")
    print("-----------------")
    print(f"{'Parameter':<18} {'State':<10} {'Value':>18} {'E.S.D.':>18}  Unit")
    print(f"{'-' * 18} {'-' * 10} {'-' * 18} {'-' * 18}  {'-' * 12}")
    for index, (name, value) in enumerate(
        zip(fit.parameter_names, fit.parameters, strict=True)
    ):
        error = None if fit.errors is None else float(fit.errors[index])
        print(
            f"{name:<18} "
            f"{fit.parameter_states[index].value:<10} "
            f"{float(value):>18.10g} "
            f"{_format_optional(error):>18}  "
            f"{PARAMETER_UNITS.get(name, '')}"
        )

    diagnostics = fit.diagnostics
    print("\nFit diagnostics")
    print("---------------")
    print(f"Observations       : {fit.n_points}")
    print(f"Free parameters    : {fit.n_parameters}")
    print(f"Degrees of freedom : {fit.dof}")
    print(f"RMSE               : {_format_optional(fit.rmse)} GPa")
    print(f"MAE                : {_format_optional(fit.mae)} GPa")
    print(f"Maximum residual   : {_format_optional(fit.max_abs_error)} GPa")
    print(f"R-squared          : {_format_optional(fit.r_squared)}")
    if diagnostics is not None:
        print(
            f"Reduced chi-square : {_format_optional(diagnostics.reduced_chi_square)}"
        )
        print(f"Iterations         : {diagnostics.n_iterations or '—'}")

        correlation = diagnostics.correlation
        if correlation is not None:
            matrix = np.asarray(correlation, dtype=np.float64)
            pairs: list[tuple[float, str, str]] = []
            for row in range(matrix.shape[0]):
                for column in range(row + 1, matrix.shape[1]):
                    value = float(matrix[row, column])
                    if np.isfinite(value):
                        pairs.append(
                            (
                                abs(value),
                                fit.parameter_names[row],
                                fit.parameter_names[column],
                            )
                        )
            if pairs:
                pairs.sort(reverse=True)
                print("\nStrongest parameter correlations")
                print("--------------------------------")
                for absolute, first, second in pairs[:5]:
                    signed = matrix[
                        fit.parameter_names.index(first),
                        fit.parameter_names.index(second),
                    ]
                    print(f"{first:<18} {second:<18} {signed: .6f}")

    warnings = [*fit.warnings, *result.warnings]
    if warnings:
        print("\nWarnings")
        print("--------")
        for warning in dict.fromkeys(warnings):
            print(f"- {warning}")


def main() -> None:
    """Read one P-V-T dataset, perform the MGD fit, and report the result."""
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA_FILE
    dataset = eos.read_input(data_file)

    normalization = eos.MGDNormalization.molar_formula_unit(
        formula=FORMULA,
    )
    model = eos.PVTModel(
        pressure_model=PRESSURE_MODEL,
        coupling="thermal-pressure",
        thermal_pressure_model=MGD_MODEL,
        mgd_normalization=normalization,
    )
    request = eos.FitRequest(
        model=model,
        domain="pvt",
        target="volume",
        constraints=_constraints(),
        options=eos.FitOptions(solver_options=_solver_options()),
        metadata={"example": "volume-only MGD API fit"},
    )

    print("Quantas volume-only MGD fit")
    print("===========================")
    print(f"Dataset            : {dataset.jobname}")
    print(f"Source             : {dataset.source}")
    print(f"Observations       : {dataset.npoints}")
    print(f"Columns            : {', '.join(dataset.column_names)}")
    print(f"Atoms per unit     : {normalization.atoms_per_unit:g}")
    print(f"Reference T        : {REFERENCE_TEMPERATURE:g} K")

    result = eos.fit(dataset, request)
    _print_result(result)


if __name__ == "__main__":
    main()
