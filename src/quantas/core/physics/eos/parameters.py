# -*- coding: utf-8 -*-

"""Physical parameters and implied-parameter rules for equations of state."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np

from .spec import EOSFamily, EOSModel, parse_eos_model

ArrayLike: TypeAlias = np.ndarray | Sequence[float]


@dataclass(frozen=True, slots=True)
class EOSParameters:
    """Store the physical parameters of an isothermal equation of state.

    Parameters
    ----------
    K0 : float
        Bulk modulus at the reference pressure.
    KP : float
        First pressure derivative of the bulk modulus at the reference state.
    KPP : float
        Second pressure derivative of the bulk modulus at the reference state.
    V0 : float
        Volume at the reference pressure.
    E0 : float or None, optional
        Energy at the reference volume for an integrated energy EOS.
    """

    K0: float
    KP: float
    KPP: float
    V0: float
    E0: float | None = None

    def __post_init__(self) -> None:
        """Validate finite and physically meaningful parameters."""
        values = [self.K0, self.KP, self.KPP, self.V0]
        if self.E0 is not None:
            values.append(self.E0)
        if not np.all(np.isfinite(np.asarray(values, dtype=np.float64))):
            raise ValueError("EOS parameters must be finite")
        if self.K0 <= 0.0:
            raise ValueError("K0 must be positive")
        if self.V0 <= 0.0:
            raise ValueError("V0 must be positive")

    def as_dict(self) -> dict[str, float | None]:
        """Return parameters as a serializable dictionary."""
        return {
            "E0": self.E0,
            "K0": self.K0,
            "KP": self.KP,
            "KPP": self.KPP,
            "V0": self.V0,
        }


def implied_kp(model: EOSModel) -> float | None:
    r"""Return the first bulk-modulus derivative implied by EOS order.

    Second-order Birch-Murnaghan and Tait equations impose
    :math:`K'_0=4`, the second-order natural-strain equation imposes
    :math:`K'_0=2`, and the optional second-order Vinet representation
    imposes :math:`K'_0=1`.

    Parameters
    ----------
    model : EOSModel
        Equation-of-state family and order.

    Returns
    -------
    float or None
        Implied value for a second-order model, otherwise ``None``.
    """
    if model.order != 2:
        return None
    if model.family in {EOSFamily.BIRCH_MURNAGHAN, EOSFamily.TAIT}:
        return 4.0
    if model.family is EOSFamily.NATURAL_STRAIN:
        return 2.0
    if model.family is EOSFamily.VINET:
        return 1.0
    return None


def implied_kpp(model: EOSModel, K0: float, KP: float) -> float:
    r"""Return the implied second pressure derivative of the bulk modulus.

    For third-order Birch-Murnaghan,

    .. math::

        K''_0 = -\frac{1}{K_0}\left[
        (3-K'_0)(4-K'_0)+\frac{35}{9}\right].

    For third-order natural strain,

    .. math::

        K''_0 = -\frac{1}{K_0}\left[
        1+(K'_0-2)+(K'_0-2)^2\right].

    The Vinet representation uses the Jeanloz relation,

    .. math::

        K''_0 = -\frac{1}{K_0}\left[
        \left(\frac{K'_0}{2}\right)^2
        +\frac{K'_0}{2}-\frac{19}{36}\right],

    while the truncated Tait form uses :math:`K''_0=-K'_0/K_0`.
    Murnaghan assumes :math:`K''_0=0`.

    Parameters
    ----------
    model : EOSModel
        Equation-of-state family and order.
    K0 : float
        Reference bulk modulus.
    KP : float
        First pressure derivative of the bulk modulus.

    Returns
    -------
    float
        Implied :math:`K''_0` in inverse-pressure units consistent with
        ``K0``.

    Raises
    ------
    ValueError
        If the EOS family is not supported.
    """
    if model.family is EOSFamily.MURNAGHAN:
        return 0.0
    if model.family is EOSFamily.BIRCH_MURNAGHAN:
        return -(((3.0 - KP) * (4.0 - KP)) + 35.0 / 9.0) / K0
    if model.family is EOSFamily.NATURAL_STRAIN:
        delta = KP - 2.0
        return -(1.0 + delta + delta**2) / K0
    if model.family is EOSFamily.VINET:
        return -((KP / 2.0) ** 2 + KP / 2.0 - 19.0 / 36.0) / K0
    if model.family is EOSFamily.TAIT:
        return -KP / K0
    raise ValueError(f"unsupported EOS family: {model.family.value}")


def resolve_energy_parameters(
    eos: str | EOSModel,
    parameters: ArrayLike | Mapping[str, float] | EOSParameters,
) -> EOSParameters:
    """Resolve free energy-fit parameters to a physical parameter set."""
    model = parse_eos_model(eos)
    if isinstance(parameters, EOSParameters):
        return parameters
    if isinstance(parameters, Mapping):
        values = {str(key).upper(): float(value) for key, value in parameters.items()}
    else:
        array = np.asarray(parameters, dtype=np.float64)
        if array.ndim != 1:
            raise ValueError("EOS parameters must form a one-dimensional array")
        names = model.energy_parameter_names
        if array.size == 5:
            names = ("E0", "K0", "KP", "KPP", "V0")
        elif array.size == 4 and len(names) != 4:
            names = ("E0", "K0", "KP", "V0")
        elif array.size != len(names):
            raise ValueError(
                f"{model.tag} energy parameters must contain {', '.join(model.energy_parameter_names)}"
            )
        values = dict(zip(names, array, strict=True))

    required = {"E0", "K0", "V0"}
    if not required.issubset(values):
        missing = ", ".join(sorted(required.difference(values)))
        raise ValueError(f"missing energy EOS parameters: {missing}")
    kp = implied_kp(model) if model.order == 2 else values.get("KP")
    if kp is None:
        raise ValueError(f"missing energy EOS parameter: KP for {model.tag}")
    if model.order == 4:
        if "KPP" not in values:
            raise ValueError(f"missing energy EOS parameter: KPP for {model.tag}")
        kpp = values["KPP"]
    else:
        kpp = implied_kpp(model, values["K0"], kp)
    return EOSParameters(
        E0=values["E0"],
        K0=values["K0"],
        KP=kp,
        KPP=kpp,
        V0=values["V0"],
    )


def resolve_pressure_parameters(
    eos: str | EOSModel,
    parameters: ArrayLike | Mapping[str, float] | EOSParameters,
) -> EOSParameters:
    """Resolve free pressure-fit parameters to a physical parameter set."""
    model = parse_eos_model(eos)
    if isinstance(parameters, EOSParameters):
        return parameters
    if isinstance(parameters, Mapping):
        values = {str(key).upper(): float(value) for key, value in parameters.items()}
    else:
        array = np.asarray(parameters, dtype=np.float64)
        if array.ndim != 1:
            raise ValueError("EOS parameters must form a one-dimensional array")
        pressure_names = model.pressure_parameter_names
        energy_names = model.energy_parameter_names if model.supports_energy else ()
        if array.size == len(pressure_names):
            names = pressure_names
        elif energy_names and array.size == len(energy_names):
            # Accept the free vector returned by the matching integrated EOS.
            names = energy_names
        elif array.size == 5:
            names = ("E0", "K0", "KP", "KPP", "V0")
        elif array.size == 4:
            names = ("K0", "KP", "KPP", "V0")
        else:
            raise ValueError(
                f"{model.tag} pressure parameters must contain "
                f"{', '.join(pressure_names)}"
            )
        values = dict(zip(names, array, strict=True))

    if not {"K0", "V0"}.issubset(values):
        raise ValueError("pressure EOS parameters must include K0 and V0")
    kp = implied_kp(model) if model.order == 2 else values.get("KP")
    if kp is None:
        raise ValueError(f"missing pressure EOS parameter: KP for {model.tag}")
    if model.order == 4:
        if "KPP" not in values:
            raise ValueError(f"missing pressure EOS parameter: KPP for {model.tag}")
        kpp = values["KPP"]
    else:
        kpp = implied_kpp(model, values["K0"], kp)
    return EOSParameters(K0=values["K0"], KP=kp, KPP=kpp, V0=values["V0"])


def resolved_energy_parameter_jacobian(
    eos: str | EOSModel,
    parameters: ArrayLike,
) -> np.ndarray:
    r"""Return the complete-parameter Jacobian for an energy EOS.

    The rows follow ``E0, K0, KP, KPP, V0`` and the columns follow the free
    energy-fit parameters of the selected EOS order.  For implied
    second-derivative parameters, analytical derivatives of the truncation
    relation are used.  The Jacobian therefore transforms a free-parameter
    covariance matrix into the covariance of the complete physical parameter
    set,

    .. math::

        \mathbf C_{\mathrm{resolved}}
        = \mathbf J\,\mathbf C_{\mathrm{free}}\,\mathbf J^{\mathsf T}.

    Parameters
    ----------
    eos : str or EOSModel
        Equation-of-state family and order.
    parameters : array-like
        Free energy-fit parameter vector.

    Returns
    -------
    ndarray
        Jacobian with shape ``(5, n_free)``.

    Raises
    ------
    ValueError
        If the free parameter vector has an invalid shape.
    """
    model = parse_eos_model(eos)
    free = np.asarray(parameters, dtype=np.float64)
    names = model.energy_parameter_names
    if free.ndim != 1 or free.size != len(names):
        raise ValueError(
            f"{model.tag} free energy parameters must contain {', '.join(names)}"
        )
    resolved = resolve_energy_parameters(model, free)
    row_names = ("E0", "K0", "KP", "KPP", "V0")
    row_index = {name: index for index, name in enumerate(row_names)}
    column_index = {name: index for index, name in enumerate(names)}
    jacobian = np.zeros((len(row_names), len(names)), dtype=np.float64)

    for name in ("E0", "K0", "KP", "KPP", "V0"):
        if name in column_index:
            jacobian[row_index[name], column_index[name]] = 1.0

    if "KPP" not in column_index:
        kpp_row = row_index["KPP"]
        k0_column = column_index["K0"]
        jacobian[kpp_row, k0_column] = -resolved.KPP / resolved.K0
        if "KP" in column_index:
            kp_column = column_index["KP"]
            if model.family is EOSFamily.BIRCH_MURNAGHAN:
                derivative = (7.0 - 2.0 * resolved.KP) / resolved.K0
            elif model.family is EOSFamily.NATURAL_STRAIN:
                derivative = -(2.0 * resolved.KP - 3.0) / resolved.K0
            elif model.family is EOSFamily.VINET:
                derivative = -(resolved.KP + 1.0) / (2.0 * resolved.K0)
            elif model.family is EOSFamily.TAIT:
                derivative = -1.0 / resolved.K0
            else:
                derivative = 0.0
            jacobian[kpp_row, kp_column] = derivative
    return jacobian


def resolved_energy_parameter_covariance(
    eos: str | EOSModel,
    parameters: ArrayLike,
    covariance: np.ndarray,
) -> np.ndarray:
    """Transform a free energy-fit covariance to physical EOS parameters.

    Parameters
    ----------
    eos : str or EOSModel
        Equation-of-state family and order.
    parameters : array-like
        Free energy-fit parameter vector.
    covariance : ndarray
        Covariance matrix of the free parameters.

    Returns
    -------
    ndarray
        Covariance matrix ordered as ``E0, K0, KP, KPP, V0``.

    Raises
    ------
    ValueError
        If ``covariance`` has an incompatible shape or non-finite values.
    """
    jacobian = resolved_energy_parameter_jacobian(eos, parameters)
    candidate = np.asarray(covariance, dtype=np.float64)
    expected = (jacobian.shape[1], jacobian.shape[1])
    if candidate.shape != expected:
        raise ValueError(
            f"free-parameter covariance must have shape {expected}, got {candidate.shape}"
        )
    if not np.all(np.isfinite(candidate)):
        raise ValueError("free-parameter covariance must be finite")
    transformed = jacobian @ candidate @ jacobian.T
    return 0.5 * (transformed + transformed.T)


def free_energy_parameters(model: EOSModel, parameters: EOSParameters) -> np.ndarray:
    """Return the free energy-fit vector for ``model``."""
    values = parameters.as_dict()
    return np.asarray(
        [values[name] for name in model.energy_parameter_names], dtype=np.float64
    )


def free_pressure_parameters(model: EOSModel, parameters: EOSParameters) -> np.ndarray:
    """Return the free pressure-fit vector for ``model``."""
    values = parameters.as_dict()
    return np.asarray(
        [values[name] for name in model.pressure_parameter_names], dtype=np.float64
    )


__all__ = [
    "EOSParameters",
    "free_energy_parameters",
    "free_pressure_parameters",
    "implied_kp",
    "implied_kpp",
    "resolve_energy_parameters",
    "resolve_pressure_parameters",
    "resolved_energy_parameter_covariance",
    "resolved_energy_parameter_jacobian",
]
