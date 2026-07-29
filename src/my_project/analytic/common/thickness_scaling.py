"""Value-only algebra for mass-preserving circular-radius scaling."""

from __future__ import annotations


def mass_preserving_denominator_squared(mu: float, eta: float) -> float:
    """Return the squared normalization denominator for radius factors."""
    return 1.0 + 2.0 * mu * eta + eta**2


def mass_preserving_tau_values_from_denominator(
    eta: float,
    denominator: float,
) -> tuple[float, float]:
    """Return the two radius factors for a caller-provided denominator."""
    return (1.0 - eta) / denominator, (1.0 + eta) / denominator


__all__ = [
    "mass_preserving_denominator_squared",
    "mass_preserving_tau_values_from_denominator",
]
