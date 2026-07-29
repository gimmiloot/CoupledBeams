"""Small algebra helpers shared by analytic model implementations."""

from .thickness_scaling import (
    mass_preserving_denominator_squared,
    mass_preserving_tau_values_from_denominator,
)

__all__ = [
    "mass_preserving_denominator_squared",
    "mass_preserving_tau_values_from_denominator",
]
