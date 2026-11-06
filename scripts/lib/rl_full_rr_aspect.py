"""Selectable rectangular aspect for the existing RL_FULL individual arms.

Only the historical study restriction aspect==4 is lifted. Section, J, arms,
energies, boundary matrices, count and Ritz remain inherited/unchanged.
"""
from dataclasses import dataclass
import math
from scripts.lib import isotropic_two_springs_thickness_asymmetry as frozen


@dataclass(frozen=True)
class Model(frozen.Model):
    def __post_init__(self):
        frozen.base.Model.__post_init__(self)
        frozen.ratios(self.mu, self.tau)
        if self.variant != "RL_FULL" or not math.isfinite(self.aspect):
            raise ValueError("Finite positive rectangular aspect and full RL inertia required")


def make_model(**parameters):
    # Preserve the exact original constructor for the accepted aspect-four case.
    return frozen.Model(**parameters) if parameters.get("aspect", 4.) == 4 else Model(**parameters)
