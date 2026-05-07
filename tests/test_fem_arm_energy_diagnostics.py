from __future__ import annotations

from scripts.lib.tracked_bending_descendant_shapes import (
    arm_energy_diagnostics,
    collect_single_branch_shape,
    radius_from_epsilon,
)


def test_desc05_beta30_fem_arm_energy_uses_direct_element_stiffness() -> None:
    epsilon = 0.0025
    shape_case = collect_single_branch_shape(
        branch_id="bending_desc_05",
        beta_deg=30.0,
        mu_value=0.0,
        radius=radius_from_epsilon(epsilon, l_total=2.0),
    )

    energy = arm_energy_diagnostics(shape_case, epsilon=epsilon)

    assert abs(float(energy["axial_energy_fraction_from_arm_energies"]) - 0.052605126) < 5e-4
    assert abs(float(energy["right_axial_share"]) - 0.5) < 1e-9
