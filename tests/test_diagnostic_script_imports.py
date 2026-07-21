import importlib
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class DiagnosticScriptImportTest(unittest.TestCase):
    def test_thickness_mismatch_diagnostic_entrypoints_import(self) -> None:
        modules = [
            "scripts.analysis.plot_diagnostic_eta_mu_beta_frequency_maps",
            "scripts.analysis.plot_lambda_mu_eta_m0p5_with_single_beam_refs",
            "scripts.analysis.plot_mode_shapes_eta_beta_scan",
            "scripts.analysis.plot_mode_shapes_eta_beta_scan_sorted_modes",
            "scripts.analysis.thickness_mismatch.maps.plot_eb_vs_timoshenko_lambda_beta_cases",
            "scripts.analysis.thickness_mismatch.maps.plot_eb_vs_timoshenko_lambda_mu_cases",
            "scripts.analysis.thickness_mismatch.audits.analyze_universal_eb_validity_parameters_stage1",
            "scripts.analysis.thickness_mismatch.audits.audit_eb_validity_fixed_epsilon_geometry_scan",
            "scripts.analysis.thickness_mismatch.audits.audit_longitudinal_suspect_modes_eb_timo",
            "scripts.analysis.thickness_mismatch.audits.audit_eb_validity_vs_timoshenko_stage1",
            "scripts.analysis.thickness_mismatch.audits.audit_eb_epsilon_baseline_thresholds",
            "scripts.analysis.thickness_mismatch.audits.run_eb_epsilon_apriori_pilot",
            "scripts.analysis.thickness_mismatch.audits.audit_timoshenko_modes_4_6_shape_diagnostics",
            "scripts.analysis.thickness_mismatch.audits.audit_timoshenko_shape_bug_thin_limit",
            "scripts.analysis.thickness_mismatch.audits.audit_timoshenko_shape_construction",
            "scripts.analysis.thickness_mismatch.audits.audit_timoshenko_modes456_visualization",
            "scripts.analysis.thickness_mismatch.postprocess.analyze_eb_safe_prefix_certification",
            "scripts.analysis.thickness_mismatch.postprocess.analyze_eb_epsilon_apriori_pilot",
            "scripts.analysis.thickness_mismatch.shapes.plot_eb_timo_full_mode_shapes_eps0p03_beta45_eta0_modes4_6",
        ]
        for module_name in modules:
            with self.subTest(module=module_name):
                module = importlib.import_module(module_name)
                self.assertTrue(callable(getattr(module, "main", None)))


if __name__ == "__main__":
    unittest.main()
