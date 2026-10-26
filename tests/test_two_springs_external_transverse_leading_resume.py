"""Leading continuation checks; tests do not solve the full dynamic model."""
import json
import math
import numpy as np
import pytest
from scripts.lib import two_springs_external_transverse_leading_resume as study

FOLDER=study.BASE/'dc0f91afbc14843c'

def test_full_and_saved_operation_guards():
    original=study.physics.prefix
    with study.forbid_full_calculations():
        for module in (study.physics,study.physics.base):
            for name in ('prefix','mode','ritz'):
                with pytest.raises(RuntimeError,match='forbidden'):getattr(module,name)()
        with pytest.raises(RuntimeError,match='forbidden'):study.api().Run.solve(None)
    assert study.physics.prefix is original
    with study.forbid_full_calculations(leading_too=True):
        with pytest.raises(RuntimeError,match='forbidden'):study.physics.Leading.solve(None,0.)

def test_diagnostic_recovery_keeps_original_failure():
    result=json.loads((FOLDER/'preflight.json').read_text())
    assert result['status']=='PASS' and result['cases']==10 and result['order_comparisons']==22
    assert result['old_failure']==1.0248928752260679e-9
    assert result['corrected']<1e-9 and result['control_xi10']<1e-9
    assert json.loads((FOLDER/'old_failure_reproduction.json').read_text())['status']=='EXPECTED_FAIL'

def test_only_actual_archived_nodes_and_correct_normalization():
    rows=study.read_csv(FOLDER/'comparison.csv');full=study.read_csv(FOLDER/'reused_full_points.csv')
    key=lambda r:(r['pair'],study.stiffness(r['k_w']),float(r['beta_deg']))
    assert len(rows)==len(full)==290
    assert {key(r) for r in rows}=={key(r) for r in full}
    zeros={key(r)[:2]:float(r['Lambda1']) for r in full if float(r['beta_deg'])==0}
    source={key(r):r for r in full}
    for r in rows:
        assert float(r['xi'])==pytest.approx(math.radians(float(r['beta_deg']))/.005,abs=1e-14)
        assert 0<=float(r['xi'])<=30
        assert float(r['F_full'])==pytest.approx((float(source[key(r)]['Lambda1'])/zeros[key(r)[:2]])**2,rel=2e-14)
        assert float(r['e'])==pytest.approx(float(r['F_full'])/float(r['F_lim'])-1,abs=1e-15)
        assert float(r['delta_percent'])==pytest.approx(100*abs(float(r['e'])),abs=1e-14)

def test_all_saved_form_gates_and_complete_prefixes():
    cfg=json.loads((FOLDER/'contract.json').read_text())['config']
    assert cfg['analytic_curvature'] and cfg['quality']['strong_ode_tol']==1e-9
    for path in (FOLDER/'leading_points').glob('*.json'):
        record=json.loads(path.read_text());identity=record['identity'];solution=record['solution']
        assert identity['analytic_curvature'] and record['status']=='PASS'
        for key in study.QUALITY:assert solution['diagnostics'][key]<=cfg['quality'][key+'_tol']
        assert np.all(np.diff(solution['zeta'])>0) and solution['zeta'][0]>0
        assert [r['negative_count'] for r in record['checks']['inertia']]==[0,1,2,3]
        for row in record['checks']['positions']:
            assert row['strong']<=1e-9 and row['natural']<=1e-9
            assert ('outer_transverse_force' in row['natural_raw'])==(identity['k_w']!='RIGID')

def test_curvature_from_mass_normalized_shape_and_controls():
    rows=study.read_csv(FOLDER/'leading_curvature.csv');assert len(rows)==10
    for r in rows:
        assert float(r['mass'])==pytest.approx(1.,abs=5e-7)
        assert float(r['C'])==pytest.approx(float(r['node0'])**2/(4*float(r['zeta0'])*float(r['mass'])),rel=1e-14)
    for r in study.read_csv(FOLDER/'local_curvature_checks.csv'):assert float(r['relative_error'])<=1e-4

def test_convergence_includes_all_maximum_errors():
    rows=study.read_csv(FOLDER/'convergence.csv')
    key=lambda r:(r['pair'],study.stiffness(r['k_w']),float(r['xi']))
    covered={key(r) for r in rows}
    for r in study.read_csv(FOLDER/'error_statistics.csv'):
        assert (r['pair'],study.stiffness(r['k_w']),float(r['max_xi'])) in covered
        assert int(r['nodes'])==29
    assert max(float(r['relative_frequency_error']) for r in rows)<1e-8

def test_monotonicity_and_no_full_recalculation():
    rows=study.read_csv(FOLDER/'monotonicity.csv')
    assert {r['parameter'] for r in rows}=={'xi','k_w'}
    assert all(float(r['min_relative_frequency_increment'])>=-1e-10 for r in rows)
    summary=json.loads((FOLDER/'numerical_summary.json').read_text())
    assert summary['status']=='PASS' and summary['full_model_solves']==0
    assert all(v['new_full_model_solves']==0 for v in json.loads((FOLDER/'executions.json').read_text()))

def test_source_archives_and_physical_modules_unchanged():
    cfg=json.loads(study.CONFIG.read_text())
    for name,value in cfg['sources'].items():
        folder=study.ROOT/value['folder']
        assert study.api().sha(folder/'manifest.json')==value['manifest_sha256']
        assert study.verify(folder)['unchanged']
    for name,digest in cfg['frozen_files'].items():assert study.api().sha(study.ROOT/name)==digest
