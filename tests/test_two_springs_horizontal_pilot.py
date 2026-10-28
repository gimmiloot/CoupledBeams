"""HS-P1 contracts, archive reuse, zero semantics and saved physical checks."""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil
import numpy as np
import pytest
from PIL import Image
from scripts.lib import two_springs_horizontal_pilot as p

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads(p.CONFIG.read_text())
PILOT=p.Pilot(CFG)
FOLDER=PILOT.folder


def read(name):
    with (FOLDER/name).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))


def test_exact_user_grid_and_policy():
    assert [(x['id'],x['k_joint'],x['kappa_horizontal']) for x in CFG['constructions']]==[
        ('H1',0,20000),('H2',.1,0),('H3',.1,1),('H4',.1,400),
        ('H5',.1,20000),('H6',.1,'RIGID'),('H7',10,20000)]
    assert CFG['angles_deg']==[0,.05,.1,.2,.5,1,2,5,10,20,30,45,60,75,90]
    policy=CFG['frequency_map']
    assert policy['frequency_map_policy']=='frequency-map-v1'
    assert policy['calculation_mode']=='fast_plot' and policy['spectrum_semantics']=='sorted_positions'
    assert (policy['K_plot'],policy['K_guard'])==(2,3)
    assert not policy['strict_audit_default']
    old=json.loads((ROOT/'data/input/isotropic_two_springs_horizontal_spring.json').read_text())
    assert CFG['policy']==old['policy']


def test_physics_verifier_and_historical_data_unchanged():
    for name,digest in CFG['frozen_files'].items():assert p.sha(ROOT/name)==digest
    source=ROOT/CFG['source_archive'];manifest=json.loads((source/'manifest.json').read_text())
    for name,digest in manifest['artifacts'].items():assert p.sha(source/name)==digest
    initial=json.loads((ROOT/CFG['initial_workspace']).read_text())
    assert len(initial['files'])==15
    for name in ['docs/theory/isotropic_two_springs_horizontal_spring.md',
                 'data/input/isotropic_two_springs_horizontal_spring.json',
                 'tests/test_two_springs_horizontal_support.py']:
        assert p.sha(ROOT/name)==initial['files'][name]


def test_full_grid_no_extra_angles_and_correct_reuse():
    manifest=PILOT.verify();rows=read('spectra.csv')
    assert manifest['status']=='PASS' and len(rows)==105
    assert manifest['new_distinct_points']==86 and manifest['archive_reused_points']==19
    assert manifest['added_angular_nodes']==[]
    assert len({(r['id'],r['beta_deg']) for r in rows})==105
    for construction in CFG['constructions']:
        assert [float(r['beta_deg']) for r in rows if r['id']==construction['id']]==CFG['angles_deg']


def test_zero_normalization_and_physical_prefixes():
    rows=read('spectra.csv')
    zero={r['id']:float(r['omega1']) for r in rows if float(r['beta_deg'])==0}
    for r in rows:
        record=json.loads((ROOT/r['point_file']).read_text())
        p.validate_record(record,p.physics.Case(**record['case']),CFG['policy'])
        for i in (1,2,3):
            w=float(r[f'omega{i}']);lam=float(r[f'Lambda{i}']);om=float(r[f'Omega{i}'])
            assert om==pytest.approx(lam**2,rel=2e-15)
            assert w==pytest.approx(.005*om,rel=2e-15,abs=1e-30)
        if r['id']=='H1':
            assert r['F_omega']=='' and r['D0_percent']==''
            assert (float(r['omega1'])==0)==(float(r['beta_deg'])==0)
        else:assert float(r['F_omega'])==pytest.approx(float(r['omega1'])/zero[r['id']],rel=1e-15)


def test_rejects_spurious_small_positive_replacement_for_zero():
    row=next(r for r in read('spectra.csv') if r['id']=='H1' and float(r['beta_deg'])==0)
    r=json.loads((ROOT/row['point_file']).read_text());r['omega'][0]=1e-8
    with pytest.raises(ArithmeticError,match='zero'):
        p.validate_record(r,p.physics.Case(**r['case']),CFG['policy'])


def test_monotonicity_and_triggered_neighbour_checks():
    data=json.loads((FOLDER/'monotonicity.json').read_text())
    assert data['status']=='PASS' and len(data['checks'])==90
    assert all(x['status']=='PASS' for x in data['checks'])
    audits=json.loads((FOLDER/'neighbour_audit.json').read_text())
    assert len(audits)==98 and sum(x['triggered'] for x in audits)==12
    assert all(x.get('ritz_frequency_error',0)<=CFG['policy']['ritz_frequency_tol'] for x in audits)


def test_energy_partition_and_normalized_selected_forms():
    rows=read('energy_distribution.csv')
    assert len(rows)==9
    for r in rows:
        assert float(r['mass'])==pytest.approx(1,abs=5e-7)
        assert sum(float(r[k+'_fraction']) for k in ['axial','bending','shear','joint','horizontal'])==pytest.approx(1,abs=2e-15)
        assert float(r['energy_error'])<=5e-8
        assert float(r['arm1_mass_fraction'])+float(r['arm2_mass_fraction'])==pytest.approx(1,abs=1e-15)
    mechanism=p.physics.mode(PILOT.model,p.physics.Case(0,k_joint=0,kappa_horizontal=20000),0)
    parts=p.energy_parts(PILOT.model,p.physics.Case(0,k_joint=0,kappa_horizontal=20000),mechanism,0)
    assert parts['energy_twice']==0 and parts['horizontal_fraction'] is None
    assert parts['rigid_rotation_residual']<1e-12


def test_local_formulas_use_radians_and_nonzero_baseline():
    co={r['id']:r for r in read('small_angle_coefficients.csv')}
    for r in read('small_angle_comparison.csv'):
        c=co[r['id']];predicted=float(c['omega0'])**2+float(c['coefficient_beta_rad_squared'])*math.radians(float(r['beta_deg']))**2
        assert float(r['omega2_approx'])==pytest.approx(predicted,rel=2e-15)
        assert float(r['signed_relative_error_lambda'])==pytest.approx(float(r['omega2_full'])/predicted-1,abs=1e-15)
    assert float(co['H1']['omega0'])==0
    assert float(co['H4']['coefficient_beta_rad_squared'])>float(co['H3']['coefficient_beta_rad_squared'])


def test_plot_only_never_calls_spectral_functions(tmp_path,monkeypatch):
    from scripts.lib.two_springs_horizontal_pilot_figures import render
    def forbidden(*args,**kwargs):raise AssertionError('Spectral calculation in plot_only')
    for name in ['prefix','count_at','mode','ritz','assembled']:monkeypatch.setattr(p.physics,name,forbidden)
    for name in ['manifest.json','analysis_manifest.json','spectra.csv','energy_distribution.csv',
                 'small_angle_comparison.csv','small_angle_coefficients.csv']:
        shutil.copy2(FOLDER/name,tmp_path/name)
    before=p.sha(tmp_path/'manifest.json');render(tmp_path)
    assert p.sha(tmp_path/'manifest.json')==before
    figmanifest=json.loads((tmp_path/'figures_manifest.json').read_text())
    assert figmanifest['new_eigenvalue_solves']==0
    assert len(list((tmp_path/'figures').glob('*.pdf')))==4
    for name in (tmp_path/'figures').glob('*.png'):
        with Image.open(name) as im:assert min(im.info['dpi'])>=300
