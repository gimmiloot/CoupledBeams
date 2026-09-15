"""Source-specific Failla Table1 / Hong Tables2-3; no parameter sweep.

This is a separate diagnostic workflow: straight segmented beams/supports,
source rounding and staged gates cannot be a preset of the angled KV pilot.
Reuses K12 augmented Newton, complex SVD and atomic serialization only.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from pathlib import Path
import sys,time,json,csv,io,hashlib,subprocess,platform
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_literature_benchmarks as lit
from scripts.analysis.laminated_beams.pilot_inplane_kelvin_voigt import atomic,write_json,clean

OUTPUT=ROOT/'results/laminated_beams/inplane_kelvin_voigt_literature_benchmarks'
PDFS={'docs/literature/pdf/failla2014.pdf':'0d10c20796b1d35368e11aaad57764ae68b239b3f493a01361ea01782e3521e2',
      'docs/literature/pdf/hong1999.pdf':'e8622f7d407ffcd353fea8578cd40f840b530c33be272a9d1a81af47ac015088'}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_gate():
    for path,digest in PDFS.items():
        if sha(ROOT/path)!=digest:raise ValueError('PDF changed: recheck transcription '+path)
    return dict(status='SOURCE_TRANSCRIBED',date='2026-09-15',pdfs=PDFS,pages=[12,20],
        report='docs/laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md',
        extraction='Git pdftotext; Windows.Data.Pdf page rendering; manual formula/table verification',
        metadata='local first pages; Hong full DOI additionally confirmed at ScienceDirect')


def preflight():
    calls=kv.Calls();checks={}
    for p in (-.3+20j,1.+13j,-2.+800j):
        h=lit.Hong();arm=kv.Arm('RLB',h.E*h.A,h.D,h.m,h.L,1/h.S,h.J)
        H=kv.state_matrix(p,arm)[np.ix_([1,2,4,5],[1,2,4,5])]
        error=np.linalg.norm(lit.hong_state(p)-lit.H_TO_PROJECT@H@lit.H_TO_PROJECT)
        checks['hong_mapping_'+str(p)]=float(error)
    for case in ('failla','hong_hh','hong_ff','hong_damped'):
        beam=lit.Beam(case,calls);z=-.2+24j
        B,Bz=beam.matrices(z,derivative=True);step=1e-5
        fd=(beam.matrices(z+step)[0]-beam.matrices(z-step)[0])/(2*step)
        checks[case+'_derivative_relative']=float(np.linalg.norm(Bz-fd)/np.linalg.norm(Bz))
        checks[case+'_conjugacy_relative']=float(np.linalg.norm(beam.matrices(z.conjugate())[0]-B.conj())/np.linalg.norm(B))
    return dict(status='PASS' if max(checks.values())<=1e-7 else 'ASSEMBLY_MISMATCH',
                checks=checks,costs=calls.snapshot(),hong_properties=lit.Hong().data())


def comparison(case,mode,z,details):
    if case=='failla':
        re,im,ratio=lit.FAILLA[mode-1];value=-1j*z
        source_key,table=lit.FAILLA_KEY,'Table 1'
    else:
        source_key=lit.HONG_KEY
        if case=='hong_damped':re,im=lit.HONG_DAMPED[mode-1];table='Table 3'
        else:re=None;im=(lit.HONG_HH if case=='hong_hh' else lit.HONG_FF)[mode-1];table='Table 2'
        value=z/lit.Hong().time
    parts={'real':lit.rounding(re,value.real),'imag':lit.rounding(im,value.imag)}
    if case=='failla':parts['ratio']=lit.rounding(ratio,details['damping_ratio'])
    passed=all(v['rounding_pass'] for v in parts.values())
    status=details['status'] if details['status']!='CONVERGED' else ('PASS' if passed else 'LITERATURE_MISMATCH')
    inactive=None
    if case=='failla' and mode==4:
        inactive=max(v[key] for v in details['interfaces'] for key in ('v_abs','moment_abs','rotation_jump_abs'))<=lit.ZERO_ATOL and abs(z.real)<=lit.ZERO_ATOL
        if not inactive:status='INACTIVE_UNRESOLVED'
    row=dict(source_key=source_key,source_table=table,case=case,mode=mode,
        printed_value=f'{re if re is not None else "0 (asserted)"} + i*{im if im is not None else "0 (asserted)"}',
        computed_real=value.real,computed_imag=value.imag,z_real=z.real,z_imag=z.imag,
        solver_residual=details['solver_residual'],sigma_ratio=details['sigma_ratio'],
        physical_residual=details['physical_residual'],rounding_pass=passed,status=status,
        inactive_confirmed=inactive,iterations=details['steps'])
    for component,data in parts.items():
        for key,val in data.items():row[key+'_'+component]=val
    return row


def save_csv(path,rows):
    if not rows:return
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
    writer.writeheader();writer.writerows(rows);atomic(path,stream.getvalue())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--preflight-only',action='store_true')
    args=parser.parse_args();out=args.output.resolve()
    # Only the new benchmark directory (or descendants) may receive scientific output.
    if out!=OUTPUT and not out.is_relative_to(OUTPUT):raise ValueError('output outside literature benchmark directory')
    start=time.perf_counter()
    initial=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    dirty=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines()
    old=ROOT/'results/laminated_beams/inplane_kelvin_voigt_pilot'
    protected={str(p.relative_to(ROOT)):sha(p) for p in sorted(old.iterdir()) if p.is_file()}
    for name in ('scripts/lib/inplane_kelvin_voigt.py','scripts/analysis/laminated_beams/pilot_inplane_kelvin_voigt.py'):
        protected[name]=sha(ROOT/name)
    # Completed runs are reused explicitly; changing code/source requires a new output subdirectory.
    versions={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'scripts/lib/inplane_kelvin_voigt_literature_benchmarks.py')}
    prior=out/'run_manifest.json'
    if prior.exists():
        data=json.loads(prior.read_text(encoding='utf-8'))
        if data.get('finished'):
            if data['source_versions']!=versions or data['pdfs']!=PDFS:
                raise ValueError('completed benchmark provenance differs; do not replace existing results')
            if any(sha(ROOT/p)!=h for p,h in data['protected_sources'].items()):
                raise ValueError('protected sources changed')
            source_gate()
            print('REUSED_COMPLETED_BENCHMARK: zero matrix/root/shape calls')
            return
    data=dict(A0=source_gate(),A1=preflight(),rows={},roots={},stages={},protected_sources=protected)
    write_json(out/'benchmark_diagnostics.json',data)
    if data['A1']['status']!='PASS' or args.preflight_only:
        print('PREFLIGHT',data['A1']['status']);return
    shapes={}
    cases=('failla','hong_hh','hong_ff','hong_damped')
    for case in cases:
        if case=='hong_damped' and any(r['status']!='PASS' for k,rs in data['rows'].items() if k.startswith('hong_') for r in rs):
            data['stages'][case]=dict(status='NOT_RUN_B1_GATE',costs=kv.Calls().snapshot(),seconds=0.)
            data['rows'][case]=[dict(source_key=lit.HONG_KEY,source_table='Table 3',case=case,mode=i+1,
                printed_value=f'{r}+i*{im}',printed_value_real=r,printed_value_imag=im,
                status='NOT_RUN_B1_GATE',rounding_pass=None) for i,(r,im) in enumerate(lit.HONG_DAMPED)]
            save_csv(out/'hong1999_table3.csv',data['rows'][case]);break
        begun=time.perf_counter();beam=lit.Beam(case);rows=[]
        if case=='failla':guesses=[lit.failla_to_project(complex(float(a),float(b or 0))) for a,b,_ in lit.FAILLA]
        elif case=='hong_damped':guesses=[complex(float(a),float(b))*beam.time for a,b in lit.HONG_DAMPED]
        else:guesses=[1j*float(x)*beam.time for x in (lit.HONG_HH if case=='hong_hh' else lit.HONG_FF)]
        for mode,guess in enumerate(guesses,1):
            z,details,y=lit.solve(beam,guess)
            row=comparison(case,mode,z,details);rows.append(row)
            data['rows'][case]=rows;data['roots'][f'{case}_{mode}']=details
            if case=='failla':shapes[f'mode_{mode:02d}']=y
            data['stages'][case]=dict(status='RUNNING',costs=beam.calls.snapshot(),seconds=time.perf_counter()-begun)
            filename='failla2014_table1.csv' if case=='failla' else ('hong1999_table3.csv' if case=='hong_damped' else 'hong1999_table2.csv')
            export=(data['rows'].get('hong_hh',[])+data['rows'].get('hong_ff',[])) if 'table2' in filename else rows
            save_csv(out/filename,export)
            write_json(out/'benchmark_diagnostics.json',data)
            print(case,mode,row['status'],complex(row['computed_real'],row['computed_imag']),flush=True)
        data['stages'][case]=dict(status='PASS' if all(r['status']=='PASS' for r in rows) else 'PARTIAL',
                                 costs=beam.calls.snapshot(),seconds=time.perf_counter()-begun)
    if shapes:
        stream=io.BytesIO();np.savez_compressed(stream,xi=np.linspace(0,1,257),**shapes)
        atomic(out/'failla2014_shapes.npz',stream.getvalue())
    write_json(out/'benchmark_diagnostics.json',data)
    flat=[r for rs in data['rows'].values() for r in rs]
    manifest=dict(initial_HEAD=initial,task_initial_dirty=['?? docs/literature/pdf/failla2014.pdf','?? docs/literature/pdf/hong1999.pdf'],
        run_initial_dirty=dirty,final_dirty=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines(),
        pdfs=PDFS,pages={'failla2014.pdf':12,'hong1999.pdf':20},source_versions=versions,
        environment=dict(executable=sys.executable,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        command=' '.join(sys.argv),seconds=time.perf_counter()-start,stages=data['stages'],preflight=data['A1'],
        root_evaluations=sum(s['costs']['B'] for s in data['stages'].values()),
        criteria=dict(kv.CRITERIA,literature_reserve=lit.ROUNDING_RESERVE,asserted_zero_atol=lit.ZERO_ATOL),
        protected_sources=protected,protected_unchanged=all(sha(ROOT/p)==h for p,h in protected.items()),
        finished=True,status='PASS_WITH_SCOPE' if all(r['status']=='PASS' for r in flat) else 'PARTIAL',
        verified_eigenvalues=sum(r.get('computed_real') is not None for r in flat),
        rounding_pass_eigenvalues=sum(r.get('rounding_pass') is True for r in flat),
        limitations=['published targets used as predictors, no completeness certification',
            'no angled joint/axial/laminate reduction validation','independent published problems; shared corrector/expm'])
    write_json(prior,manifest)
    print(manifest['status'],manifest['verified_eigenvalues'],'eigenvalues',manifest['seconds'],'seconds')


if __name__=='__main__':main()
