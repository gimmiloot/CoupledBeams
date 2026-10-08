"""Bounded full-family linear rectangular solid validation; frozen 1D physics.

One new monolithic C3D10 workflow reuses existing modal helpers. No nonlinear
step, old-result mutation, tuning or unrestricted eigenpair/mesh search.
"""
from __future__ import annotations
import argparse,csv,ctypes,datetime,hashlib,importlib.metadata,json,os,shutil,subprocess,sys,time
from pathlib import Path
if __name__ == '__main__':
    for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.analysis import solid_fem_single_rod_fixed_fixed as single
from scripts.analysis.thickness_mismatch.audits.audit_full_spectrum_3d_fem_smoke_extraction import parse_calculix_frequency_table
CONFIG=ROOT/'data/input/nlsp_linear_rectangular_3d_fem.json'
OUTPUT=ROOT/'results/nlsp_linear_rectangular_3d_fem'

"""Temporary reusable FEM-1 analytic preflight; imports existing frozen physics."""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import yartsev_ch2_monoclinic_rod as book
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section
from scripts.lib.weakly_nonlinear_spatial_rod import RodCoefficients

FAMILY_FIELD = {'axial_mh':(0,6), 'inplane_bending':(1,5), 'outplane_bending':(2,4)}


def build_preflight(config=None):
    config = {} if config is None else config
    geometry=config.get('geometry', {'L':1.,'b':.20,'h':.10})
    material=config.get('material', {'E':1.,'rho':1.,'nu':.3,'kappa':5/6})
    L,b,h=(float(geometry[name]) for name in ('L','b','h'))
    E,rho,nu,kappa=(float(material[name]) for name in ('E','rho','nu','kappa'))
    if (L,b,E,rho,nu,kappa)!=(1.,.20,1.,1.,.3,5/6) or h not in (.10,.12):
        raise ValueError('Only the authorized fixed control and its one backup thickness are supported')
    guard_multiplier=float(config.get('frequency_guard_multiplier',1.10))
    if guard_multiplier!=1.10:
        raise ValueError('Predeclared bounded guard is 10 percent')
    scan_ceiling=4.0
    policy={'scan_intervals':800,'root_xtol':1e-12,'root_rtol':1e-13}
    G=E/(2*(1+nu))
    isotropic=book.BookMaterial(E1_real=E,E2_real=E,G12_real=G,G13_real=G,G23_real=G,
        nu12=nu,rho=rho,eta1=0.,eta2=0.,eta12=0.,eta13=0.,eta23=0.)
    point=book.make_rod_point(0.,geometry=book.Geometry(a=b,b=h,length=L,shear_factor=kappa),
        material=isotropic,material_mode='elastic')
    if point.properties.Sbar16!=0 or point.torsion.C_T!=point.torsion.Cbar or point.torsion.C_T.imag!=0:
        raise ArithmeticError('Generalized isotropic C_T gate failed')
    ct=float(point.torsion.C_T.real)
    coefficients=RodCoefficients.rectangular(E,rho,nu,b,h,ct,kappa).values()
    models={
        'axial_mh': mh.project_jang_reduced_rectangular(rectangular_section(E=E,nu=nu,rho=rho,width=b,thickness=h,K=kappa)),
        'inplane_bending': mh.project_jang_reduced_rectangular(rectangular_section(E=E,nu=nu,rho=rho,width=b,thickness=h,K=kappa)),
        # Swapping section dimensions only calls the existing isotropic Tim block with I_perp.
        # The original MH contraction gradient and inertia remain the inplane values.
        'outplane_bending': mh.project_jang_reduced_rectangular(rectangular_section(E=E,nu=nu,rho=rho,width=h,thickness=b,K=kappa))}
    coefficient_mapping={
        'm':models['axial_mh'].coefficients['m'], 'jp':models['axial_mh'].coefficients['j'],
        'jb':models['outplane_bending'].coefficients['r'], 'C':models['axial_mh'].coefficients['C'],
        'H':models['axial_mh'].coefficients['H'], 'S':models['inplane_bending'].coefficients['S'],
        'Bp':models['inplane_bending'].coefficients['B'],'Bb':models['outplane_bending'].coefficients['B'],
        'CT':ct,'nu':nu}
    coeff_error=max(abs(coefficients[k]-v)/max(1.,abs(v)) for k,v in coefficient_mapping.items())
    reference_state_error=[]
    for omega in (.6,1.,3.15):
        old=book.state_matrix(omega,point).real[np.ix_((0,1,3,4),(0,1,3,4))]
        new=mh.harmonic_state_matrix(models['outplane_bending'],omega,'timoshenko')
        reference_state_error.append(float(np.max(np.abs(old-new)/np.maximum(1.,np.abs(old)))))
    if coeff_error>1e-14 or max(reference_state_error)>1e-14:
        raise ArithmeticError('Existing coefficient/operator mapping failed')
    root_sets,searches={},{}
    for family,model in models.items():
        block='mh' if family=='axial_mh' else 'timoshenko'
        root_sets[family],searches[family]=mh.finite_roots(model,L,block,1e-5,scan_ceiling,policy)
    if len(root_sets['axial_mh'])<1 or len(root_sets['inplane_bending'])<2 or len(root_sets['outplane_bending'])<2:
        raise ArithmeticError('Bounded analytic screening does not cover required family modes')
    first_torsion=math.pi/L*math.sqrt(ct/(coefficients['jp']+coefficients['jb']))
    required=max(root_sets['axial_mh'][0]['omega'],root_sets['inplane_bending'][1]['omega'],
                 root_sets['outplane_bending'][1]['omega'],first_torsion)
    ceiling=guard_multiplier*required
    if ceiling>scan_ceiling:
        raise ArithmeticError('Authorized initial scan ceiling does not cover predeclared guarded range')
    profiles,inventory,certificates={},[],{}
    grid=np.linspace(0.,L,401)
    for family,roots in root_sets.items():
        model=models[family]
        block='mh' if family=='axial_mh' else 'timoshenko'
        active=[r for r in roots if r['omega']<=ceiling]
        certificate=mh.finite_count_upper_bound(model,L,ceiling,block)
        certificate['found_distinct_modes']=len(active)
        certificate['saturated']=len(active)==certificate['upper_count']
        certificates[family]=certificate
        if not certificate['saturated']:
            raise ArithmeticError('Existing min-max completeness certificate was not saturated')
        for local,row in enumerate(active,1):
            omega=row['omega']
            mode=mh.finite_mode(model,L,omega,block,order=180)
            shooting,shotmeta=mh.transfer_boundary_matrix(model,L,omega,block)
            singular=np.linalg.svd(shooting,compute_uv=False)
            reference_residual=float(singular[-1]/max(singular[0],1e-30))
            diag=dict(mode['diagnostics'])
            diag.update({'analytic_singular_ratio':row['singular_ratio'],
                'independent_state_shooting_singular_ratio':reference_residual,
                'independent_state_shooting_method':shotmeta['method'],
                'independent_state_shooting_steps':shotmeta['steps']})
            if max(diag['boundary_scaled_residual'],diag['equation_scaled_residual'],
                   diag['energy_relative_error'],diag['analytic_singular_ratio'],reference_residual)>1e-8:
                raise ArithmeticError('Analytic eigenpair/reference/energy gate failed')
            values=mh.finite_state_basis(model,L,omega,grid,block)@mode['coefficients']
            q=np.zeros((len(grid),7))
            q[:,FAMILY_FIELD[family]]=values[:,:2]
            # Kinetic character is independent of arbitrary shape normalization.
            second_mass=model.coefficients['j' if block=='mh' else 'r']
            axial_fraction=float(mode['weights']@(model.coefficients['m']*mode['values'][:,0]**2))
            second_fraction=float(mode['weights']@(second_mass*mode['values'][:,1]**2))
            contraction_label=('acoustic_axial' if axial_fraction>=.5 else 'contraction_dominated') if block=='mh' else None
            record={'family':family,'local_mode':local,'omega':omega,'frequency_cycles_per_time':omega/(2*math.pi),
                'branch_character':contraction_label,'kinetic_primary_fraction':axial_fraction,
                'kinetic_secondary_fraction':second_fraction,'diagnostics':diag}
            inventory.append(record)
            profiles[f'{family}:{local}']={'x':grid,'q':q,'state':values,
                'analytic_basis_coefficients':mode['coefficients'],'family':family,'local_mode':local,'omega':omega}
    torsion_count=int(math.floor(ceiling/first_torsion))
    certificates['torsion']={'upper_count':torsion_count,'found_distinct_modes':torsion_count,
        'saturated':True,'method':'Exact scalar Dirichlet spectrum omega_n=n*pi/L*sqrt(C_T/(rho*Ip))'}
    for local in range(1,torsion_count+1):
        omega=local*first_torsion
        q=np.zeros((len(grid),7))
        q[:,3]=math.sqrt(2/(L*(coefficients['jp']+coefficients['jb'])))*np.sin(local*math.pi*grid/L)
        inventory.append({'family':'torsion','local_mode':local,'omega':omega,'frequency_cycles_per_time':omega/(2*math.pi),
            'branch_character':'generalized_condensed_warping','diagnostics':{'exact_scalar_spectrum':True}})
        profiles[f'torsion:{local}']={'x':grid,'q':q,'family':'torsion','local_mode':local,'omega':omega}
    inventory.sort(key=lambda row:row['omega'])
    for sorted_index,row in enumerate(inventory,1): row['sorted_index']=sorted_index
    axial=next(row for row in inventory if row['family']=='axial_mh' and row['branch_character']=='acoustic_axial')
    gate=axial['sorted_index']<=15
    result={'geometry':{'L':L,'b':b,'h':h},'material':{'E':E,'rho':rho,'nu':nu,'kappa':kappa},
        'section':{'A0':b*h,'I_parallel':b*h**3/12,'I_perp':h*b**3/12,'Ip':b*h*(b*b+h*h)/12},
        'coefficients':coefficients,'C_T_provenance':{'status':'PASS','source':'scripts/lib/yartsev_ch2_monoclinic_rod.py:generalized_torsional_stiffness',
            'equations':'Yartsev 2.8-2.10 existing generalized reduction','book_geometry':{'a':b,'b':h,'length':L},
            'isotropic_input':{'E1':E,'E2':E,'G12':G,'G13':G,'G23':G,'nu12':nu,'rho':rho,'theta_deg':0.},
            'Sbar16':0.,'Cbar':float(point.torsion.Cbar.real),'C_T':ct,
            'terms_used':point.torsion.terms_used,'series_sum':float(point.torsion.series_sum.real),
            'estimated_relative_tail':point.torsion.estimated_relative_tail,
            'substituted_G_Ip':False,'independent_dynamic_warping_field':False},
        'boundary_conditions':{'endpoints':[0.,L],'essential_fields':['u','w','v','Phi','psi','theta','c'],
            'values':0.,'rotation_clamp':'section_rotation','source_book_slope_clamp_used':False},
        'coefficient_mapping_error':coeff_error,'outplane_yartsev_state_mapping_max_error':max(reference_state_error),
        'screen_searches':searches,'screen_ceiling_omega':scan_ceiling,'search_policy':policy,
        'frequency_window':{'omega_required_max':required,'guard_multiplier':guard_multiplier,'omega_max':ceiling,
            'frequency_max_cycles_per_time':ceiling/(2*math.pi),'inventory_count':len(inventory),'counts_by_family':{k:v['found_distinct_modes'] for k,v in certificates.items()},
            'initial_requested_3d_eigenpairs':24,'maximum_one_extension_3d_eigenpairs':36},
        'geometry_gate':'PASS' if gate else 'GEOMETRY_GATE_PARTIAL','first_axial_sorted_index':axial['sorted_index'],
        'geometry_choice_basis':'Only 1D seven-field inventory before any 3D result',
        'merged_spectrum':inventory,'completeness':{'status':'PASS','method':'Distinct independent analytical roots saturate existing min-max upper counts blockwise; exact scalar torsion count',
            'certificates':certificates,'total_found':len(inventory),'total_upper_count':sum(c['upper_count'] for c in certificates.values())},
        'contraction_dominated_window':'OUTSIDE_BOUNDED_WINDOW','contraction_optical_cutoff_omega':2*math.pi*mh.blocks(models['axial_mh'])[0].cutoff_hz,
        'scope':'Linear accepted seven-field limit; no nonlinear trajectories; no 3D calls'}
    return result,profiles



"""Temporary geometry/audit/integration fragment for the scoped FEM-1 CLI."""
from itertools import combinations
from pathlib import Path
import re
import numpy as np

TET10_EDGES = ((0,1),(1,2),(0,2),(0,3),(1,3),(2,3))
TET4_FACES = ((0,1,2),(0,1,3),(0,2,3),(1,2,3))


def rectangular_geo(L: float, b: float, h: float, target_size: float) -> str:
    """Single centered global Box; NLSP B=diag(1,-1,-1), positive w=-Y, v=-Z."""
    if not all(np.isfinite(v) and v > 0 for v in (L,b,h,target_size)):
        raise ValueError('Finite positive geometry and target_size required')
    return f'''// Monolithic NLSP FEM-1 solid; no internal joint.
// Global X longitudinal, Y thickness, Z width. NLSP B=diag(1,-1,-1).
SetFactory("OpenCASCADE");
L = {L:.17g};
b = {b:.17g};
h = {h:.17g};
size = {target_size:.17g};
tol = 1e-8;
Box(1) = {{0, -h/2, -b/2, L, h, b}};
left_faces[] = Surface In BoundingBox {{-tol,-h/2-tol,-b/2-tol,tol,h/2+tol,b/2+tol}};
right_faces[] = Surface In BoundingBox {{L-tol,-h/2-tol,-b/2-tol,L+tol,h/2+tol,b/2+tol}};
Physical Volume("SOLID",1) = {{1}};
Physical Surface("FIXED_LEFT",2) = left_faces[];
Physical Surface("FIXED_RIGHT",3) = right_faces[];
Mesh.CharacteristicLengthMin = size;
Mesh.CharacteristicLengthMax = size;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 0;
Mesh.SaveAll = 0;
Mesh.MshFileVersion = 4.1;
'''


def tet10_quadrature() -> tuple[np.ndarray,np.ndarray]:
    """Positive symmetric 14-point rule exact through degree five, volume 1/6.

    Exact quadratic displacement mass products on straight affine C3D10
    geometry. Avoids negative corner weights of quadratic nodal lumping.
    """
    bary, weights = [], []
    for a, other, weight in ((.7217942490673264,.0927352503108912,.0122488405193937),
                              (.0673422422100982,.3108859192633006,.0187813209530026)):
        for i in range(4):
            row = np.full(4, other)
            row[i] = a
            bary.append(row)
            weights.append(weight)
    for pair in combinations(range(4),2):
        row = np.full(4,.04550370412564965)
        row[list(pair)] = .45449629587435035
        bary.append(row)
        weights.append(.0070910034628469)
    return np.asarray(bary), np.asarray(weights)


def tet10_shape(bary: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    """C3D10 N[point,node], dN[point,node,(r,s,t)] in Abaqus connectivity order."""
    bary = np.atleast_2d(np.asarray(bary,dtype=float))
    if bary.shape[1] != 4 or not np.allclose(bary.sum(axis=1),1):
        raise ValueError('Four barycentric coordinates summing to one required')
    dl = np.asarray(((-1.,-1.,-1.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)))
    N = np.empty((len(bary),10))
    dN = np.empty((len(bary),10,3))
    N[:,:4] = bary*(2*bary-1)
    dN[:,:4] = (4*bary[:,:,None]-1)*dl
    for k,(i,j) in enumerate(TET10_EDGES,4):
        N[:,k] = 4*bary[:,i]*bary[:,j]
        dN[:,k] = 4*(bary[:,i,None]*dl[j]+bary[:,j,None]*dl[i])
    return N,dN


def mesh_arrays(mesh) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """Sorted node IDs, xyz, sorted element IDs, zero-based connectivity."""
    ids = np.asarray(sorted(mesh.nodes),dtype=int)
    xyz = np.asarray([mesh.nodes[int(i)] for i in ids],dtype=float)
    eids = np.asarray(sorted(mesh.solid_elements),dtype=int)
    lookup = {int(n):i for i,n in enumerate(ids)}
    conn = np.asarray([[lookup[n] for n in mesh.solid_elements[int(e)]] for e in eids],dtype=int)
    if conn.ndim != 2 or conn.shape[1] != 10:
        raise ValueError('Every solid element must have ten C3D10 nodes')
    return ids,xyz,eids,conn


def quadrature_arrays(mesh,rho: float = 1.) -> dict[str,np.ndarray]:
    """Exact mass samples for straight-sided quadratic geometry.

    Uq=np.einsum('qi,eic->eqc',data['N'],Unodal[data['conn']]).
    xyz=[element,point,3], weights=[element,point], N=[point,10].
    """
    ids,xyz,eids,conn = mesh_arrays(mesh)
    points = xyz[conn]
    mids = np.asarray([(points[:,i]+points[:,j])/2 for i,j in TET10_EDGES]).transpose(1,0,2)
    if np.max(np.abs(points[:,4:]-mids)) > 1e-9*max(float(np.ptp(xyz,axis=0).max()),1.):
        raise ValueError('Nonaffine midside geometry needs higher order mass rule')
    bary,weights = tet10_quadrature()
    N,dN = tet10_shape(bary)
    dets = np.linalg.det(np.einsum('eic,qij->eqcj',points,dN))
    if np.any(dets <= 0):
        raise ValueError('Nonpositive quadratic element Jacobian')
    return {'node_ids':ids,'element_ids':eids,'conn':conn,'N':N,
            'xyz':np.einsum('qi,eic->eqc',N,points),
            'weights':rho*dets*weights[None,:],'barycentric':bary}


def _solid_types(path: Path) -> list[str]:
    out = []
    for line in path.read_text(encoding='utf-8',errors='strict').splitlines():
        if line.lstrip().upper().startswith('*ELEMENT'):
            match = re.search(r'\bTYPE\s*=\s*([^,\s]+)',line,re.I)
            if match and match[1].upper().startswith('C3D'):
                out.append(match[1].upper())
    return out


def _face_connectivity(conn: np.ndarray) -> dict:
    parent = list(range(len(conn)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    faces, quadratic_face_nodes = {}, {}
    for e,row in enumerate(conn):
        for face in TET4_FACES:
            key = tuple(sorted(int(row[i]) for i in face))
            faces.setdefault(key,[]).append(e)
            face_nodes = set(int(row[i]) for i in face)
            face_nodes.update(int(row[k]) for k,(i,j) in enumerate(TET10_EDGES,4) if i in face and j in face)
            quadratic_face_nodes.setdefault(key,[]).append(face_nodes)
    for owners in faces.values():
        for other in owners[1:]:
            parent[find(other)] = find(owners[0])
    return {'volume_face_connected_components':len({find(i) for i in range(len(conn))}),
            'boundary_triangle_count':sum(len(v)==1 for v in faces.values()),
            'maximum_face_incidence':max(map(len,faces.values())),
            'nonconforming_quadratic_face_count':sum(any(v!=sets[0] for v in sets[1:]) for sets in quadratic_face_nodes.values())}


def _resolution_metrics(corners: np.ndarray,axis: int,L: float,extent: float,other_extent: float) -> dict:
    """Actual corner-tetra intersections of nine interior transverse rays.

    Counts are intersected tetra cells, not claims of structured layer count.
    Offset rays avoid symmetry-plane/internal-face coincidences.
    """
    inv = np.linalg.inv((corners[:,1:]-corners[:,:1]).transpose(0,2,1))
    other_axis = 3-axis
    start,direction = np.zeros(3),np.zeros(3)
    direction[axis] = 1
    slope3 = np.einsum('eij,j->ei',inv,direction)
    slopes = np.column_stack((-slope3.sum(axis=1),slope3))
    lines = []
    for xf in (.233,.507,.773):
        for of in (-.271,.113,.319):
            start[0],start[other_axis] = xf*L,of*other_extent
            at3 = np.einsum('eij,ej->ei',inv,start[None,:]-corners[:,0])
            intercepts = np.column_stack((1-at3.sum(axis=1),at3))
            lo,hi = np.full(len(corners),-extent/2),np.full(len(corners),extent/2)
            good = np.ones(len(corners),dtype=bool)
            for k in range(4):
                positive,negative = slopes[:,k]>1e-12,slopes[:,k]<-1e-12
                flat = ~(positive|negative)
                lo[positive] = np.maximum(lo[positive],-intercepts[positive,k]/slopes[positive,k])
                hi[negative] = np.minimum(hi[negative],-intercepts[negative,k]/slopes[negative,k])
                good[flat&(intercepts[:,k]<-1e-10)] = False
            good &= hi-lo > extent*1e-9
            intervals = sorted(zip(lo[good],hi[good]))
            covered,cursor = 0.,-extent/2
            for left,right in intervals:
                covered += max(0.,right-max(left,cursor))
                cursor = max(cursor,right)
            lengths = hi[good]-lo[good]
            lines.append({'x':float(start[0]),'other_transverse':float(start[other_axis]),
                          'intersected_tetra_count':int(good.sum()),'covered_fraction':float(covered/extent),
                          'median_chord':float(np.median(lengths)) if len(lengths) else None})
    counts = [s['intersected_tetra_count'] for s in lines]
    return {'method':'nine interior rays intersect corner tetrahedra; counts are tetra cells',
            'counts_min_median_max':[min(counts),float(np.median(counts)),max(counts)],'rays':lines}


def _edge_resolution(xyz: np.ndarray,corner_ids: np.ndarray,axis: int,L: float,extent: float,other_extent: float) -> dict:
    """Actual outer-edge segment counts from vertices, excluding midside nodes."""
    vertices = xyz[np.unique(corner_ids)]
    other_axis,tol = 3-axis,1e-8*max(L,extent,other_extent)
    edges = []
    for x in (0.,L):
        for other in (-other_extent/2,other_extent/2):
            mask = (np.abs(vertices[:,0]-x)<tol)&(np.abs(vertices[:,other_axis]-other)<tol)
            coords = np.unique(np.round(vertices[mask,axis],12))
            edges.append({'x':x,'other_transverse':other,'corner_vertex_count':int(len(coords)),
                          'linear_segment_count':max(0,int(len(coords))-1),'coordinates':coords.tolist()})
    return {'method':'actual corner vertices on four box edges; midsides excluded','edges':edges}


def audit_rectangular_mesh(inp_path: Path,L: float,b: float,h: float,rho: float,reader,gmsh_logs: str = '') -> tuple[dict,object]:
    """Enforce C3D10/geometry/connectivity/Jacobian gates before eigenanalysis."""
    failures,types = [],_solid_types(inp_path)
    if not types or set(types) != {'C3D10'}:
        raise ValueError(f'Expected only C3D10 solid blocks; found {types}')
    mesh = reader(inp_path)
    if not mesh.nodes or not mesh.solid_elements:
        raise ValueError('Empty nodes or solid elements')
    ids,xyz,eids,conn = mesh_arrays(mesh)
    if not np.all(np.isfinite(xyz)):
        raise ValueError('Nonfinite mesh coordinates')
    points = xyz[conn]
    if len(np.unique(conn)) != len(ids):
        failures.append('Mesh contains nodes unused by the solid volume')
    if any(len(set(row))!=10 for row in conn):
        failures.append('Duplicate connectivity node in C3D10')
    expected_min,expected_max = np.asarray((0.,-h/2,-b/2)),np.asarray((L,h/2,b/2))
    bbox_min,bbox_max,tol = xyz.min(axis=0),xyz.max(axis=0),1e-8*max(L,b,h)
    bbox_ok = np.allclose(bbox_min,expected_min,atol=tol,rtol=0) and np.allclose(bbox_max,expected_max,atol=tol,rtol=0)
    if not bbox_ok:
        failures.append('Bounding box differs from frozen geometry')
    left,right = ids[np.abs(xyz[:,0])<=tol],ids[np.abs(xyz[:,0]-L)<=tol]
    ends_ok = len(left)>=6 and len(right)>=6 and not np.intersect1d(left,right).size
    if not ends_ok:
        failures.append('Empty/degenerate/overlapping fixed-face sets')
    rank_ok = True
    for label,fixed in (('left',left),('right',right)):
        face = xyz[np.isin(ids,fixed),1:]
        if len(face)<3 or np.linalg.matrix_rank(face-face.mean(axis=0),tol=tol)!=2:
            rank_ok = False
            failures.append(f'{label} fixed face lacks three noncollinear points')
    connectivity = _face_connectivity(conn)
    if connectivity['volume_face_connected_components']!=1:
        failures.append('Solid is not one face-connected volume')
    if connectivity['maximum_face_incidence']>2:
        failures.append('Nonmanifold tetra face')
    if connectivity['nonconforming_quadratic_face_count']:
        failures.append('Shared tetra face has nonconforming quadratic midside nodes')
    bary,weights = tet10_quadrature()
    test_bary = np.vstack((np.eye(4),np.full((1,4),.25),bary))
    _,dN = tet10_shape(test_bary)
    dets = np.linalg.det(np.einsum('eic,qij->eqcj',points,dN))
    bad = np.any(dets<=0,axis=1)
    if np.any(bad):
        failures.append('Nonpositive sampled quadratic Jacobian')
    mids = np.asarray([(points[:,i]+points[:,j])/2 for i,j in TET10_EDGES]).transpose(1,0,2)
    midpoint_error = float(np.max(np.abs(points[:,4:]-mids)))
    straight = midpoint_error<=1e-9*max(L,b,h)
    if not straight:
        failures.append('Nonaffine quadratic midsides need stronger certification')
    volume = float(np.sum(dets[:,5:]*weights))
    expected_volume = L*b*h
    volume_error = abs(volume-expected_volume)/expected_volume
    if volume_error>1e-7:
        failures.append('Integrated volume/mass differs from frozen geometry')
    corners,quality_values = points[:,:4],[]
    for anchor in range(4):
        edges = corners[:,[i for i in range(4) if i!=anchor]]-corners[:,anchor,None]
        quality_values.append(np.abs(np.linalg.det(edges))/np.maximum(np.prod(np.linalg.norm(edges,axis=2),axis=1),np.finfo(float).tiny))
    quality = np.min(quality_values,axis=0)
    warnings = [s.strip() for s in gmsh_logs.splitlines() if re.search(r'warning|error',s,re.I)]
    if any(re.search(r'\berror\b|negative.*jacob|invalid.*element',s,re.I) for s in warnings):
        failures.append('Critical Gmsh error/invalid-element warning')
    resolutions = {}
    if not np.any(bad):
        for axis,name,extent,other in ((1,'thickness_h',h,b),(2,'width_b',b,h)):
            resolutions[name] = {'interior_ray_intersections':_resolution_metrics(corners,axis,L,extent,other),
                                  'box_edge_resolution':_edge_resolution(xyz,conn[:,:4],axis,L,extent,other)}
            if any(abs(s['covered_fraction']-1)>1e-7 for s in resolutions[name]['interior_ray_intersections']['rays']):
                failures.append(f'{name} ray does not cover complete section')
    return {'status':'FAIL' if failures else 'PASS','failures':failures,
            'nodes':len(ids),'c3d10_elements':len(eids),'solid_element_types':sorted(set(types)),
            'bbox_min':bbox_min.tolist(),'bbox_max':bbox_max.tolist(),'bbox_matches':bool(bbox_ok),
            'volume':volume,'expected_volume':expected_volume,'volume_relative_error':volume_error,
            'mass':rho*volume,'expected_mass':rho*expected_volume,
            'fixed_left_count':len(left),'fixed_right_count':len(right),
            'fixed_left_ids':left.tolist(),'fixed_right_ids':right.tolist(),
            'rigid_body_kinematic_constraint_check':'PASS' if ends_ok and rank_ok else 'FAIL',
            **connectivity,'straight_midpoint_geometry':bool(straight),'maximum_midpoint_error':midpoint_error,
            'quadratic_jacobian_sample_count_per_element':len(test_bary),
            'minimum_quadratic_jacobian':float(dets.min()),'negative_or_zero_jacobian_elements':int(bad.sum()),
            'jacobian_evidence':'affine J confirmed by straight midsides; corners/centroid/14 quadrature samples',
            'quality_metric':'minimum corner normalized determinant (absolute triple product / edge-norm product)',
            'quality_min_median_max':[float(quality.min()),float(np.median(quality)),float(quality.max())],
            'actual_resolution':resolutions,'gmsh_warnings_errors':warnings},mesh


"""FEM-1 diagnostics fragment; no solver calls or physical operator edits."""
import numpy as np
from scipy.optimize import linear_sum_assignment

NLSP_FRAME_SIGNS = np.array((1., -1., -1.))
NLSP_MODE_FAMILIES = ('axial_mh', 'inplane_bending', 'outplane_bending', 'torsion')


def nlsp_local_vectors(values):
    """Global Box coordinates to B=[t,n,t cross n]=diag(1,-1,-1)."""
    a = np.asarray(values, dtype=float)
    if a.shape[-1] != 3 or not np.all(np.isfinite(a)):
        raise ValueError('Expected finite Cartesian vectors')
    return a * NLSP_FRAME_SIGNS


def nlsp_lift_fields(profile_x, fields, global_points):
    """Linear mass-kinematic shape lift q=(u,w,v,Phi,psi,theta,c).

    Local rotation a=(Phi,-psi,theta), c scales thickness director only.
    This diagnostic lift does not assert equivalence to 3D constitutive theory.
    """
    grid, q = np.asarray(profile_x, float), np.asarray(fields, float)
    p = nlsp_local_vectors(global_points)
    if grid.ndim != 1 or len(grid)<2 or np.any(np.diff(grid)<=0):
        raise ValueError('1D profile grid must be strictly increasing')
    if q.shape!=(len(grid),7) or not np.all(np.isfinite(q)):
        raise ValueError('Expected seven finite canonical fields')
    xyz = p.reshape(-1,3)
    if np.min(xyz[:,0])<grid[0]-1e-10 or np.max(xyz[:,0])>grid[-1]+1e-10:
        raise ValueError('Samples outside 1D profile interval')
    sampled = np.column_stack([np.interp(xyz[:,0],grid,q[:,k]) for k in range(7)])
    d = sampled[:,:3].copy()
    rotation = sampled[:,3:6]*np.array((1.,-1.,1.))
    r = xyz.copy()
    r[:,0] = 0.
    d += np.cross(rotation,r)
    d[:,1] += sampled[:,6]*xyz[:,1]
    return (d*NLSP_FRAME_SIGNS).reshape(p.shape)


def nlsp_weighted_mac(first, second, weights):
    """Normalization/sign independent full-vector quadrature MAC."""
    a,b,w = np.asarray(first,float),np.asarray(second,float),np.asarray(weights,float)
    if a.shape!=b.shape or a.shape[:-1]!=w.shape:
        raise ValueError('MAC dimensions differ')
    if np.any(w<=0) or not np.all(np.isfinite(w)):
        raise ValueError('MAC needs positive finite mass weights')
    aa,bb,ab = [float(np.sum(w*np.sum(v,axis=-1))) for v in (a*a,b*b,a*b)]
    if min(aa,bb)<=np.finfo(float).tiny:
        return 0.
    return float(np.clip(ab*ab/(aa*bb),0.,1.))


def nlsp_validate_nodal_modes(node_ids, mode_shapes, expected_mode_ids):
    """Require actual complete FRD vectors; never impute missing nodal values."""
    ids = np.asarray(node_ids,dtype=int)
    if len(np.unique(ids))!=len(ids):
        raise ValueError('Duplicate mesh node IDs')
    expected=list(expected_mode_ids)
    if len(set(expected))!=len(expected):
        raise ValueError('Duplicate expected FEM mode IDs')
    missing=[int(k) for k in expected if k not in mode_shapes]
    if missing:
        raise ValueError(f'Missing FEM eigenvector modes: {missing}')
    nodal={}
    for k in expected:
        shape=mode_shapes[k]
        absent=set(ids.tolist())-set(shape)
        extra=set(shape)-set(ids.tolist())
        if absent or extra:
            raise ValueError(f'FEM mode {k}: missing_nodes={len(absent)}, extra_nodes={len(extra)}')
        a=np.array([shape[int(node)] for node in ids],dtype=float)
        if a.shape!=(len(ids),3) or not np.all(np.isfinite(a)) or not np.any(a):
            raise ValueError(f'Incomplete/nonfinite/zero FEM eigenvector mode {k}')
        nodal[int(k)]=a
    return nodal


def nlsp_evaluate_tet10_displacements(nodal_displacement, connectivity, shape_values):
    nodal,conn,basis=np.asarray(nodal_displacement,float),np.asarray(connectivity,int),np.asarray(shape_values,float)
    if conn.ndim!=2 or conn.shape[1]!=10 or basis.ndim!=2 or basis.shape[1]!=10:
        raise ValueError('Expected complete C3D10 connectivity and shape functions')
    return np.einsum('qi,eic->eqc',basis,nodal[conn],optimize=True)


def _nlsp_rigid_section_design(p):
    design=np.zeros((len(p),3,6))
    design[:,:,:3]=np.eye(3)
    design[:,0,4]=p[:,2]
    design[:,0,5]=-p[:,1]
    design[:,1,3]=-p[:,2]
    design[:,2,3]=p[:,1]
    return design


def nlsp_project_sections(global_points, global_displacement, mass_weights, length,
                          m, jp, jb, section_count=41, dominance=0.70, residual_limit=0.35):
    """Mass-weighted rigid-section fit D+Omega cross r plus explicit residual.

    Each disjoint x bin uses local linear variation of D and physical Omega;
    ordinary variation in x therefore is not mislabeled sectional warping.
    Profiles are evaluated at volume-mass centroids. Endpoints are zero from
    the audited face clamp. Effective c=<y Ures_y>/<y^2> is a thickness strain
    diagnostic, not an FEM DOF and not an asserted identity with M-H c.
    """
    xyz=nlsp_local_vectors(global_points).reshape(-1,3)
    v=nlsp_local_vectors(global_displacement).reshape(-1,3)
    weight=np.asarray(mass_weights,float).reshape(-1)
    if xyz.shape!=v.shape or len(weight)!=len(xyz) or np.any(weight<=0):
        raise ValueError('Invalid section arrays or nonpositive mass weights')
    if min(m,jp,jb,length)<=0 or section_count<4:
        raise ValueError('Invalid section physical coefficients/count')
    bins=np.minimum((xyz[:,0]/length*section_count).astype(int),section_count-1)
    profiles,rows=[],[]
    total=float(np.sum(weight*np.sum(v*v,axis=1)))
    residual_norm=axial_warp_norm=trans_residual_norm=0.
    families=np.zeros(4)
    for k in range(section_count):
        sel=bins==k
        if np.count_nonzero(sel)<12:
            raise ValueError(f'Insufficient volume samples in section bin {k}')
        p,u,w=xyz[sel],v[sel],weight[sel]
        xc=float(np.average(p[:,0],weights=w))
        design=_nlsp_rigid_section_design(p)
        dx=(p[:,0]-xc)/length
        mat=np.concatenate((design,design*dx[:,None,None]),axis=2).reshape(-1,12)
        sw=np.repeat(np.sqrt(w),3)
        fitted,_,rank,_=np.linalg.lstsq(mat*sw[:,None],u.reshape(-1)*sw,rcond=1e-12)
        if rank!=12:
            raise ValueError(f'Rank-deficient section bin {k}: {rank}/12')
        residual=u-(mat@fitted).reshape(-1,3)
        resid=float(np.sum(w*np.sum(residual*residual,axis=1)))
        residual_norm+=resid
        axial_warp_norm+=float(np.sum(w*residual[:,0]**2))
        trans_residual_norm+=float(np.sum(w*np.sum(residual[:,1:]**2,axis=1)))
        dy,dz=float(np.sum(w*p[:,1]**2)),float(np.sum(w*p[:,2]**2))
        cy=float(np.sum(w*p[:,1]*residual[:,1])/dy) if dy else 0.
        cz=float(np.sum(w*p[:,2]*residual[:,2])/dz) if dz else 0.
        D,omega=fitted[:3],fitted[3:6]
        profiles.append(np.array((D[0],D[1],D[2],omega[0],-omega[1],omega[2],cy)))
        families+=(float(np.sum(w))/m)*np.array((m*D[0]**2,m*D[1]**2+jp*omega[2]**2,
                                               m*D[2]**2+jb*omega[1]**2,(jp+jb)*omega[0]**2))
        rows.append({'x':xc,'mass_weight':float(np.sum(w)), 'effective_thickness_contraction':cy,
                     'effective_width_strain':cz,'residual_norm':resid,'fit_rank':int(rank)})
    frac=families/max(float(np.sum(families)),np.finfo(float).tiny)
    residual_fraction=residual_norm/max(total,np.finfo(float).tiny)
    dominant=int(np.argmax(frac))
    family=('cross_section_or_local' if residual_fraction>residual_limit else
            'mixed_or_ambiguous' if frac[dominant]<dominance else NLSP_MODE_FAMILIES[dominant])
    return {'x':np.array([0.]+[r['x'] for r in rows]+[length]),
            'fields':np.vstack((np.zeros(7),profiles,np.zeros(7))), 'section_rows':rows,
            'family':family,'family_fractions':dict(zip(NLSP_MODE_FAMILIES,frac.tolist())),
            'residual_fraction':float(residual_fraction),
            'axial_warp_fraction':float(axial_warp_norm/max(total,np.finfo(float).tiny)),
            'transverse_residual_fraction':float(trans_residual_norm/max(total,np.finfo(float).tiny)),
            'mass_norm':total,'weight_source':'positive C3D10 degree-5 volume quadrature times rho',
            'contraction_status':'DIAGNOSTIC_EFFECTIVE_THICKNESS_STRAIN_NOT_MH_DOF'}


def nlsp_section_profile_mac(first, second, length, m, jp, jb, common_count=241):
    """Common-grid translation/rotation mass MAC; diagnostic c_eff is excluded."""
    x=np.linspace(0.,length,common_count)
    metric=np.sqrt(np.array((m,m,m,jp+jb,jb,jp)))
    a=np.column_stack([np.interp(x,first['x'],first['fields'][:,k]) for k in range(6)])*metric
    b=np.column_stack([np.interp(x,second['x'],second['fields'][:,k]) for k in range(6)])*metric
    w=np.ones(common_count)*length/(common_count-1)
    w[[0,-1]]*=.5
    return nlsp_weighted_mac(a,b,w)


def nlsp_shape_assignment(mac_matrix, min_mac=0.70, min_margin=0.08):
    """One-to-one shape-only assignment; retain raw independent-max conflicts.

    Frequency values are not accepted. Duplicated candidates, low MAC/margin
    and assignments displaced from independent best remain diagnostic only.
    """
    a=np.asarray(mac_matrix,float)
    if a.ndim!=2 or not np.all(np.isfinite(a)) or np.any((a<0)|(a>1+1e-12)):
        raise ValueError('Invalid shape MAC matrix')
    if a.shape[1]==0:
        return [{'row':k,'column':None,'status':'MISSING_FEM_MODE','mac':None,'margin':None,
                 'independent_best':None,'conflict':False} for k in range(a.shape[0])]
    best=np.argmax(a,axis=1)
    counts=np.bincount(best,minlength=a.shape[1])
    rows,columns=linear_sum_assignment(-a)
    assignment=dict(zip(rows.tolist(),columns.tolist()))
    out=[]
    for row in range(a.shape[0]):
        col=assignment.get(row)
        ordered=np.sort(a[row])
        margin=float(ordered[-1]-(ordered[-2] if len(ordered)>1 else 0.))
        conflict=bool(counts[best[row]]>1)
        mac=float(a[row,col]) if col is not None else None
        status=('MISSING_FEM_MODE' if col is None else 'DUPLICATE_INDEPENDENT_MATCH' if conflict else
                'ASSIGNMENT_NOT_INDEPENDENT_BEST' if col!=best[row] else 'LOW_MAC' if mac<min_mac else
                'AMBIGUOUS_SHAPE_OR_SUBSPACE' if margin<min_margin else 'MATCHED')
        out.append({'row':row,'column':col,'mac':mac,'margin':margin,'independent_best':int(best[row]),
                    'conflict':conflict,'status':status,'low_mac':mac is None or mac<min_mac,'ambiguous':margin<min_margin})
    return out


def nlsp_subspace_overlap(first_modes, second_modes, weights):
    """Principal-angle diagnostic, without asserting individual identities."""
    a,b=np.asarray(first_modes,float),np.asarray(second_modes,float)
    w=np.asarray(weights,float)
    if a.shape[1:]!=b.shape[1:] or a.shape[1:-1]!=w.shape:
        raise ValueError('Subspace samples/weights differ')
    sw=np.repeat(np.sqrt(w.reshape(-1)),a.shape[-1])
    aa=a.reshape(a.shape[0],-1).T*sw[:,None]
    bb=b.reshape(b.shape[0],-1).T*sw[:,None]
    ua,sa,_=np.linalg.svd(aa,full_matrices=False)
    ub,sb,_=np.linalg.svd(bb,full_matrices=False)
    if not len(sa) or not len(sb) or min(sa[0],sb[0])==0:
        return {'status':'ZERO_SUBSPACE','principal_mac':[]}
    ua,ub=ua[:,sa>1e-10*sa[0]],ub[:,sb>1e-10*sb[0]]
    singular=np.linalg.svd(ua.T@ub,compute_uv=False)
    return {'status':'DIAGNOSTIC_SUBSPACE_ONLY','first_rank':int(ua.shape[1]),'second_rank':int(ub.shape[1]),
            'principal_mac':np.clip(singular**2,0.,1.).tolist()}


def nlsp_relative_frequency_difference(omega_1d,omega_3d):
    if min(omega_1d,omega_3d)<=0:
        raise ValueError('Angular frequencies must be positive')
    signed=(float(omega_1d)-float(omega_3d))/float(omega_3d)
    return {'signed_relative_difference':signed,'absolute_relative_difference':abs(signed)}



def nlsp_read_frd_modes(path):
    """Read modal DISP fixed-width records without merging node ID and Ux.

    Historical regex tokenization is unsafe when positive Ux is adjacent to
    the ten-column node ID. Scope this strict reader to the new workflow;
    existing FEM scripts/results remain unchanged.
    """
    from pathlib import Path
    modes={}
    current=None
    collecting=False
    for line_number,raw in enumerate(Path(path).read_text(encoding='utf-8',errors='strict').splitlines(),1):
        s=raw.strip()
        if s.startswith('1PMODE'):
            current=int(s.split()[-1])
            collecting=False
        elif s.startswith('-4') and 'DISP' in s.upper():
            if current is None:
                raise ValueError(f'DISP block without modal ID at FRD line {line_number}')
            if current in modes:
                raise ValueError(f'Duplicate FEM mode DISP block {current}')
            modes[current]={}
            collecting=True
        elif collecting and raw[:3].strip()=='-1':
            if len(raw)<49:
                raise ValueError(f'Incomplete FRD displacement record at line {line_number}')
            try:
                node=int(raw[3:13])
                payload=raw[13:].replace('D','E').replace('d','E')
                import re
                # E12.5E3 negative values may overflow to 13 chars.
                number=re.compile(r'[-+]?(?:\d+\.\d*|\.\d+)[Ee][-+]\d{3}')
                tokens=number.findall(payload)
                if len(tokens)==3 and not number.sub('',payload).strip():
                    u=tuple(float(token) for token in tokens)
                else:
                    # Conventional E12.5 fixed-width export records.
                    u=tuple(float(raw[start:start+12].replace('D','E').replace('d','E')) for start in (13,25,37))
            except ValueError as exc:
                raise ValueError(f'Invalid fixed-width FRD record at line {line_number}') from exc
            if node in modes[current]:
                raise ValueError(f'Duplicate node {node} in FEM mode {current}')
            if not np.all(np.isfinite(u)):
                raise ValueError(f'Nonfinite FRD displacement at line {line_number}')
            modes[current][node]=u
        elif collecting and (s.startswith('-3') or s.startswith('-4') or s.startswith('100C') or s.startswith('1P')):
            collecting=False
    if not modes or any(not mode for mode in modes.values()):
        raise ValueError('Missing or empty FEM modal DISP blocks')
    return modes





def json_value(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)

def read_json(path):return json.loads(Path(path).read_text(encoding='utf8'))
def write_json(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(value,indent=2,default=json_value)+'\n',encoding='utf8')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def identity(config_path=CONFIG):
    c=read_json(config_path)
    item={'schema':c['schema'],'config':c,'config_sha256':sha(config_path),
          'code_sha256':sha(__file__),'model_sha256':{p:sha(ROOT/p) for p in c['model_hashes']},
          'executables':{k:sha(c[k]) for k in ('gmsh_exe','ccx_exe')},
          'python':sys.version,'dependencies':{k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item

def validate_config(c):
    if c['geometry']!={'L':1.,'b':.2,'h':.1} or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}:raise ValueError('Frozen pre-FEM chosen geometry/material mismatch')
    if c['element']!='C3D10' or c['initial_eigenpairs']!=24 or c['one_extension_eigenpairs']!=36:raise ValueError('Bounded mesh/eigenpair policy mismatch')
    if [v['name'] for v in c['mesh_levels']]!=['coarse','medium','fine']:raise ValueError('Exactly three predefined mesh levels')
    for v,n in zip(c['mesh_levels'],(2,3,4)):
        if v['target_size']!=c['geometry']['h']/n:raise ValueError('Mesh resolution changed after freeze')
    for path,digest in c['model_hashes'].items():
        if sha(ROOT/path)!=digest:raise ValueError('Frozen model source changed: '+path)
    if c['numerical_budget_seconds']>3600 or not c['semantics']['linear_only']:raise ValueError('Unauthorized computation policy')
    return c

def artifact_manifest(bundle,item):
    paths=sorted(p for p in bundle.rglob('*') if p.is_file() and p.name!='manifest.json')
    return {'identity':item,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'artifact_hashes':{p.relative_to(bundle).as_posix():sha(p) for p in paths}}

def validate_cache(bundle,item=None):
    b=Path(bundle);m=read_json(b/'manifest.json')
    if item is not None and m['identity']!=item:raise ValueError('Cache identity mismatch')
    for p,digest in m['artifact_hashes'].items():
        if sha(b/p)!=digest:raise ValueError('Cache artifact hash mismatch: '+p)
    return read_json(b/'summary.json') if (b/'summary.json').exists() else read_json(b/'preflight.json')

def write_preflight(c,bundle):
    start=time.perf_counter();validate_config(c);result,profiles=build_preflight(c)
    if result['geometry_gate']!='PASS' or abs(result['frequency_window']['omega_max']-c['frozen_omega_window'])>1e-12:raise ArithmeticError('Pre-FEM geometry/window gate not recovered')
    result['numerical_seconds']=time.perf_counter()-start
    result['configuration_frozen_before_fem']=True
    write_json(bundle/'preflight.json',result);write_json(bundle/'frozen_config.json',c)
    arrays={key.replace(':','_')+'_'+name:value for key,profile in profiles.items() for name,value in profile.items() if isinstance(value,np.ndarray)}
    np.savez_compressed(bundle/'one_D_profiles.npz',**arrays)
    return result,profiles

def load_profiles(bundle,pre):
    profiles={}
    with np.load(bundle/'one_D_profiles.npz',allow_pickle=False) as data:
        for row in pre['merged_spectrum']:
            key=f"{row['family']}:{row['local_mode']}";prefix=key.replace(':','_')
            profiles[key]={'x':data[prefix+'_x'].copy(),'q':data[prefix+'_q'].copy(),**row}
    return profiles

class WinMemory(ctypes.Structure):
    _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(k,ctypes.c_size_t) for k in ('peak_ws','ws','pool_peak','pool','nonpaged_peak','nonpaged','pagefile','peak_pagefile','private')]

def run_job(command,cwd,timeout,memory_limit,log_prefix,env):
    begin=time.perf_counter();peak=0;failure=None
    outpath=Path(str(log_prefix)+'.stdout.txt');errpath=Path(str(log_prefix)+'.stderr.txt')
    with outpath.open('w',encoding='utf8') as out,errpath.open('w',encoding='utf8') as err:
        process=subprocess.Popen(command,cwd=cwd,stdout=out,stderr=err,env=env,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        while process.poll() is None:
            if os.name=='nt':
                counters=WinMemory();counters.cb=ctypes.sizeof(counters)
                if ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(int(process._handle)),ctypes.byref(counters),counters.cb):peak=max(peak,int(counters.peak_ws))
            if time.perf_counter()-begin>timeout:failure='JOB_TIMEOUT'
            if peak>memory_limit:failure='JOB_MEMORY_LIMIT'
            if failure:process.kill();process.wait();break
            time.sleep(.1)
    stdout,stderr=outpath.read_text(encoding='utf8',errors='replace'),errpath.read_text(encoding='utf8',errors='replace')
    stats={'command':command,'returncode':process.returncode,'seconds':time.perf_counter()-begin,'peak_working_set_bytes':peak,'failure':failure}
    return subprocess.CompletedProcess(command,process.returncode,stdout,stderr),stats

def new_case_paths(directory,stem='modal'):
    return single.CasePaths(case_dir=directory,geo=directory/'rod.geo',msh=directory/'rod.msh',gmsh_inp=directory/'rod.inp',
        ccx_mesh_inp=directory/'solid_mesh.inp',ccx_modal_inp=directory/(stem+'.inp'),ccx_dat=directory/(stem+'.dat'),
        ccx_frd=directory/(stem+'.frd'),ccx_stdout=directory/(stem+'.stdout.txt'),ccx_stderr=directory/(stem+'.stderr.txt'))

def inspect_modal(c,pre,profiles,paths,audit,mesh,bundle):
    rows=parse_calculix_frequency_table(paths.ccx_dat)
    if not rows:raise ArithmeticError('No genuine DAT eigenfrequency table')
    for row in rows:
        if row['angular_frequency']<=0 or not np.isfinite(row['angular_frequency']):raise ArithmeticError('Nonpositive/nonfinite FEM frequency')
    units=[]
    for row in rows:
        omega=row['angular_frequency'];lam=row['eigenvalue'];freq=row['cyclic_frequency']
        if lam<=0 or freq<=0:raise ArithmeticError('Nonpositive independent DAT frequency columns')
        units.append(max(abs(omega**2/lam-1),abs(2*math.pi*freq/omega-1)))
    if max(units)>1e-5:raise ArithmeticError('DAT omega/sqrt(lambda)/2pi-f unit gate failed')
    modes=nlsp_read_frd_modes(paths.ccx_frd);quad=quadrature_arrays(mesh,c['material']['rho'])
    ids=[int(r['raw_mode_number']) for r in rows];nodal=nlsp_validate_nodal_modes(quad['node_ids'],modes,ids)
    fixed=np.isin(quad['node_ids'],audit['fixed_left_ids']+audit['fixed_right_ids'])
    coeff=pre['coefficients'];L=c['geometry']['L'];fits={};mac=np.zeros((len(pre['merged_spectrum']),len(rows)));sign=np.zeros_like(mac);shape_norm=np.zeros(len(rows));clamp=[]
    one=[nlsp_lift_fields(profiles[f"{r['family']}:{r['local_mode']}"]['x'],profiles[f"{r['family']}:{r['local_mode']}"]['q'],quad['xyz']) for r in pre['merged_spectrum']]
    one_norm=[float(np.sum(quad['weights']*np.sum(v*v,axis=-1))) for v in one]
    for j,modeid in enumerate(ids):
        values=nodal[modeid];res=float(np.max(abs(values[fixed]))/max(np.max(abs(values)),1e-30));clamp.append(res)
        if res>1e-8:raise ArithmeticError('FEM eigenvector violates fixed faces')
        disp=nlsp_evaluate_tet10_displacements(values,quad['conn'],quad['N']);shape_norm[j]=float(np.sum(quad['weights']*np.sum(disp*disp,axis=-1)))
        fits[modeid]=nlsp_project_sections(quad['xyz'],disp,quad['weights'],L,coeff['m'],coeff['jp'],coeff['jb'],c['matching']['section_count'],c['matching']['dominance'],c['matching']['section_residual_limit'])
        for i,predict in enumerate(one):
            mac[i,j]=nlsp_weighted_mac(predict,disp,quad['weights']);sign[i,j]=float(np.sum(quad['weights']*np.sum(predict*disp,axis=-1)))
    assignment=nlsp_shape_assignment(mac,c['matching']['minimum_mac'],c['matching']['minimum_margin']);matched=[]
    for a in assignment:
        row=dict(pre['merged_spectrum'][a['row']]);row.update(a);j=a['column']
        if j is not None:
            fr=rows[j];modeid=ids[j];fit=fits[modeid]
            row.update(fem_mode=modeid,omega_3d=float(fr['angular_frequency']),fem_family=fit['family'],section_residual_fraction=fit['residual_fraction'],axial_warp_fraction=fit['axial_warp_fraction'],full_shape_sign=float(np.sign(sign[a['row'],j])),one_D_lift_mass_norm=one_norm[a['row']],fem_mass_norm=float(shape_norm[j]))
            row.update(nlsp_relative_frequency_difference(row['omega'],row['omega_3d']))
            if row['status']=='MATCHED' and fit['family']!=row['family']:row['status']='SECTION_CHARACTER_UNRESOLVED'
        matched.append(row)
    matched_ids={r['fem_mode'] for r in matched if r['status']=='MATCHED'}
    additional=[{**r,'character':fits[int(r['raw_mode_number'])]['family'],'section_residual_fraction':fits[int(r['raw_mode_number'])]['residual_fraction']} for r in rows if r['angular_frequency']<=pre['frequency_window']['omega_max'] and r['raw_mode_number'] not in matched_ids]
    result={'raw_frequencies':rows,'maximum_printed_frequency_unit_error':max(units),'unit_gate_relative':1e-5,'frequency_unit':'rad per normalized time','eigenpairs':len(rows),'parsed_eigenvectors':len(nodal),'maximum_clamp_relative_residual':max(clamp),'matches':matched,'additional_modes_in_window':additional,'full_frequency_window_covered':rows[-1]['angular_frequency']>=pre['frequency_window']['omega_max'],'all_eight_identified':all(r['status']=='MATCHED' for r in matched),'matching_uses_frequencies':False,'legacy_FRD_reader_used':False,'strict_FRD_adapter_reason':'Confirmed fixed-width ID/payload concatenation defect in existing regex reader; frozen old code unchanged'}
    arrays={'node_ids':quad['node_ids'],'nodes':np.array([mesh.nodes[int(k)] for k in quad['node_ids']]),'MAC':mac,'full_shape_sign':sign,'clamp_residuals':np.array(clamp)}
    for modeid,value in nodal.items():arrays['U_'+str(modeid)]=value;arrays['profile_x_'+str(modeid)]=fits[modeid]['x'];arrays['profile_q_'+str(modeid)]=fits[modeid]['fields']
    np.savez_compressed(paths.case_dir/'modal_vectors.npz',**arrays)
    write_json(paths.case_dir/'section_diagnostics.json',{str(k):{name:v for name,v in fit.items() if name not in ('x','fields')} for k,fit in fits.items()})
    write_json(paths.case_dir/'modal_analysis.json',result)
    return result

def run_fem(c,bundle,pre,profiles):
    started=time.perf_counter();deadline=started+c['numerical_budget_seconds']-pre['numerical_seconds'];summary={'preflight':pre,'config':c,'meshes':{},'job_calls':{'gmsh':0,'ccx':0},'extensions_used':0,'statuses':{'NLSP_FEM1_GEOMETRY_SCREEN':'PASS','NLSP_FEM1_LINEAR_1D_SPECTRUM':'PASS'}}
    single.L=c['geometry']['L'];single.E=c['material']['E'];single.RHO=c['material']['rho'];single.NU=c['material']['nu'];single.SOLID_MODES_REQUESTED=c['initial_eigenpairs']
    env=dict(os.environ);env.update(OMP_NUM_THREADS='1',NUMBER_OF_CPUS='1')
    for level in c['mesh_levels']:
        name=level['name'];record={};summary['meshes'][name]=record
        if name=='fine':
            previous=summary['meshes']['medium'];ratio=(level['through_h_nominal']/3)**3
            estimate={'seconds':previous['total_seconds']*ratio**c['fine_forecast']['time_node_power'],'peak_memory_bytes':previous['peak_working_set_bytes']*ratio**c['fine_forecast']['memory_node_power'],'basis':'medium actual costs and predeclared size-volume scaling; checked before fine mesh'}
            summary['fine_forecast']=estimate;write_json(bundle/'fine_forecast.json',estimate)
            if estimate['seconds']>(deadline-time.perf_counter())*c['fine_forecast']['remaining_budget_fraction'] or estimate['peak_memory_bytes']>c['job_memory_limit_bytes']:
                record.update(status='NOT_RUN_RESOURCE_FORECAST');break
        case=bundle/'meshes'/name;case.mkdir(parents=True,exist_ok=True);paths=new_case_paths(case);paths.geo.write_text(rectangular_geo(**c['geometry'],target_size=level['target_size'])+'\nGeneral.NumThreads=1;\nMesh.RandomSeed=1;\n',encoding='utf8')
        jobs=[];original=single.subprocess.run
        def monitored(command,**kwargs):
            label='gmsh' if Path(command[0]).name.lower().startswith('gmsh') else 'ccx';summary['job_calls'][label]+=1
            remaining=deadline-time.perf_counter()
            if remaining<=0:raise TimeoutError('TOTAL_NUMERICAL_BUDGET')
            result,stats=run_job(command,kwargs.get('cwd',case),min(c['job_timeout_seconds'],remaining),c['job_memory_limit_bytes'],case/(label+'_job'+str(summary['job_calls'][label])),env);jobs.append(stats)
            if kwargs.get('check') and result.returncode:raise subprocess.CalledProcessError(result.returncode,command,result.stdout,result.stderr)
            return result
        single.subprocess.run=monitored
        stage=time.perf_counter()
        try:
            generated,message,warnings=single.generate_mesh_with_gmsh_cli(paths,c['gmsh_exe'])
            if not generated:raise ArithmeticError(message)
            logs='\n'.join(p.read_text(encoding='utf8',errors='replace') for p in case.glob('gmsh_job*.txt'))
            audit,mesh=audit_rectangular_mesh(paths.gmsh_inp,**c['geometry'],rho=c['material']['rho'],reader=single.read_gmsh_inp_mesh_data,gmsh_logs=logs);write_json(case/'mesh_audit.json',audit)
            record.update(mesh_audit=audit,target_size=level['target_size'])
            if audit['status']!='PASS':raise ArithmeticError('SOLID_MESH_GATE_FAILED')
            single.write_calculix_template(paths,mesh)
            outcome=single.run_calculix(paths,c['ccx_exe'],0.)
            if not outcome.success:raise ArithmeticError(outcome.message)
            modal=inspect_modal(c,pre,profiles,paths,audit,mesh,bundle);record.update(modal=modal)
            if not modal['full_frequency_window_covered'] and not summary['extensions_used']:
                summary['extensions_used']=1;single.SOLID_MODES_REQUESTED=c['one_extension_eigenpairs'];paths=new_case_paths(case,'modal_extended');single.write_calculix_template(paths,mesh)
                outcome=single.run_calculix(paths,c['ccx_exe'],0.)
                if not outcome.success:raise ArithmeticError(outcome.message)
                modal=inspect_modal(c,pre,profiles,paths,audit,mesh,bundle);record['modal']=modal
            record['status']='PASS' if modal['full_frequency_window_covered'] and modal['all_eight_identified'] else 'PARTIAL'
        except Exception as exc:record.update(status='FAIL',failure=str(exc))
        finally:
            single.subprocess.run=original;record.update(jobs=jobs,total_seconds=time.perf_counter()-stage,peak_working_set_bytes=max((j['peak_working_set_bytes'] for j in jobs),default=0));write_json(case/'case.json',record);write_json(bundle/'summary.json',summary)
        if record['status']!='PASS':break
    complete=all(summary['meshes'].get(n,{}).get('status')=='PASS' for n in ('coarse','medium','fine'))
    comparisons=[]
    if complete:
        for i,row in enumerate(pre['merged_spectrum']):
            matches=[summary['meshes'][n]['modal']['matches'][i] for n in ('coarse','medium','fine')];omegas=[r['omega_3d'] for r in matches];delta21=abs(omegas[1]-omegas[0])/omegas[1];delta32=abs(omegas[2]-omegas[1])/omegas[2]
            comparisons.append({'sorted_1d':row['sorted_index'],'family':row['family'],'local_mode':row['local_mode'],'omega_1d':row['omega'],'omega_coarse':omegas[0],'omega_medium':omegas[1],'omega_fine':omegas[2],'fem_mode_fine':matches[2]['fem_mode'],'shape_MAC_fine':matches[2]['mac'],'coarse_medium_relative':delta21,'medium_fine_relative':delta32,**nlsp_relative_frequency_difference(row['omega'],omegas[2]),'mesh_status':'PASS' if delta32<=c['numerical_mesh_convergence_relative'] and delta32<=delta21 else 'MESH_UNRESOLVED'})
    convergence=complete and all(r['mesh_status']=='PASS' for r in comparisons)
    summary['comparisons']=comparisons;summary['statuses'].update(NLSP_FEM1_3D_MESH_QUALITY='PASS' if complete else 'PARTIAL',NLSP_FEM1_3D_MODAL_EXECUTION='PASS' if complete else 'PARTIAL',NLSP_FEM1_3D_MESH_CONVERGENCE='PASS' if convergence else 'PARTIAL',NLSP_FEM1_MODE_IDENTIFICATION='PASS' if complete else 'PARTIAL',NLSP_FEM1_ALL_FAMILY_COMPARISON='PASS' if convergence else 'PARTIAL')
    summary['runtime']={'numerical_seconds':pre['numerical_seconds']+time.perf_counter()-started,'budget_seconds':3600};summary['qualification']='Linear only; finite 3D discretization is not exact truth and does not validate nonlinear V0 coefficients'
    write_json(bundle/'summary.json',summary)
    if comparisons:
        with (bundle/'all_family_comparison.csv').open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    return summary


def plot_only(bundle):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b=Path(bundle);s=validate_cache(b);pre=s['preflight'];fig=b/'figures';fig.mkdir(exist_ok=True)
    def finish(f,name):
        for ext in ('pdf','png'):f.savefig(fig/(name+'.'+ext),dpi=200,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(f)
    if not s['comparisons']:return {'solver_calls':0,'preflight_calls':0,'figures':0}
    rows=s['comparisons'];families=NLSP_MODE_FAMILIES;colors=dict(zip(families,('tab:red','tab:blue','tab:green','tab:orange')))
    f,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    for family in families:
        rr=[r for r in rows if r['family']==family];x=[r['sorted_1d'] for r in rr];ax.scatter(x,[r['omega_1d'] for r in rr],marker='o',facecolors='none',edgecolors=colors[family],label=family+' 1D');ax.scatter(x,[r['omega_fine'] for r in rr],marker='x',color=colors[family],label=family+' 3D')
    ax.set(xlabel='1D sorted position (shape matched)',ylabel='omega, rad / normalized time');ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2,frameon=False);finish(f,'combined_linear_spectrum')
    f,axes=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
    for family in families:
        rr=[r for r in rows if r['family']==family];x=[r['sorted_1d'] for r in rr];axes[0].plot(x,[100*r['absolute_relative_difference'] for r in rr],'o',label=family);axes[1].plot(x,[100*r['medium_fine_relative'] for r in rr],'o',label=family)
    axes[0].set(ylabel='1D / fine 3D difference, %');axes[1].set(ylabel='medium / fine change, %')
    for ax in axes:ax.set(xlabel='1D sorted position');ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
    finish(f,'frequency_difference_and_mesh_convergence')
    profiles=load_profiles(b,pre);f,axes=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    with np.load(b/'meshes/fine/modal_vectors.npz') as data:
        for family,field,ax in zip(families,(0,1,2,3),axes.flat):
            r=next(row for row in s['meshes']['fine']['modal']['matches'] if row['family']==family);p=profiles[family+':1'];x=data['profile_x_'+str(r['fem_mode'])];q=data['profile_q_'+str(r['fem_mode'])]
            ax.plot(p['x'],p['q'][:,field]/math.sqrt(r['one_D_lift_mass_norm']),label='1D');ax.plot(x,r['full_shape_sign']*q[:,field]/math.sqrt(r['fem_mass_norm']),'--',label='3D section fit');ax.set(title=family,xlabel='x/L',ylabel=('u','w','v','Phi')[field]);ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(f,'representative_all_family_shapes')
    return {'solver_calls':0,'preflight_calls':0,'figures':3}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);m=parser.add_mutually_exclusive_group(required=True);m.add_argument('--preflight',action='store_true');m.add_argument('--run-fem',action='store_true');m.add_argument('--report-only',type=Path);m.add_argument('--plot-only',type=Path);parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--output-dir',type=Path,default=OUTPUT);a=parser.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_cache(b)
        if a.plot_only:plot_only(b)
        print(json.dumps({'bundle':str(b),'statuses':s.get('statuses'),'new_solver_eigen_BVP_calls':0},indent=2));return s
    key,item=identity(a.config);b=a.output_dir/key;b.mkdir(parents=True,exist_ok=True)
    if (b/'manifest.json').exists():
        cached=validate_cache(b,item)
        if a.preflight or (b/'summary.json').exists():print(json.dumps({'bundle':str(b),'cache_hit':True,'new_solver_eigen_BVP_calls':0,'statuses':cached.get('statuses')},indent=2));return cached
    c=validate_config(read_json(a.config))
    if (b/'preflight.json').exists():pre=read_json(b/'preflight.json');profiles=load_profiles(b,pre)
    else:pre,profiles=write_preflight(c,b)
    if a.preflight:write_json(b/'manifest.json',artifact_manifest(b,item));print(json.dumps({'bundle':str(b),'geometry':pre['geometry'],'first_axial_sorted_index':pre['first_axial_sorted_index'],'frequency_window':pre['frequency_window']},indent=2));return pre
    (b/'execution_code').mkdir(exist_ok=True);shutil.copyfile(__file__,b/'execution_code'/Path(__file__).name)
    s=run_fem(c,b,pre,profiles);write_json(b/'manifest.json',artifact_manifest(b,item));plot_only(b);write_json(b/'manifest.json',artifact_manifest(b,item));print(json.dumps({'bundle':str(b),'statuses':s['statuses'],'runtime':s['runtime']},indent=2));return s

if __name__=='__main__':main()
