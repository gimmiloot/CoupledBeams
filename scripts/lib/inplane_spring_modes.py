"""Two-arm EB/RLB diagnostics for the bounded laminate/length check.

Physical equations and spring assembly are delegated to existing modules.
Unlike the identical-arm EB helper, unequal arms are never parity-projected.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm

from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_rlb as rlb
from scripts.lib import inplane_rotational_spring_eb_modes as modes
from scripts.lib import reddy_symmetric_laminated_beam as laminate

REFERENCE = eb.EBArm(.20*.05, .20*.05**3/12, .20*.05, 1.)
FREQUENCY_SCALE = np.sqrt(REFERENCE.m/REFERENCE.D)


def section(contrast=.4):
    """Native ply integration; contrast=0 is a constitutive control only."""
    if contrast not in (0., .4):
        raise ValueError('only the fixed laminate and its constitutive baseline')
    plies = []
    for label, factor in zip('HLLH', (1+contrast,1-contrast,1-contrast,1+contrast)):
        material = laminate.OrthotropicLamina(
            E1=1.1*factor, E2=.9*factor, nu12=.3,
            G12=factor/2.6, G13=factor/2.6, G23=factor/2.6, rho=1., name=label)
        plies.append(laminate.Ply(material,0.,.05/4,label))
    integrated = laminate.integrate_laminate(plies)
    return integrated, laminate.reduce_to_beam_properties(integrated,width=.20,K=5/6)


@dataclass(frozen=True)
class Arm:
    model: str
    properties: laminate.BeamProperties
    L: float

    def __post_init__(self):
        if self.model not in ('EB','RLB'):
            raise ValueError('model must be EB or RLB')
        if not np.isfinite(self.L) or self.L <= 0:
            raise ValueError('positive finite length required')

    @property
    def native(self):
        p = self.properties
        return eb.EBArm(p.A,p.D,p.m,self.L) if self.model == 'EB' else rlb.LimitArm(p,self.L,1.)

    @property
    def rotational_mass(self):
        return self.properties.J if self.model == 'RLB' else 0.

    def matrix(self,omega):
        return (eb if self.model == 'EB' else rlb).state_matrix(omega,self.native)

    def scale(self):
        return (eb if self.model == 'EB' else rlb).state_scale(self.native)


def arms(model,mu,properties):
    if not np.isfinite(mu) or abs(mu)>=1:
        raise ValueError('finite |mu|<1 required')
    return Arm(model,properties,1-mu),Arm(model,properties,1+mu)


def assembly(omega,pair,beta,joint):
    if pair[0].model != pair[1].model:
        raise ValueError('both arms use the same model')
    module = eb if pair[0].model == 'EB' else rlb
    a,b = pair[0].native,pair[1].native
    if pair[0] == pair[1]:
        b = a
    return module.boundary_assembly(omega,a,b,beta,joint,REFERENCE)


def states_along_arm(omega,arm,reactions,nodes=129,check=False):
    """EB analytic functions, RLB constant-H recurrence in scaled variables."""
    reactions=np.asarray(reactions)
    if reactions.shape!=(3,) or not np.all(np.isfinite(reactions)):
        raise ValueError('three finite physical clamp reactions required')
    xi,_ = modes.quadrature(nodes)
    checks = []
    calls = dict(analytic=0,expm=0)
    if arm.model == 'EB':
        states = modes.arm_states(omega,arm.native,reactions,xi)
        calls['analytic'] = 1
    else:
        scale = arm.scale()
        h = arm.matrix(omega)*scale[None,:]/scale[:,None]
        step = expm(h*arm.L/(nodes-1))
        calls['expm'] += 1
        initial = np.r_[np.zeros(3),reactions]/scale
        work = initial.copy()
        states = np.empty((nodes,6),dtype=np.result_type(reactions,float))
        for i in range(nodes):
            states[i] = work*scale
            work = step@work
    if check:
        scale = arm.scale()
        h = arm.matrix(omega)*scale[None,:]/scale[:,None]
        initial = np.r_[np.zeros(3),reactions]/scale
        for index in (nodes//4,nodes//2,nodes-1):
            exact = expm(h*arm.L*xi[index])@initial
            calls['expm'] += 1
            error = np.linalg.norm(states[index]/scale-exact)
            checks.append(float(error/max(np.linalg.norm(exact),1e-30)))
    return states,calls,checks


def physical_vector(states,pair,weights):
    """One physical modal mass, with each actual L and (RLB only) J."""
    components = []
    for y,arm in zip(states,pair):
        density = np.array([arm.properties.m,arm.properties.m,arm.rotational_mass])
        components.append((y[:,:3]*np.sqrt(arm.L*weights[:,None]*density)).ravel())
    return np.concatenate(components)


def common_vector(states,weights):
    """Comparison metric, not RLB modal mass: fixed reference m,L and u,w."""
    vector = (np.asarray(states)[...,:2]*np.sqrt(REFERENCE.m*REFERENCE.L*weights)[None,:,None]).ravel()
    return vector/np.linalg.norm(vector)


def recover(assembled,omega,pair,beta,joint,parity=None,nodes=129,check=False):
    if parity is not None and pair[0] != pair[1]:
        raise ValueError('parity projection forbidden for unequal arms')
    endpoint = eb.endpoint_diagnostics(assembled,beta,joint,REFERENCE)
    if endpoint['nullity']<1 or endpoint['sigma_ratio']>1e-9:
        raise ValueError('FULL_MATRIX_ROOT_GATE')
    _,weights = modes.quadrature(nodes)
    possibilities = []
    for record in endpoint['vectors']:
        reactions = np.asarray(record['physical_clamp_reactions']).reshape(2,3)
        if parity is not None:
            projected = (reactions+parity*reactions[::-1]*[1,-1,-1])/2
            if np.linalg.norm(projected)<1e-6*np.linalg.norm(reactions):
                continue
            reactions = projected
        values = [states_along_arm(omega,a,r,nodes,check) for a,r in zip(pair,reactions)]
        states = np.array([v[0] for v in values])
        vector = physical_vector(states,pair,weights)
        mass = float(np.vdot(vector,vector).real)
        states,reactions,vector = states/np.sqrt(mass),reactions/np.sqrt(mass),vector/np.sqrt(mass)
        common = common_vector(states,weights)
        mirrored = common_vector(modes.reflect(states),weights)
        reflection = float(np.vdot(common,mirrored).real)
        eta = (1 if reflection >= 0 else -1) if pair[0] == pair[1] else None
        ends = states[:,-1,:].ravel()
        units = np.tile([1,1,1,REFERENCE.D,REFERENCE.D,REFERENCE.D],2)
        amplitude = max(abs(ends/units))
        residual = abs(eb.scalar_joint_residuals(ends/amplitude,beta,joint)/assembled.row_units)
        hat = reactions.ravel()/assembled.reaction_scales
        null = float(np.linalg.norm(assembled.dimensionless@hat)/(np.linalg.norm(assembled.dimensionless)*np.linalg.norm(hat)))
        delta = float(states[0,-1,2]-states[1,-1,2])
        failures = []
        if max(residual)>1e-9 or max(residual[:2])>1e-10:
            failures.append('RECONSTRUCTED_PHYSICAL_GATE')
        if null>1e-9:
            failures.append('RECONSTRUCTED_NULL_GATE')
        if eta is not None and np.linalg.norm(mirrored-eta*common)>1e-6:
            failures.append('SYMMETRY_UNRESOLVED')
        possibilities.append(dict(states=states,reactions=reactions,vector=vector,common=common,
            M=float(np.vdot(vector,vector).real),mass_before_normalization=mass,
            symmetry_class=eta,reflection_weight_plus=(1+reflection)/2,
            psi1=float(states[0,-1,2]),psi2=float(states[1,-1,2]),Delta_psi=delta,
            s=REFERENCE.D*delta**2/omega**2,
            physical_residuals=residual.tolist(),null_residual=null,sigma_ratio=endpoint['sigma_ratio'],
            detected_nullity=endpoint['nullity'],failures=failures,
            reconstruction_calls={key:sum(v[1][key] for v in values) for key in ('analytic','expm')},
            transfer_check_errors=[error for v in values for error in v[2]]))
    if not possibilities:
        raise ValueError('NO_CLASS_IN_NULLSPACE')
    return min(possibilities,key=lambda v:max(v['physical_residuals']))


def match(left,right,*,physical=False,restrict_symmetry=False):
    """Global assignment; an old symmetry label never constrains unequal arms."""
    key = 'vector' if physical else 'common'
    if restrict_symmetry and any(r['symmetry_class'] is None for r in left+right):
        raise ValueError('cannot restrict broken symmetry')
    classes = ([r['symmetry_class'] for r in left],[r['symmetry_class'] for r in right]) if restrict_symmetry else (None,None)
    indices,mac,margins = modes.assign([r[key] for r in left],[r[key] for r in right],*classes)
    return indices,mac[np.arange(len(left)),indices],margins
