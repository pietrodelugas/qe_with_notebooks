"""
ph_input.py — object-oriented layer for the QE phonon workflow.

Classes
-------
PhInputph        &INPUTPH namelist object  (ph.x)
  Modes (mutually exclusive, set via class-methods):
    PhInputph.single_q(prefix, qpoint, ...)   — one explicit q-vector
    PhInputph.dispersion(prefix, nq1,nq2,nq3) — uniform ldisp grid
    PhInputph.qplot(prefix, qpoints, ...)      — explicit list of q-points

  Sub-workflow activators (fluent API):
    .with_dielectric()    — epsil + zeu (Born charges)
    .with_raman()         — lraman (+ elop optional)
    .with_eph(method)     — electron_phonon=method + el_ph_sigma etc.
    .with_ahc(...)        — electron_phonon='ahc' + ahc_* parameters
    .with_recover()       — recover=.true.
    .parallel(start_irr, last_irr, start_q, last_q)  — job splitting

Q2rInput         &INPUT namelist object  (q2r.x)
  Q2rInput(fildyn, flfrc, zasr='crystal')

MatdynInput      &INPUT namelist object  (matdyn.x)
  Modes:
    MatdynInput.dispersion(flfrc, qpath, ...)  — phonon band structure
    MatdynInput.dos(flfrc, nk1,nk2,nk3, ...)  — phonon DOS

PhononWorkflow   top-level container: ph.x → q2r.x → matdyn.x
  .write(directory)  — writes ph.in, q2r.in, matdyn.in
  .run_commands()    — returns the shell commands to execute the workflow

QPointPath       helper for building q-point paths (band form or explicit)
"""

from __future__ import annotations
import copy
import os
from typing import Any

from ph_namelists import PH_INPUTPH, Q2R_INPUT, MATDYN_INPUT, DYNMAT_INPUT

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt(v: Any) -> str:
    """Format a Python value as a Fortran literal."""
    if isinstance(v, bool):
        return '.true.' if v else '.false.'
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, (list, tuple)):
        return ', '.join(_fmt(x) for x in v)
    return str(v)


def _validate(ref: dict, key: str, value: Any) -> None:
    if key not in ref:
        raise KeyError(f"Unknown parameter '{key}' — not in reference dict.")
    valid = ref[key].get('valid', [])
    if valid and value not in valid:
        raise ValueError(f"'{key}' = {value!r} not in valid choices: {valid}")


# ---------------------------------------------------------------------------
# Base namelist
# ---------------------------------------------------------------------------

class _Namelist:
    _ref:  dict = {}
    _name: str  = ''

    def __init__(self, **kwargs):
        self._params: dict = {}
        self.update(**kwargs)

    def update(self, **kwargs) -> '_Namelist':
        for k, v in kwargs.items():
            _validate(self._ref, k, v)
            self._params[k] = v
        return self

    def set(self, key: str, value: Any) -> '_Namelist':
        return self.update(**{key: value})

    def get(self, key: str) -> Any:
        if key in self._params:
            return self._params[key]
        return self._ref.get(key, {}).get('default')

    def to_dict(self) -> dict:
        return copy.deepcopy(self._params)

    def to_string(self) -> str:
        lines = [f'&{self._name}']
        for k, v in self._params.items():
            lines.append(f'  {k} = {_fmt(v)}')
        lines.append('/')
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._params})'


# ============================================================================
# PhInputph  — &INPUTPH
# ============================================================================

class PhInputph(_Namelist):
    """
    &INPUTPH namelist for ph.x.

    Do not instantiate directly — use one of the three class-method constructors:

        PhInputph.single_q(prefix, qpoint, ...)
        PhInputph.dispersion(prefix, nq1, nq2, nq3, ...)
        PhInputph.qplot(prefix, qpoints, ...)

    Then chain sub-workflow activators:

        ph = (PhInputph.dispersion('silicon', 4, 4, 4, fildyn='si.dyn')
                        .with_dielectric()
                        .with_recover())

    Parameters stored in self._params are written to the &INPUTPH block.
    The q-point specification (for single_q and qplot modes) is stored
    separately in self._qcard and appended after the namelist.
    """

    _ref  = PH_INPUTPH
    _name = 'INPUTPH'

    # --- internal init -------------------------------------------------------

    def __init__(self, **kwargs):
        self._params: dict = {}
        self._qcard: str   = ''      # text appended after the namelist
        self._mode: str    = ''      # 'single_q' | 'dispersion' | 'qplot'
        self.update(**kwargs)

    # =========================================================================
    # Class-method constructors — three mutually exclusive modes
    # =========================================================================

    @classmethod
    def single_q(cls,
                 prefix:   str,
                 qpoint:   tuple | list,
                 fildyn:   str   = 'matdyn',
                 outdir:   str   = './',
                 tr2_ph:   float = 1e-12,
                 **kwargs) -> 'PhInputph':
        """
        Compute phonons at a single explicit q-vector.

        Parameters
        ----------
        prefix  : str   — must match the pw.x prefix
        qpoint  : (qx, qy, qz) in units of 2π/a (tpiba)
        fildyn  : str   — dynamical matrix output file name
        outdir  : str
        tr2_ph  : float — SCF convergence threshold
        **kwargs        — any other &INPUTPH parameters

        Example
        -------
        >>> ph = PhInputph.single_q('silicon', (0, 0, 0), fildyn='si.dynG')
        """
        obj = cls.__new__(cls)
        obj._params = {}
        obj._mode   = 'single_q'

        q = list(qpoint)
        if len(q) != 3:
            raise ValueError('qpoint must be a 3-vector (qx, qy, qz)')

        obj.update(
            prefix=prefix,
            outdir=outdir,
            fildyn=fildyn,
            tr2_ph=tr2_ph,
            **kwargs,
        )
        # q-card: one q-point, weight 1
        obj._qcard = f'{q[0]:12.6f}  {q[1]:12.6f}  {q[2]:12.6f}  1\n'
        return obj

    @classmethod
    def dispersion(cls,
                   prefix: str,
                   nq1:    int,
                   nq2:    int,
                   nq3:    int,
                   fildyn: str   = 'matdyn',
                   outdir: str   = './',
                   tr2_ph: float = 1e-12,
                   lshift_q: bool = False,
                   **kwargs) -> 'PhInputph':
        """
        Compute phonons on a uniform nq1×nq2×nq3 Monkhorst-Pack q-grid.

        Parameters
        ----------
        prefix     : str
        nq1,nq2,nq3 : int — q-mesh dimensions
        fildyn     : str  — root name for dynamical matrix files
                            (files will be <fildyn>1, <fildyn>2, …)
        outdir     : str
        tr2_ph     : float
        lshift_q   : bool — shift grid to avoid Γ
        **kwargs          — any other &INPUTPH parameters

        Example
        -------
        >>> ph = PhInputph.dispersion('silicon', 4, 4, 4, fildyn='si.dyn')
        """
        obj = cls.__new__(cls)
        obj._params = {}
        obj._mode   = 'dispersion'
        obj._qcard  = ''

        obj.update(
            prefix=prefix,
            outdir=outdir,
            fildyn=fildyn,
            tr2_ph=tr2_ph,
            ldisp=True,
            nq1=nq1,
            nq2=nq2,
            nq3=nq3,
            **({'lshift_q': True} if lshift_q else {}),
            **kwargs,
        )
        return obj

    @classmethod
    def qplot(cls,
              prefix:  str,
              qpoints: list,
              fildyn:  str   = 'matdyn',
              outdir:  str   = './',
              tr2_ph:  float = 1e-12,
              q_in_cryst_coord: bool = False,
              q_in_band_form:   bool = False,
              **kwargs) -> 'PhInputph':
        """
        Compute phonons at an explicit list of q-points.

        Parameters
        ----------
        prefix   : str
        qpoints  : list of (qx, qy, qz) or (qx, qy, qz, npoints_to_next)
                   npoints_to_next is used with q_in_band_form=True.
                   If not provided it defaults to 1.
        fildyn   : str
        outdir   : str
        tr2_ph   : float
        q_in_cryst_coord : bool — use crystallographic coordinates
        q_in_band_form   : bool — band-path form (each row = endpoint + npoints)
        **kwargs

        Example
        -------
        >>> path = [(0,0,0,10), (0.5,0,0.5,10), (0.5,0.25,0.75,1)]
        >>> ph = PhInputph.qplot('silicon', path, q_in_band_form=True,
        ...                      q_in_cryst_coord=True)
        """
        obj = cls.__new__(cls)
        obj._params = {}
        obj._mode   = 'qplot'

        kw = dict(
            prefix=prefix,
            outdir=outdir,
            fildyn=fildyn,
            tr2_ph=tr2_ph,
            qplot=True,
        )
        if q_in_cryst_coord:
            kw['q_in_cryst_coord'] = True
        if q_in_band_form:
            kw['q_in_band_form'] = True
        kw.update(kwargs)
        obj.update(**kw)

        # Build the q-card block
        lines = [str(len(qpoints))]
        for q in qpoints:
            q = list(q)
            if len(q) == 3:
                q.append(1)
            lines.append(f'  {q[0]:12.8f}  {q[1]:12.8f}  {q[2]:12.8f}  {q[3]}')
        obj._qcard = '\n'.join(lines) + '\n'
        return obj

    # =========================================================================
    # Sub-workflow activators  (fluent — each returns self)
    # =========================================================================

    def with_dielectric(self,
                        epsil: bool = True,
                        zeu:   bool = True,
                        zue:   bool = False) -> 'PhInputph':
        """
        Activate dielectric tensor and Born effective charge calculation.

        Parameters
        ----------
        epsil : compute macroscopic dielectric tensor ε∞ (default True)
        zeu   : compute Born effective charges Z* via E-field perturbation (default True)
        zue   : also compute Z* via atomic displacement (cross-check, default False)

        Notes
        -----
        Only meaningful at q=Γ (single_q mode with qpoint=(0,0,0), or
        automatically included in dispersion mode which always visits Γ).
        """
        self._params['epsil'] = epsil
        self._params['zeu']   = zeu
        if zue:
            self._params['zue'] = True
        return self

    def with_raman(self,
                   elop: bool = False) -> 'PhInputph':
        """
        Activate non-resonant Raman tensor calculation.

        Requires epsil=.true. (calls with_dielectric() automatically if not set).

        Parameters
        ----------
        elop : also compute electro-optic tensor (default False)
        """
        if not self._params.get('epsil', False):
            self.with_dielectric()
        self._params['lraman'] = True
        if elop:
            self._params['elop'] = True
        return self

    def with_eph(self,
                 method:     str   = 'lambda_tetra',
                 el_ph_sigma: float = 0.02,
                 el_ph_nsigma: int  = 10,
                 fildvscf:   str   = 'dvscf',
                 nk1: int = 0, nk2: int = 0, nk3: int = 0) -> 'PhInputph':
        """
        Activate electron-phonon coupling calculation.

        Parameters
        ----------
        method       : one of 'simple', 'interpolated', 'lambda_tetra',
                       'gamma_tetra', 'yambo'
        el_ph_sigma  : Fermi-surface smearing (Ry)
        el_ph_nsigma : number of smearing values
        fildvscf     : root name for the δV_scf output (needed by EPW/Wannier)
        nk1,nk2,nk3  : dense k-mesh for tetrahedra methods (0 = use pw.x mesh)
        """
        valid_methods = ['simple', 'interpolated', 'lambda_tetra',
                         'gamma_tetra', 'yambo']
        if method not in valid_methods:
            raise ValueError(f"method must be one of {valid_methods}")
        self._params['electron_phonon'] = method
        self._params['el_ph_sigma']     = el_ph_sigma
        self._params['el_ph_nsigma']    = el_ph_nsigma
        self._params['fildvscf']        = fildvscf
        if nk1:
            self._params['nk1'] = nk1
        if nk2:
            self._params['nk2'] = nk2
        if nk3:
            self._params['nk3'] = nk3
        return self

    def with_ahc(self,
                 ahc_dir:      str = './ahc_dir',
                 ahc_nbnd:     int = 0,
                 ahc_nbndskip: int = 0,
                 fildvscf:     str = 'dvscf',
                 skip_upperfan:  bool = False,
                 ldoubledelta:   bool = False) -> 'PhInputph':
        """
        Activate Allen-Heine-Cardona electron-phonon self-energy calculation.

        Parameters
        ----------
        ahc_dir      : output directory for AHC files
        ahc_nbnd     : bands included in self-energy sum (0 = all)
        ahc_nbndskip : lowest bands to exclude
        fildvscf     : δV_scf file root name
        skip_upperfan : omit upper Fan term
        ldoubledelta  : use double-delta approximation
        """
        self._params['electron_phonon'] = 'ahc'
        self._params['ahc_dir']         = ahc_dir
        self._params['ahc_nbnd']        = ahc_nbnd
        self._params['ahc_nbndskip']    = ahc_nbndskip
        self._params['fildvscf']        = fildvscf
        if skip_upperfan:
            self._params['skip_upperfan'] = True
        if ldoubledelta:
            self._params['ldoubledelta'] = True
        return self

    def with_recover(self) -> 'PhInputph':
        """Set recover=.true. to restart from a previous incomplete run."""
        self._params['recover'] = True
        return self

    def parallel(self,
                 start_irr: int = None,
                 last_irr:  int = None,
                 start_q:   int = None,
                 last_q:    int = None) -> 'PhInputph':
        """
        Restrict this run to a subset of irreducible representations / q-points,
        for splitting a large calculation across multiple jobs.

        Parameters
        ----------
        start_irr, last_irr : irrep range (applies at each q-point)
        start_q,   last_q   : q-point range (dispersion mode only)
        """
        if start_irr is not None:
            self._params['start_irr'] = start_irr
        if last_irr is not None:
            self._params['last_irr'] = last_irr
        if start_q is not None:
            if self._mode != 'dispersion':
                raise ValueError("start_q / last_q only meaningful in dispersion mode.")
            self._params['start_q'] = start_q
        if last_q is not None:
            if self._mode != 'dispersion':
                raise ValueError("start_q / last_q only meaningful in dispersion mode.")
            self._params['last_q'] = last_q
        return self

    def only_init(self) -> 'PhInputph':
        """Set only_init=.true.: initialise data structures then stop (parallel prep)."""
        self._params['only_init'] = True
        return self

    # =========================================================================
    # Rendering
    # =========================================================================

    def to_string(self) -> str:
        lines = ['&INPUTPH']
        for k, v in self._params.items():
            if isinstance(v, dict):
                # dvscf_star / drho_star: render as separate logical variables
                # QE reads them as namelist group: dvscf_star%open = .true. etc.
                for sub_k, sub_v in v.items():
                    lines.append(f'  {k}%{sub_k} = {_fmt(sub_v)}')
            else:
                lines.append(f'  {k} = {_fmt(v)}')
        lines.append('/')
        if self._qcard:
            lines.append(self._qcard.rstrip())
        return '\n'.join(lines)

    @property
    def mode(self) -> str:
        return self._mode

    def __repr__(self) -> str:
        return (f'PhInputph(mode={self._mode!r}, '
                f"prefix={self._params.get('prefix','?')!r}, "
                f"fildyn={self._params.get('fildyn','?')!r})")


# ============================================================================
# QPointPath  — helper for building q-point paths for matdyn.x
# ============================================================================

class QPointPath:
    """
    Build a q-point path for matdyn.x in band form or as explicit points.

    Band form (recommended):
        path = QPointPath(cryst=True)
        path.add(0,   0,   0,   name='Γ', npoints=40)
        path.add(0.5, 0,   0.5, name='X', npoints=40)
        path.add(0.5, 0.25,0.75,name='W', npoints=1)   # last point: npoints=1

    Explicit list:
        path = QPointPath(cryst=True, band_form=False)
        path.add(0, 0, 0)
        path.add(0.1, 0, 0.1)
        ...
    """

    def __init__(self, cryst: bool = True, band_form: bool = True):
        """
        Parameters
        ----------
        cryst     : use crystallographic (reciprocal lattice) coordinates
        band_form : True = high-symmetry endpoint + npoints per segment
                    False = explicit q-point list
        """
        self.cryst     = cryst
        self.band_form = band_form
        self._points: list = []        # (qx, qy, qz, npoints, name)
        self._labels: dict = {}        # for documentation only

    def add(self,
            qx: float, qy: float, qz: float,
            npoints: int = 10,
            name: str = '') -> 'QPointPath':
        """
        Add a q-point to the path.

        Parameters
        ----------
        qx, qy, qz : coordinates
        npoints     : number of points from this q to the next segment
                      (use 1 for the final point in band_form)
        name        : optional label (printed as a comment)
        """
        self._points.append((qx, qy, qz, npoints, name))
        return self

    def to_string(self) -> str:
        """Render the q-point card block for matdyn.x input."""
        lines = [str(len(self._points))]
        for qx, qy, qz, npoints, name in self._points:
            comment = f'  ! {name}' if name else ''
            if self.band_form:
                lines.append(f'  {qx:10.6f}  {qy:10.6f}  {qz:10.6f}  {npoints}{comment}')
            else:
                lines.append(f'  {qx:10.6f}  {qy:10.6f}  {qz:10.6f}{comment}')
        return '\n'.join(lines)

    def __len__(self) -> int:
        return len(self._points)

    def __repr__(self) -> str:
        return (f'QPointPath(cryst={self.cryst}, band_form={self.band_form}, '
                f'npoints={len(self._points)})')


# ============================================================================
# Q2rInput  — &INPUT for q2r.x
# ============================================================================

class Q2rInput(_Namelist):
    """
    &INPUT namelist for q2r.x.

    Usage
    -----
    q2r = Q2rInput(fildyn='si.dyn', flfrc='si.fc', zasr='crystal')
    print(q2r.to_string())
    """

    _ref  = Q2R_INPUT
    _name = 'INPUT'

    def __init__(self,
                 fildyn: str,
                 flfrc:  str,
                 zasr:   str  = 'crystal',
                 loto_2d: bool = False,
                 **kwargs):
        """
        Parameters
        ----------
        fildyn  : root name of ph.x dynamical matrix files (required)
        flfrc   : output interatomic force constants file (required)
        zasr    : acoustic sum rule: 'no','simple','crystal','one-dim','zero-dim'
        loto_2d : 2D LO-TO splitting (polar 2D materials)
        """
        super().__init__(fildyn=fildyn, flfrc=flfrc, zasr=zasr, **kwargs)
        if loto_2d:
            self._params['loto_2d'] = True

    def __repr__(self) -> str:
        return (f"Q2rInput(fildyn={self._params.get('fildyn')!r}, "
                f"flfrc={self._params.get('flfrc')!r}, "
                f"zasr={self._params.get('zasr')!r})")


# ============================================================================
# MatdynInput  — &INPUT for matdyn.x
# ============================================================================

class MatdynInput(_Namelist):
    """
    &INPUT namelist for matdyn.x.

    Use one of the two class-method constructors:

        MatdynInput.dispersion(flfrc, qpath, ...)   — phonon band structure
        MatdynInput.dos(flfrc, nk1, nk2, nk3, ...)  — phonon DOS

    The q-point path (for dispersion mode) is stored in self._qpath
    (a QPointPath object) and appended after the namelist.
    """

    _ref  = MATDYN_INPUT
    _name = 'INPUT'

    def __init__(self, **kwargs):
        self._params: dict  = {}
        self._qpath: QPointPath | None = None
        self.update(**kwargs)

    # =========================================================================
    # Class-method constructors
    # =========================================================================

    @classmethod
    def dispersion(cls,
                   flfrc:   str,
                   qpath:   QPointPath,
                   asr:     str  = 'crystal',
                   flfrq:   str  = 'matdyn.freq',
                   flvec:   str  = 'matdyn.modes',
                   loto_disable: bool = False,
                   eigen_similarity: bool = True,
                   **kwargs) -> 'MatdynInput':
        """
        Compute phonon dispersion along a q-point path.

        Parameters
        ----------
        flfrc    : interatomic force constants file from q2r.x (required)
        qpath    : QPointPath object defining the band path
        asr      : acoustic sum rule ('no','simple','crystal',…)
        flfrq    : output file for phonon frequencies
        flvec    : output file for phonon eigenvectors
        loto_disable     : disable LO-TO splitting (non-polar materials)
        eigen_similarity : sort branches by eigenvector continuity
        **kwargs         : any other &INPUT parameters

        Example
        -------
        >>> path = (QPointPath(cryst=True)
        ...         .add(0,0,0,       name='Γ', npoints=40)
        ...         .add(0.5,0,0.5,   name='X', npoints=40)
        ...         .add(0.5,0.25,0.75,name='W',npoints=1))
        >>> md = MatdynInput.dispersion('si.fc', path)
        """
        if not isinstance(qpath, QPointPath):
            raise TypeError('qpath must be a QPointPath instance')

        obj = cls.__new__(cls)
        obj._params = {}
        obj._qpath  = qpath

        kw = dict(
            flfrc=flfrc,
            asr=asr,
            flfrq=flfrq,
            flvec=flvec,
            q_in_band_form=qpath.band_form,
            q_in_cryst_coord=qpath.cryst,
            eigen_similarity=eigen_similarity,
        )
        if loto_disable:
            kw['loto_disable'] = True
        kw.update(kwargs)
        obj.update(**kw)
        return obj

    @classmethod
    def dos(cls,
            flfrc:  str,
            nk1:    int,
            nk2:    int,
            nk3:    int,
            asr:    str   = 'crystal',
            fldos:  str   = 'matdyn.dos',
            deltaE: float = 1.0,
            degauss: float = 0.0,
            nosym:  bool  = False,
            **kwargs) -> 'MatdynInput':
        """
        Compute phonon density of states on a q-mesh.

        Parameters
        ----------
        flfrc        : interatomic force constants file from q2r.x (required)
        nk1,nk2,nk3  : q-mesh dimensions (should be denser than the ph.x grid)
        asr          : acoustic sum rule
        fldos        : output DOS file
        deltaE       : frequency bin width (cm⁻¹)
        degauss      : Gaussian smearing for DOS (cm⁻¹)
        nosym        : disable symmetry reduction of q-mesh
        **kwargs

        Example
        -------
        >>> md = MatdynInput.dos('si.fc', 16, 16, 16, deltaE=2.0)
        """
        obj = cls.__new__(cls)
        obj._params = {}
        obj._qpath  = None

        kw = dict(
            flfrc=flfrc,
            asr=asr,
            dos=True,
            nk1=nk1, nk2=nk2, nk3=nk3,
            fldos=fldos,
            deltaE=deltaE,
        )
        if degauss:
            kw['degauss'] = degauss
        if nosym:
            kw['nosym'] = True
        kw.update(kwargs)
        obj.update(**kw)
        return obj

    # =========================================================================
    # Rendering
    # =========================================================================

    def to_string(self) -> str:
        lines = ['&INPUT']
        for k, v in self._params.items():
            lines.append(f'  {k} = {_fmt(v)}')
        lines.append('/')
        if self._qpath is not None:
            lines.append(self._qpath.to_string())
        return '\n'.join(lines)

    @property
    def mode(self) -> str:
        return 'dos' if self._params.get('dos', False) else 'dispersion'

    def __repr__(self) -> str:
        return (f"MatdynInput(mode={self.mode!r}, "
                f"flfrc={self._params.get('flfrc')!r})")


# ============================================================================
# PhononWorkflow  — top-level container
# ============================================================================

class PhononWorkflow:
    """
    Complete ph.x → q2r.x → matdyn.x workflow container.

    Usage
    -----
    ph = (PhInputph.dispersion('silicon', 4, 4, 4, fildyn='si.dyn')
                   .with_dielectric())

    q2r = Q2rInput(fildyn='si.dyn', flfrc='si.fc', zasr='crystal')

    path = (QPointPath(cryst=True)
            .add(0,   0,   0,    name='Γ', npoints=40)
            .add(0.5, 0,   0.5,  name='X', npoints=40)
            .add(0.5, 0.25,0.75, name='W', npoints=1))

    md_disp = MatdynInput.dispersion('si.fc', path)
    md_dos  = MatdynInput.dos('si.fc', 16, 16, 16)

    wf = PhononWorkflow(ph, q2r, matdyn_disp=md_disp, matdyn_dos=md_dos)
    wf.write('/content/phonon')
    for cmd in wf.run_commands():
        print(cmd)
    """

    def __init__(self,
                 ph:           PhInputph,
                 q2r:          Q2rInput       = None,
                 matdyn_disp:  MatdynInput    = None,
                 matdyn_dos:   MatdynInput    = None,
                 conda_env:    str            = 'qe_env'):
        """
        Parameters
        ----------
        ph           : PhInputph  (required)
        q2r          : Q2rInput   (optional; skip if you only want ph.x)
        matdyn_disp  : MatdynInput in dispersion mode (optional)
        matdyn_dos   : MatdynInput in DOS mode (optional)
        conda_env    : conda environment name used in run_commands()
        """
        self.ph          = ph
        self.q2r         = q2r
        self.matdyn_disp = matdyn_disp
        self.matdyn_dos  = matdyn_dos
        self.conda_env   = conda_env

        # Consistency checks
        if q2r is not None and ph._params.get('fildyn') != q2r._params.get('fildyn'):
            import warnings
            warnings.warn(
                f"ph.x fildyn={ph._params.get('fildyn')!r} does not match "
                f"q2r.x fildyn={q2r._params.get('fildyn')!r}"
            )
        if matdyn_disp is not None and matdyn_disp.mode != 'dispersion':
            raise ValueError("matdyn_disp must be a dispersion-mode MatdynInput.")
        if matdyn_dos is not None and matdyn_dos.mode != 'dos':
            raise ValueError("matdyn_dos must be a dos-mode MatdynInput.")

    def write(self, directory: str = '.') -> None:
        """Write all input files to the given directory."""
        os.makedirs(directory, exist_ok=True)

        files = {'ph.in': self.ph.to_string()}
        if self.q2r is not None:
            files['q2r.in'] = self.q2r.to_string()
        if self.matdyn_disp is not None:
            files['matdyn_disp.in'] = self.matdyn_disp.to_string()
        if self.matdyn_dos is not None:
            files['matdyn_dos.in'] = self.matdyn_dos.to_string()

        for fname, content in files.items():
            path = os.path.join(directory, fname)
            with open(path, 'w') as f:
                f.write(content + '\n')
            print(f'Written: {path}')

    def run_commands(self, nproc: int = 1) -> list[str]:
        """
        Return a list of shell commands to run the full workflow.

        Parameters
        ----------
        nproc : number of MPI processes (default 1 = serial)
        """
        def _cmd(exe, inp, nproc=nproc, env=self.conda_env):
            if nproc > 1:
                return (f'conda run -n {env} mpirun -np {nproc} '
                        f'{exe} -input {inp}')
            return f'conda run -n {env} {exe} -input {inp}'

        cmds = [_cmd('ph.x', 'ph.in')]
        if self.q2r is not None:
            cmds.append(_cmd('q2r.x', 'q2r.in', nproc=1))  # q2r is serial
        if self.matdyn_disp is not None:
            cmds.append(_cmd('matdyn.x', 'matdyn_disp.in', nproc=1))
        if self.matdyn_dos is not None:
            cmds.append(_cmd('matdyn.x', 'matdyn_dos.in', nproc=1))
        return cmds

    def to_string(self) -> str:
        """Return all input files concatenated with separators."""
        parts = [f'# ── ph.in ──\n{self.ph.to_string()}']
        if self.q2r is not None:
            parts.append(f'# ── q2r.in ──\n{self.q2r.to_string()}')
        if self.matdyn_disp is not None:
            parts.append(f'# ── matdyn_disp.in ──\n{self.matdyn_disp.to_string()}')
        if self.matdyn_dos is not None:
            parts.append(f'# ── matdyn_dos.in ──\n{self.matdyn_dos.to_string()}')
        return '\n\n'.join(parts)

    def __repr__(self) -> str:
        parts = ['PhononWorkflow(', f'  ph      = {self.ph!r}']
        if self.q2r:
            parts.append(f'  q2r     = {self.q2r!r}')
        if self.matdyn_disp:
            parts.append(f'  matdyn_disp = {self.matdyn_disp!r}')
        if self.matdyn_dos:
            parts.append(f'  matdyn_dos  = {self.matdyn_dos!r}')
        parts.append(')')
        return '\n'.join(parts)


# ============================================================================
# DynmatInput  — /input/ namelist for dynmat.x
# ============================================================================

class DynmatInput(_Namelist):
    """
    /input/ namelist for dynmat.x.

    Usage
    -----
    # Minimal — just read and diagonalise:
    dm = DynmatInput('si.dyn', asr='crystal')

    # With LO-TO splitting along [1 0 0]:
    dm = DynmatInput('si.dyn', asr='crystal', q=(0.1, 0, 0))

    # Compute ionic dielectric permittivity:
    dm = DynmatInput('si.dyn', asr='crystal', q=(0.1, 0, 0), lperm=True)

    Parameters stored in self._params are rendered to the /input/ block.
    The Fortran array variables q and amass are written as indexed elements:
        q(1) = …  q(2) = …  q(3) = …
        amass(1) = …  amass(2) = …  …
    """

    _ref  = DYNMAT_INPUT
    _name = 'input'

    def __init__(self,
                 fildyn: str,
                 asr:    str  = 'no',
                 q:      tuple | list | None = None,
                 lperm:  bool = False,
                 **kwargs):
        """
        Parameters
        ----------
        fildyn : str
            Dynamical matrix file (root name for ldisp runs; exact name for
            single-q runs).
        asr : str
            Acoustic sum rule: 'no', 'simple', 'crystal', 'one-dim', 'zero-dim'.
        q : (qx, qy, qz) | None
            LO-TO splitting direction in Cartesian 2π/a units.  If None or
            (0, 0, 0), no non-analytic correction is applied (TO frequencies
            only).  Written as q(1) = …, q(2) = …, q(3) = … in the namelist.
        lperm : bool
            Compute dielectric permittivity tensor (ε∞ + ionic correction).
            Requires ε∞ and Z* to be stored in fildyn.
        **kwargs
            Any other /input/ parameters (filout, filmol, filxsf, fileig,
            axis, lplasma, loto_2d, amass, remove_interaction_blocks).
        """
        super().__init__(fildyn=fildyn, asr=asr, **kwargs)
        if q is not None:
            q = list(q)
            if len(q) != 3:
                raise ValueError('q must be a 3-vector (qx, qy, qz)')
            self._params['q'] = q
        if lperm:
            self._params['lperm'] = True

    # =========================================================================
    # Rendering  — override to handle Fortran array-indexed variables
    # =========================================================================

    def to_string(self) -> str:
        lines = [f'&{self._name}']
        for k, v in self._params.items():
            if k == 'q':
                for i, qi in enumerate(v, 1):
                    lines.append(f'  q({i}) = {_fmt(qi)}')
            elif k == 'amass':
                for i, mi in enumerate(list(v), 1):
                    lines.append(f'  amass({i}) = {_fmt(mi)}')
            else:
                lines.append(f'  {k} = {_fmt(v)}')
        lines.append('/')
        return '\n'.join(lines)

    def __repr__(self) -> str:
        q = self._params.get('q')
        q_str = f'({", ".join(str(x) for x in q)})' if q else 'None'
        return (
            f"DynmatInput(fildyn={self._params.get('fildyn')!r}, "
            f"asr={self._params.get('asr')!r}, "
            f"q={q_str}, "
            f"lperm={self._params.get('lperm', False)})"
        )
