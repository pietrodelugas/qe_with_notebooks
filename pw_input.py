"""
pw_input.py — object-oriented layer for building pw.x input files.

Classes
-------
ControlNamelist   &CONTROL
SystemNamelist    &SYSTEM  (ibrav-aware; can be populated from an ASE Atoms object)
ElectronsNamelist &ELECTRONS
IonsNamelist      &IONS
CellNamelist      &CELL
AtomicSpeciesCard ATOMIC_SPECIES
AtomicPositionsCard ATOMIC_POSITIONS  (can be populated from an ASE Atoms object)
KPointsAutoCard   K_POINTS {automatic}  (ibrav-aware constructor)
PWInput           top-level container → renders a complete pw.x input file

Design principles
-----------------
* Every namelist class stores only parameters that differ from the QE default,
  so the rendered output is minimal and readable.
* set() / update() accept keyword arguments validated against the reference dict.
* to_string() renders the Fortran namelist block or card block.
* from_atoms(atoms, ...) class-methods handle the ASE→QE translation.
* KPointsAutoCard.__init__ accepts only the celldm parameters relevant to the
  chosen ibrav, raising TypeError for incompatible ones.
"""

from __future__ import annotations
import copy
import math
from typing import Any

# Reference dicts (schema: default, type, unit, description, valid)
from pw_namelists import (
    CONTROL  as _REF_CONTROL,
    SYSTEM   as _REF_SYSTEM,
    ELECTRONS as _REF_ELECTRONS,
    IONS     as _REF_IONS,
    CELL     as _REF_CELL,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_BOHR_TO_ANG = 0.529177210903
_ANG_TO_BOHR = 1.0 / _BOHR_TO_ANG


def _fmt_value(v: Any) -> str:
    """Format a Python value as a Fortran literal."""
    if isinstance(v, bool):
        return '.true.' if v else '.false.'
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, (list, tuple)):
        return ', '.join(_fmt_value(x) for x in v)
    return str(v)


def _validate(ref: dict, key: str, value: Any) -> None:
    if key not in ref:
        raise KeyError(f"Unknown parameter '{key}'")
    valid = ref[key].get('valid', [])
    if valid and value not in valid:
        raise ValueError(f"'{key}' = {value!r} not in valid choices: {valid}")


# ---------------------------------------------------------------------------
# Base namelist class
# ---------------------------------------------------------------------------

class _Namelist:
    """
    Base class for a Fortran &NAMELIST block.

    Subclasses set:
        _ref   : the reference dict (e.g. CONTROL)
        _name  : the Fortran namelist name (e.g. 'CONTROL')
    """
    _ref: dict = {}
    _name: str = ''

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
        return self._ref[key]['default']

    def to_dict(self) -> dict:
        """Return the current (non-default) parameters as a plain dict."""
        return copy.deepcopy(self._params)

    def to_string(self) -> str:
        lines = [f'&{self._name}']
        for k, v in self._params.items():
            lines.append(f'  {k} = {_fmt_value(v)}')
        lines.append('/')
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._params})'


# ---------------------------------------------------------------------------
# Concrete namelist classes
# ---------------------------------------------------------------------------

class ControlNamelist(_Namelist):
    _ref  = _REF_CONTROL
    _name = 'CONTROL'


class ElectronsNamelist(_Namelist):
    _ref  = _REF_ELECTRONS
    _name = 'ELECTRONS'


class IonsNamelist(_Namelist):
    _ref  = _REF_IONS
    _name = 'IONS'


class CellNamelist(_Namelist):
    _ref  = _REF_CELL
    _name = 'CELL'


# ---------------------------------------------------------------------------
# ibrav metadata
# ---------------------------------------------------------------------------

# For each ibrav: (label, needed celldm indices, needed angle cosines if any)
# celldm indices are 1-based as in QE documentation.
# 'free' means ibrav=0 → CELL_PARAMETERS card required, no celldm.
_IBRAV_INFO = {
    0:  dict(label='Free (CELL_PARAMETERS required)',
             celldm=[], angles=[]),
    1:  dict(label='Simple cubic (sc)',
             celldm=[1], angles=[]),
    2:  dict(label='Face-centred cubic (fcc)',
             celldm=[1], angles=[]),
    3:  dict(label='Body-centred cubic (bcc)',
             celldm=[1], angles=[]),
   -3:  dict(label='BCC (primitive, alternative setting)',
             celldm=[1], angles=[]),
    4:  dict(label='Hexagonal / trigonal',
             celldm=[1, 3], angles=[]),          # c/a = celldm(3)
    5:  dict(label='Trigonal (rhombohedral, 3-fold axis along (111))',
             celldm=[1, 4], angles=[]),           # cos(α) = celldm(4)
   -5:  dict(label='Trigonal (rhombohedral, alternative orientation)',
             celldm=[1, 4], angles=[]),
    6:  dict(label='Simple tetragonal (st)',
             celldm=[1, 3], angles=[]),           # c/a
    7:  dict(label='Body-centred tetragonal (bct)',
             celldm=[1, 3], angles=[]),
    8:  dict(label='Simple orthorhombic',
             celldm=[1, 2, 3], angles=[]),        # b/a, c/a
    9:  dict(label='Base-centred orthorhombic (C)',
             celldm=[1, 2, 3], angles=[]),
   -9:  dict(label='Base-centred orthorhombic (A, alternative)',
             celldm=[1, 2, 3], angles=[]),
   91:  dict(label='Base-centred orthorhombic (A)',
             celldm=[1, 2, 3], angles=[]),
   10:  dict(label='Face-centred orthorhombic',
             celldm=[1, 2, 3], angles=[]),
   11:  dict(label='Body-centred orthorhombic',
             celldm=[1, 2, 3], angles=[]),
   12:  dict(label='Simple monoclinic (unique axis c)',
             celldm=[1, 2, 3, 4], angles=[]),     # cos(γ) = celldm(4)
  -12:  dict(label='Simple monoclinic (unique axis b)',
             celldm=[1, 2, 3, 5], angles=[]),     # cos(β) = celldm(5)
   13:  dict(label='Base-centred monoclinic (unique axis c)',
             celldm=[1, 2, 3, 4], angles=[]),
  -13:  dict(label='Base-centred monoclinic (unique axis b)',
             celldm=[1, 2, 3, 5], angles=[]),
   14:  dict(label='Triclinic',
             celldm=[1, 2, 3, 4, 5, 6], angles=[]),
}

# Human-readable celldm index labels
_CELLDM_LABEL = {
    1: 'a (bohr)',
    2: 'b/a',
    3: 'c/a',
    4: 'cos(α)  [or cos(γ) for ibrav=12/-13, or cos(α) for ibrav=5/-5]',
    5: 'cos(β)  [monoclinic unique-b]  or  cos(β) [triclinic]',
    6: 'cos(γ)  [triclinic]',
}


def ibrav_info(ibrav: int) -> dict:
    """Return the metadata dict for a given ibrav value."""
    if ibrav not in _IBRAV_INFO:
        raise ValueError(f'ibrav={ibrav} is not a recognised QE Bravais lattice index.')
    return _IBRAV_INFO[ibrav]


def ibrav_celldm_params(ibrav: int) -> list[int]:
    """Return which celldm indices (1-based) are required for this ibrav."""
    return _IBRAV_INFO[ibrav]['celldm']


# ---------------------------------------------------------------------------
# SystemNamelist  (ibrav-aware)
# ---------------------------------------------------------------------------

class SystemNamelist(_Namelist):
    """
    &SYSTEM namelist with ibrav-awareness.

    Parameters
    ----------
    ibrav : int
        Bravais lattice index (required).  Sets which celldm_* kwargs are
        accepted; others raise TypeError.
    celldm_1 … celldm_6 : float, optional
        Lattice constants.  Only those compatible with ibrav are accepted.
    **kwargs :
        Any other &SYSTEM parameter.

    Class methods
    -------------
    from_atoms(atoms, ibrav=0, ecutwfc=..., **kwargs)
        Populate ibrav, nat, ntyp, celldm(1), and optionally celldm(2/3/4…)
        from an ASE Atoms object.  ibrav=0 stores nothing (use CELL_PARAMETERS).
    """
    _ref  = _REF_SYSTEM
    _name = 'SYSTEM'

    def __init__(self, ibrav: int, **kwargs):
        self._params: dict = {}
        if ibrav not in _IBRAV_INFO:
            raise ValueError(f'ibrav={ibrav} is not valid. '
                             f'Valid values: {sorted(_IBRAV_INFO)}')
        self._params['ibrav'] = ibrav
        self._ibrav = ibrav

        # Separate celldm_N kwargs from the rest
        celldm_kwargs = {k: v for k, v in kwargs.items() if k.startswith('celldm_')}
        other_kwargs  = {k: v for k, v in kwargs.items() if not k.startswith('celldm_')}

        self._set_celldm(**celldm_kwargs)
        if other_kwargs:
            self.update(**other_kwargs)

    # ---- celldm handling ---------------------------------------------------

    def _celldm_list(self) -> list:
        """Return the current celldm list (expanding sparse entries)."""
        cd = self._params.get('celldm', [None] * 6)
        if len(cd) < 6:
            cd = list(cd) + [None] * (6 - len(cd))
        return cd

    def _set_celldm(self, **kwargs) -> None:
        """
        Accept celldm_1 … celldm_6 kwargs, validate against ibrav,
        and store as celldm list.
        """
        allowed = ibrav_celldm_params(self._ibrav)
        cd = self._celldm_list()
        for k, v in kwargs.items():
            if not k.startswith('celldm_'):
                raise TypeError(f"Expected 'celldm_N' keyword, got '{k}'")
            idx = int(k.split('_')[1])  # 1-based
            if idx not in allowed:
                allowed_str = ', '.join(f'celldm_{i}' for i in allowed)
                raise TypeError(
                    f"celldm_{idx} is not compatible with ibrav={self._ibrav} "
                    f"({_IBRAV_INFO[self._ibrav]['label']}). "
                    f"Allowed: {allowed_str}"
                )
            cd[idx - 1] = float(v)
        # Trim trailing Nones
        while cd and cd[-1] is None:
            cd.pop()
        if cd:
            self._params['celldm'] = cd

    def set_celldm(self, **kwargs) -> 'SystemNamelist':
        """Public method: set celldm_1 … celldm_6 with ibrav validation."""
        self._set_celldm(**kwargs)
        return self

    @property
    def ibrav(self) -> int:
        return self._ibrav

    @ibrav.setter
    def ibrav(self, value: int) -> None:
        if value not in _IBRAV_INFO:
            raise ValueError(f'ibrav={value} is not valid.')
        # Check existing celldm is still compatible
        new_allowed = ibrav_celldm_params(value)
        cd = self._celldm_list()
        for i, v in enumerate(cd):
            if v is not None and (i + 1) not in new_allowed:
                raise ValueError(
                    f"Cannot change ibrav to {value}: celldm_{i+1} is set but "
                    f"not compatible with ibrav={value}. Clear celldm first."
                )
        self._ibrav = value
        self._params['ibrav'] = value

    def celldm_info(self) -> str:
        """Print which celldm parameters are needed for the current ibrav."""
        info = _IBRAV_INFO[self._ibrav]
        lines = [f"ibrav={self._ibrav}  {info['label']}",
                 "Required celldm parameters:"]
        for idx in info['celldm']:
            lines.append(f"  celldm_{idx}  —  {_CELLDM_LABEL[idx]}")
        return '\n'.join(lines)

    # ---- from ASE Atoms ----------------------------------------------------

    @classmethod
    def from_atoms(cls,
                   atoms,
                   ibrav: int = 0,
                   ecutwfc: float = None,
                   ecutrho: float = None,
                   **kwargs) -> 'SystemNamelist':
        """
        Build a SystemNamelist from an ASE Atoms object.

        Parameters
        ----------
        atoms   : ase.Atoms
        ibrav   : int
            0  → free cell (no celldm set; CELL_PARAMETERS card needed).
            2  → FCC, only celldm(1)=a extracted.
            4  → hexagonal, celldm(1)=a and celldm(3)=c/a extracted.
            … any ibrav: a=celldm(1) is always extracted from atoms.cell;
            higher celldm are extracted where the geometry uniquely defines them.
        ecutwfc : float (Ry), optional
        ecutrho : float (Ry), optional
        **kwargs : other &SYSTEM parameters
        """
        from ase.data import atomic_numbers  # noqa: F401 — just to confirm ASE present

        cell = atoms.cell  # ASE Cell object (Angstrom)
        symbols = atoms.get_chemical_symbols()
        species = list(dict.fromkeys(symbols))  # ordered unique species

        kw = {}
        kw['nat']  = len(atoms)
        kw['ntyp'] = len(species)
        if ecutwfc is not None:
            kw['ecutwfc'] = ecutwfc
        if ecutrho is not None:
            kw['ecutrho'] = ecutrho
        kw.update(kwargs)

        obj = cls(ibrav=ibrav, **kw)

        if ibrav == 0:
            # No celldm — caller must provide a CellParametersCard
            return obj

        # Extract lattice parameters from the ASE cell (in bohr)
        lengths = [v * _ANG_TO_BOHR for v in cell.lengths()]  # a, b, c in bohr
        angles  = cell.angles()                                 # α, β, γ in degrees
        a, b, c = lengths
        alpha, beta, gamma = angles

        allowed = ibrav_celldm_params(ibrav)
        cd: dict[int, float] = {}

        if 1 in allowed:
            cd[1] = a

        # Ratios — only set if the ibrav actually uses them
        if 2 in allowed and a > 0:
            cd[2] = b / a
        if 3 in allowed and a > 0:
            cd[3] = c / a

        # Angle cosines
        if ibrav in (5, -5) and 4 in allowed:
            cd[4] = math.cos(math.radians(alpha))
        if ibrav == 12 and 4 in allowed:
            cd[4] = math.cos(math.radians(gamma))
        if ibrav == -12 and 5 in allowed:
            cd[5] = math.cos(math.radians(beta))
        if ibrav in (13, -13):
            if 4 in allowed:
                cd[4] = math.cos(math.radians(gamma))
        if ibrav == 14:
            if 4 in allowed:
                cd[4] = math.cos(math.radians(alpha))
            if 5 in allowed:
                cd[5] = math.cos(math.radians(beta))
            if 6 in allowed:
                cd[6] = math.cos(math.radians(gamma))

        celldm_kwargs = {f'celldm_{i}': v for i, v in cd.items()}
        obj._set_celldm(**celldm_kwargs)
        return obj

    # ---- rendering ---------------------------------------------------------

    def to_string(self) -> str:
        lines = ['&SYSTEM']
        for k, v in self._params.items():
            if k == 'celldm':
                for i, val in enumerate(v, start=1):
                    if val is not None:
                        lines.append(f'  celldm({i}) = {_fmt_value(val)}')
            else:
                lines.append(f'  {k} = {_fmt_value(v)}')
        lines.append('/')
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Cards
# ---------------------------------------------------------------------------

class AtomicSpeciesCard:
    """
    ATOMIC_SPECIES card.

    Usage
    -----
    card = AtomicSpeciesCard()
    card.add('Si', 28.0855, 'Si.pbe-n-kjpaw_psl.1.0.0.UPF')
    card.add('O',  15.9994, 'O.pbe-n-kjpaw_psl.1.0.0.UPF')
    """

    def __init__(self):
        self._species: list[tuple] = []  # (label, mass, pseudo)

    def add(self, label: str, mass: float, pseudo: str) -> 'AtomicSpeciesCard':
        self._species.append((label, float(mass), pseudo))
        return self

    @classmethod
    def from_atoms(cls, atoms, pseudos: dict) -> 'AtomicSpeciesCard':
        """
        Build from an ASE Atoms object.

        Parameters
        ----------
        atoms   : ase.Atoms
        pseudos : dict  {symbol: (mass_amu, pseudo_filename)}
                  or    {symbol: pseudo_filename}  (mass from ASE data)
        """
        from ase.data import atomic_masses, atomic_numbers
        obj = cls()
        seen = []
        for sym in atoms.get_chemical_symbols():
            if sym in seen:
                continue
            seen.append(sym)
            val = pseudos[sym]
            if isinstance(val, (list, tuple)):
                mass, pseudo = val[0], val[1]
            else:
                pseudo = val
                Z = atomic_numbers[sym]
                mass = atomic_masses[Z]
            obj.add(sym, mass, pseudo)
        return obj

    def to_string(self) -> str:
        lines = ['ATOMIC_SPECIES']
        for label, mass, pseudo in self._species:
            lines.append(f'  {label:<6s}  {mass:10.4f}  {pseudo}')
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return f'AtomicSpeciesCard({self._species})'


class AtomicPositionsCard:
    """
    ATOMIC_POSITIONS card.

    Usage
    -----
    card = AtomicPositionsCard(units='crystal')
    card.add('Si', 0.0, 0.0, 0.0)
    card.add('Si', 0.25, 0.25, 0.25)

    # Or from ASE Atoms:
    card = AtomicPositionsCard.from_atoms(atoms, units='angstrom')
    """

    VALID_UNITS = ('alat', 'bohr', 'angstrom', 'crystal', 'crystal_sg')

    def __init__(self, units: str = 'angstrom'):
        units = units.lower()
        if units not in self.VALID_UNITS:
            raise ValueError(f"units must be one of {self.VALID_UNITS}")
        self.units = units
        self._positions: list[tuple] = []  # (label, x, y, z, if1, if2, if3)

    def add(self,
            label: str,
            x: float, y: float, z: float,
            fix: tuple[int, int, int] = (1, 1, 1)) -> 'AtomicPositionsCard':
        """
        Add an atom.

        Parameters
        ----------
        label : species label
        x, y, z : coordinates in the chosen units
        fix : (if1, if2, if3) — 1=free, 0=fixed. Default all free.
        """
        self._positions.append((label, x, y, z, fix[0], fix[1], fix[2]))
        return self

    @classmethod
    def from_atoms(cls, atoms, units: str = 'angstrom',
                   constraints: dict = None) -> 'AtomicPositionsCard':
        """
        Build from an ASE Atoms object.

        Parameters
        ----------
        atoms       : ase.Atoms
        units       : 'angstrom' (default), 'crystal', 'bohr', or 'alat'
        constraints : dict {atom_index: (if1,if2,if3)} for fixed atoms.
                      Default: all atoms free.
        """
        import numpy as np
        if constraints is None:
            constraints = {}

        obj = cls(units=units)
        symbols = atoms.get_chemical_symbols()

        if units == 'angstrom':
            coords = atoms.get_positions()          # Å, Cartesian
        elif units == 'bohr':
            coords = atoms.get_positions() * _ANG_TO_BOHR
        elif units == 'crystal':
            coords = atoms.get_scaled_positions()   # fractional
        elif units == 'alat':
            # Cartesian in units of a = |a1|
            a = atoms.cell.lengths()[0] * _ANG_TO_BOHR  # bohr
            coords = atoms.get_positions() * _ANG_TO_BOHR / a
        else:
            raise ValueError(f"units '{units}' not supported in from_atoms")

        for i, (sym, pos) in enumerate(zip(symbols, coords)):
            fix = constraints.get(i, (1, 1, 1))
            obj.add(sym, pos[0], pos[1], pos[2], fix=fix)
        return obj

    def to_string(self) -> str:
        lines = [f'ATOMIC_POSITIONS {{{self.units}}}']
        for label, x, y, z, if1, if2, if3 in self._positions:
            fixed = f'  {if1} {if2} {if3}' if (if1, if2, if3) != (1, 1, 1) else ''
            lines.append(f'  {label:<6s}  {x:16.10f}  {y:16.10f}  {z:16.10f}{fixed}')
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return f'AtomicPositionsCard(units={self.units!r}, n={len(self._positions)})'


# ---------------------------------------------------------------------------
# K_POINTS {automatic}  — ibrav-aware constructor
# ---------------------------------------------------------------------------

# For each ibrav, the mesh has a natural symmetry constraint.
# This table encodes which mesh dimensions must be equal.
# Format: list of tuples of dimension indices (0=nk1, 1=nk2, 2=nk3) that must match.
_IBRAV_KMESH_CONSTRAINTS = {
    0:  [],                      # free cell — no constraint enforced
    1:  [(0, 1, 2)],             # sc  → nk1 = nk2 = nk3
    2:  [(0, 1, 2)],             # fcc → nk1 = nk2 = nk3
    3:  [(0, 1, 2)],             # bcc → nk1 = nk2 = nk3
   -3:  [(0, 1, 2)],
    4:  [(0, 1)],                # hex/trig → nk1 = nk2, nk3 independent
    5:  [(0, 1, 2)],             # rhomb → nk1 = nk2 = nk3
   -5:  [(0, 1, 2)],
    6:  [(0, 1)],                # st → nk1 = nk2, nk3 independent
    7:  [(0, 1)],                # bct → nk1 = nk2, nk3 independent
    8:  [],                      # ortho → all independent
    9:  [],
   -9:  [],
   91:  [],
   10:  [],
   11:  [],
   12:  [],                      # mono → all independent
  -12:  [],
   13:  [],
  -13:  [],
   14:  [],                      # triclinic → all independent
}

# Which kwargs are accepted for each ibrav (maps to (nk1,nk2,nk3) argument names)
_IBRAV_KMESH_ARGS = {
    0:  ('nk1', 'nk2', 'nk3'),
    1:  ('nk',),                 # cubic: single value
    2:  ('nk',),
    3:  ('nk',),
   -3:  ('nk',),
    4:  ('nk1', 'nk3'),          # hex: in-plane + out-of-plane
    5:  ('nk',),                 # rhomb: single value
   -5:  ('nk',),
    6:  ('nk1', 'nk3'),          # st: in-plane + out-of-plane
    7:  ('nk1', 'nk3'),          # bct
    8:  ('nk1', 'nk2', 'nk3'),   # ortho: three independent
    9:  ('nk1', 'nk2', 'nk3'),
   -9:  ('nk1', 'nk2', 'nk3'),
   91:  ('nk1', 'nk2', 'nk3'),
   10:  ('nk1', 'nk2', 'nk3'),
   11:  ('nk1', 'nk2', 'nk3'),
   12:  ('nk1', 'nk2', 'nk3'),
  -12:  ('nk1', 'nk2', 'nk3'),
   13:  ('nk1', 'nk2', 'nk3'),
  -13:  ('nk1', 'nk2', 'nk3'),
   14:  ('nk1', 'nk2', 'nk3'),
}

_KMESH_ARG_DOC = {
    'nk' : 'Grid size along all three directions (cubic / rhombohedral symmetry).',
    'nk1': 'Grid size along b1 (and b2 for hex/tet where they must be equal).',
    'nk2': 'Grid size along b2.',
    'nk3': 'Grid size along b3 (out-of-plane for hex/tet).',
}


class KPointsAutoCard:
    """
    K_POINTS {automatic} card with ibrav-aware constructor.

    The accepted keyword arguments depend on ibrav:

    Cubic  (ibrav=1,2,3,-3)   : KPointsAutoCard(ibrav, nk=8)
    Hex/ST/BCT (ibrav=4,6,7)  : KPointsAutoCard(ibrav, nk1=8, nk3=6)
    Rhombohedal (ibrav=5,-5)  : KPointsAutoCard(ibrav, nk=8)
    Ortho/Mono/Tri (ibrav≥8)  : KPointsAutoCard(ibrav, nk1=6, nk2=8, nk3=4)
    Free cell (ibrav=0)        : KPointsAutoCard(ibrav, nk1=6, nk2=6, nk3=6)

    Shift flags sk1/sk2/sk3:
        1 = shift the grid by half a step (off-Γ)
        0 = Γ-centred (default)

    Examples
    --------
    >>> k = KPointsAutoCard(2, nk=8)          # FCC, 8×8×8
    >>> k = KPointsAutoCard(4, nk1=8, nk3=6) # hexagonal, 8×8×6
    >>> k = KPointsAutoCard(8, nk1=6, nk2=8, nk3=4)  # orthorhombic
    >>> k = KPointsAutoCard(2, nk=8, sk1=1, sk2=1, sk3=1)  # shifted FCC
    """

    def __init__(self, ibrav: int, sk1: int = 0, sk2: int = 0, sk3: int = 0, **kwargs):
        if ibrav not in _IBRAV_INFO:
            raise ValueError(f'ibrav={ibrav} not recognised.')

        self._ibrav = ibrav
        self.sk1 = sk1
        self.sk2 = sk2
        self.sk3 = sk3

        accepted = _IBRAV_KMESH_ARGS[ibrav]
        bad = [k for k in kwargs if k not in accepted and k not in ('nk1','nk2','nk3','nk')]
        if bad:
            raise TypeError(f"Unexpected keyword(s): {bad}")

        # Check no incompatible args were given
        for k in kwargs:
            if k not in accepted:
                raise TypeError(
                    f"'{k}' is not a valid mesh argument for ibrav={ibrav} "
                    f"({_IBRAV_INFO[ibrav]['label']}). "
                    f"Expected: {accepted}"
                )

        # Expand to nk1, nk2, nk3
        if 'nk' in accepted:
            nk = kwargs.get('nk')
            if nk is None:
                raise TypeError(f"'nk' is required for ibrav={ibrav} "
                                f"({_IBRAV_INFO[ibrav]['label']})")
            self.nk1 = self.nk2 = self.nk3 = int(nk)
        elif accepted == ('nk1', 'nk3'):
            nk1 = kwargs.get('nk1')
            nk3 = kwargs.get('nk3')
            if nk1 is None or nk3 is None:
                raise TypeError(f"'nk1' and 'nk3' are required for ibrav={ibrav} "
                                f"({_IBRAV_INFO[ibrav]['label']})")
            self.nk1 = self.nk2 = int(nk1)
            self.nk3 = int(nk3)
        else:
            for arg in ('nk1', 'nk2', 'nk3'):
                if arg not in kwargs:
                    raise TypeError(f"'{arg}' is required for ibrav={ibrav} "
                                    f"({_IBRAV_INFO[ibrav]['label']})")
            self.nk1 = int(kwargs['nk1'])
            self.nk2 = int(kwargs['nk2'])
            self.nk3 = int(kwargs['nk3'])

    @classmethod
    def info(cls, ibrav: int) -> str:
        """Describe the accepted mesh arguments for a given ibrav."""
        accepted = _IBRAV_KMESH_ARGS[ibrav]
        lines = [f"ibrav={ibrav}  {_IBRAV_INFO[ibrav]['label']}",
                 f"KPointsAutoCard accepts: {accepted}"]
        for a in accepted:
            lines.append(f"  {a:5s}  {_KMESH_ARG_DOC.get(a,'')}")
        return '\n'.join(lines)

    def to_string(self) -> str:
        return (
            'K_POINTS {automatic}\n'
            f'  {self.nk1} {self.nk2} {self.nk3}  {self.sk1} {self.sk2} {self.sk3}'
        )

    def __repr__(self) -> str:
        return (f'KPointsAutoCard(ibrav={self._ibrav}, '
                f'{self.nk1}×{self.nk2}×{self.nk3}, '
                f'shift={self.sk1}{self.sk2}{self.sk3})')


# ---------------------------------------------------------------------------
# CellParametersCard  (for ibrav=0)
# ---------------------------------------------------------------------------

class CellParametersCard:
    """
    CELL_PARAMETERS card.  Required when ibrav=0.

    Usage
    -----
    card = CellParametersCard('angstrom', [5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43])

    From ASE:
    card = CellParametersCard.from_atoms(atoms)
    """

    VALID_UNITS = ('bohr', 'angstrom', 'alat')

    def __init__(self, units: str, v1, v2, v3):
        units = units.lower()
        if units not in self.VALID_UNITS:
            raise ValueError(f"units must be one of {self.VALID_UNITS}")
        self.units = units
        self.v1 = [float(x) for x in v1]
        self.v2 = [float(x) for x in v2]
        self.v3 = [float(x) for x in v3]

    @classmethod
    def from_atoms(cls, atoms, units: str = 'angstrom') -> 'CellParametersCard':
        cell = atoms.cell
        if units == 'bohr':
            mat = cell[:] * _ANG_TO_BOHR
        else:
            mat = cell[:]  # angstrom
        return cls(units, mat[0], mat[1], mat[2])

    def to_string(self) -> str:
        def row(v):
            return '  ' + '  '.join(f'{x:16.10f}' for x in v)
        return (f'CELL_PARAMETERS {{{self.units}}}\n'
                + row(self.v1) + '\n'
                + row(self.v2) + '\n'
                + row(self.v3))

    def __repr__(self) -> str:
        return f'CellParametersCard(units={self.units!r})'


# ---------------------------------------------------------------------------
# PWInput  — top-level container
# ---------------------------------------------------------------------------

class PWInput:
    """
    Complete pw.x input file container.

    Attributes
    ----------
    control   : ControlNamelist
    system    : SystemNamelist
    electrons : ElectronsNamelist
    ions      : IonsNamelist       (optional)
    cell      : CellNamelist       (optional)
    atomic_species    : AtomicSpeciesCard
    atomic_positions  : AtomicPositionsCard
    k_points          : KPointsAutoCard  (or CellParametersCard-compatible object)
    cell_parameters   : CellParametersCard  (required if ibrav=0)

    Usage
    -----
    inp = PWInput(
        control   = ControlNamelist(calculation='scf', prefix='silicon'),
        system    = SystemNamelist.from_atoms(atoms, ibrav=2, ecutwfc=40),
        electrons = ElectronsNamelist(conv_thr=1e-8),
        atomic_species   = AtomicSpeciesCard.from_atoms(atoms, pseudos={'Si': 'Si.upf'}),
        atomic_positions = AtomicPositionsCard.from_atoms(atoms, units='crystal'),
        k_points  = KPointsAutoCard(2, nk=8),
    )
    print(inp.to_string())
    inp.write('scf.in')
    """

    def __init__(self,
                 control:   ControlNamelist,
                 system:    SystemNamelist,
                 electrons: ElectronsNamelist,
                 atomic_species:   AtomicSpeciesCard,
                 atomic_positions: AtomicPositionsCard,
                 k_points,
                 ions:      IonsNamelist    = None,
                 cell:      CellNamelist    = None,
                 cell_parameters: CellParametersCard = None):

        self.control   = control
        self.system    = system
        self.electrons = electrons
        self.ions      = ions
        self.cell      = cell
        self.atomic_species    = atomic_species
        self.atomic_positions  = atomic_positions
        self.k_points          = k_points
        self.cell_parameters   = cell_parameters

        # Sanity checks
        calc = control.get('calculation')
        if calc in ('relax', 'md', 'vc-relax', 'vc-md') and ions is None:
            import warnings
            warnings.warn(f"calculation='{calc}' but no IonsNamelist provided.")
        if calc in ('vc-relax', 'vc-md') and cell is None:
            import warnings
            warnings.warn(f"calculation='{calc}' but no CellNamelist provided.")
        if system.ibrav == 0 and cell_parameters is None:
            raise ValueError("ibrav=0 requires a CellParametersCard.")

    def to_string(self) -> str:
        blocks = [
            self.control.to_string(),
            self.system.to_string(),
            self.electrons.to_string(),
        ]
        if self.ions is not None:
            blocks.append(self.ions.to_string())
        if self.cell is not None:
            blocks.append(self.cell.to_string())
        blocks += [
            self.atomic_species.to_string(),
            self.atomic_positions.to_string(),
        ]
        if self.cell_parameters is not None:
            blocks.append(self.cell_parameters.to_string())
        blocks.append(self.k_points.to_string())
        return '\n\n'.join(blocks) + '\n'

    def write(self, filename: str) -> None:
        """Write the input to a file."""
        with open(filename, 'w') as f:
            f.write(self.to_string())
        print(f'Written: {filename}')

    def __repr__(self) -> str:
        calc = self.control.get('calculation')
        prefix = self.control.get('prefix')
        ibrav = self.system.ibrav
        return f'PWInput(calculation={calc!r}, prefix={prefix!r}, ibrav={ibrav})'


# ---------------------------------------------------------------------------
# Convenience factory
# ---------------------------------------------------------------------------

def pw_input_from_atoms(atoms,
                        ibrav: int,
                        ecutwfc: float,
                        pseudos: dict,
                        calculation: str = 'scf',
                        prefix: str = 'pwscf',
                        pseudo_dir: str = '/content/pseudo',
                        outdir: str = '/content/out',
                        k_points: KPointsAutoCard = None,
                        pos_units: str = 'crystal',
                        ecutrho: float = None,
                        **system_kwargs) -> PWInput:
    """
    One-shot factory: build a complete PWInput from an ASE Atoms object.

    Parameters
    ----------
    atoms       : ase.Atoms
    ibrav       : int
    ecutwfc     : float (Ry)
    pseudos     : dict  {symbol: pseudo_filename}  or  {symbol: (mass, filename)}
    calculation : str   default 'scf'
    prefix      : str
    pseudo_dir  : str
    outdir      : str
    k_points    : KPointsAutoCard  (required)
    pos_units   : str  atomic positions units, default 'crystal'
    ecutrho     : float (Ry), optional
    **system_kwargs : extra &SYSTEM parameters (e.g. occupations, degauss, …)

    Returns
    -------
    PWInput
    """
    if k_points is None:
        raise ValueError("k_points (a KPointsAutoCard) is required.")

    control = ControlNamelist(
        calculation=calculation,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=outdir,
        tprnfor=True,
        tstress=(calculation in ('vc-relax', 'vc-md')),
    )

    system = SystemNamelist.from_atoms(
        atoms, ibrav=ibrav, ecutwfc=ecutwfc, ecutrho=ecutrho, **system_kwargs
    )

    electrons = ElectronsNamelist()

    species_card   = AtomicSpeciesCard.from_atoms(atoms, pseudos)
    positions_card = AtomicPositionsCard.from_atoms(atoms, units=pos_units)

    cell_card = CellParametersCard.from_atoms(atoms) if ibrav == 0 else None

    return PWInput(
        control=control,
        system=system,
        electrons=electrons,
        atomic_species=species_card,
        atomic_positions=positions_card,
        k_points=k_points,
        cell_parameters=cell_card,
    )
