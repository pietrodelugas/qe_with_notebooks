def build_mgo_input(atoms, ecutwfc, nk, prefix):
    # &CONTROL: scf calculation, set prefix, pseudo_dir, outdir; enable tprnfor and tstress
    control   = ControlNamelist(...)

    # &SYSTEM: ibrav=2 (FCC rocksalt), set ecutwfc
    system    = SystemNamelist.from_atoms(...)

    # &ELECTRONS: use defaults
    electrons = ElectronsNamelist()

    # ATOMIC_SPECIES card
    species   = AtomicSpeciesCard.from_atoms(...)

    # ATOMIC_POSITIONS in crystal (fractional) coordinates
    positions = AtomicPositionsCard.from_atoms(...)

    # K_POINTS {automatic}: ibrav=2 → all lattice vectors equivalent
    kpoints   = KPointsAutoCard(...)

    return PWInput(
        control=control, system=system, electrons=electrons,
        atomic_species=species, atomic_positions=positions, k_points=kpoints,
    )
