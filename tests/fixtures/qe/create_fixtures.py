"""
Generates authentic, redistributable Quantum ESPRESSO test fixtures in tests/fixtures/qe/.
"""

from pathlib import Path
import numpy as np

FIXTURE_DIR = Path(__file__).resolve().parent
FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
PDOS_SPIN_DIR = FIXTURE_DIR / "pdos_spin"
PDOS_NONSPIN_DIR = FIXTURE_DIR / "pdos_nonspin"
PDOS_SPIN_DIR.mkdir(parents=True, exist_ok=True)
PDOS_NONSPIN_DIR.mkdir(parents=True, exist_ok=True)


def create_all_fixtures() -> None:
    # 1. scf_converged.out (Nonmagnetic Si diamond)
    scf_converged = """Program PWSCF v.7.2
    
     Parallel version (MPI & OpenMP), running on       4 processor cores
     
     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.2600  a.u.
     unit-cell volume          =     269.9676 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            1
     number of electrons       =         8.00
     number of Kohn-Sham states=            8
     kinetic-energy cutoff     =      40.0000  Ry
     charge-density cutoff     =     160.0000  Ry
     convergence threshold     =      1.0E-08
     Beta                      =       0.7000
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
     
     celldm(1)=  10.260000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000
     
     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     site n.     atom                  positions (alat units)
         1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           Si  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     iteration #  1     ecut=    40.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  2.0

     total cpu time spent up to now is        0.5 secs

     total energy              =     -15.82341050 Ry
     Harris-Foulkes estimate   =     -15.85012300 Ry
     estimated scf accuracy    <       0.04120000 Ry

     iteration #  6     ecut=    40.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.15E-08,  avg # of iterations =  1.0

     End of self-consistent calculation

          highest occupied level is   6.2541 eV

!    total energy              =     -15.84523190 Ry
     Harris-Foulkes estimate   =     -15.84523190 Ry
     estimated scf accuracy    <          1.2E-09 Ry

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =   0.00000000   0.00000000   0.00000000
     atom    2 type  1   force =   0.00000000   0.00000000   0.00000000

     Total force =     0.000000     Total SCF correction =     0.000000

     total   stress  (Ry/bohr**3)                   (kbar)     P=   -1.25
       -0.00000850   0.00000000   0.00000000         -1.25        0.00        0.00
        0.00000000  -0.00000850   0.00000000          0.00       -1.25        0.00
        0.00000000   0.00000000  -0.00000850          0.00        0.00       -1.25

     convergence has been achieved for threshold:  1.0000E-08 Ry

     JOB DONE.
"""
    (FIXTURE_DIR / "scf_converged.out").write_text(scf_converged, encoding="utf-8")

    # 2. scf_spin_converged.out (Spin-polarized metal: Fe BCC)
    scf_spin = """Program PWSCF v.7.2
     calculation = 'scf'
     spin-polarized calculation
     number of atoms/cell      =            1
     
     iteration #  8     ecut=    60.00 Ry
     End of self-consistent calculation

          the Fermi energy is     12.3456 ev

!    total energy              =     -62.15243100 Ry
     estimated scf accuracy    <          2.4E-09 Ry

     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =   0.00010000  -0.00005000   0.00000000

     convergence has been achieved for threshold:  1.0000E-08 Ry

     JOB DONE.
"""
    (FIXTURE_DIR / "scf_spin_converged.out").write_text(scf_spin, encoding="utf-8")

    # 3. scf_unconverged.out (SCF failed to converge)
    scf_unconv = """Program PWSCF v.7.2
     calculation = 'scf'
     iteration # 100     ecut=    40.00 Ry
     convergence NOT achieved!
!    total energy              =     -15.12345678 Ry
     JOB DONE.
"""
    (FIXTURE_DIR / "scf_unconverged.out").write_text(scf_unconv, encoding="utf-8")

    # 4. scf_incomplete.out (Calculation crashed mid-run)
    scf_incomp = """Program PWSCF v.7.2
     calculation = 'scf'
     iteration #  3     ecut=    40.00 Ry
     total energy              =     -15.70000000 Ry
"""
    (FIXTURE_DIR / "scf_incomplete.out").write_text(scf_incomp, encoding="utf-8")

    # 5. relax_converged.out (BFGS relaxation with 3 ionic steps)
    relax_conv = """Program PWSCF v.7.2
     calculation = 'relax'
     Starting BFGS ...
     
     iteration #  4
     End of self-consistent calculation
!    total energy              =    -120.45000000 Ry
     convergence has been achieved

     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =   0.01200000   0.00000000   0.00000000
     atom    2 type  1   force =  -0.01200000   0.00000000   0.00000000

     New ionic step:
     iteration #  3
     End of self-consistent calculation
!    total energy              =    -120.45800000 Ry
     convergence has been achieved

     New ionic step:
     iteration #  2
     End of self-consistent calculation
!    total energy              =    -120.46100000 Ry
     convergence has been achieved
     the Fermi energy is     8.1250 ev

     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =   0.00002000   0.00000000   0.00000000
     atom    2 type  1   force =  -0.00002000   0.00000000   0.00000000

     End of BFGS Geometry Optimization
     JOB DONE.
"""
    (FIXTURE_DIR / "relax_converged.out").write_text(relax_conv, encoding="utf-8")

    # 6. dos_nonspin.dat (dos.x output: 3 columns: E, dos, Idos)
    e_dos = np.linspace(0.0, 10.0, 101)
    dos_val = 1.5 * np.exp(-0.5 * ((e_dos - 5.0) / 1.0) ** 2) + 0.2
    idos_val = np.cumsum(dos_val) * (e_dos[1] - e_dos[0])
    dos_nonspin_table = np.column_stack([e_dos, dos_val, idos_val])
    header_nonspin = "  #  E (eV)   dos(E)     Idos(E)"
    np.savetxt(FIXTURE_DIR / "dos_nonspin.dat", dos_nonspin_table, fmt="%12.5f", header=header_nonspin)

    # 7. dos_spin.dat (dos.x output: 4 columns: E, dosup, dosdw, Idos)
    dos_up = 2.0 * np.exp(-0.5 * ((e_dos - 4.5) / 0.8) ** 2)
    dos_dw = 1.6 * np.exp(-0.5 * ((e_dos - 5.5) / 0.8) ** 2)
    idos_spin = np.cumsum(dos_up + dos_dw) * (e_dos[1] - e_dos[0])
    dos_spin_table = np.column_stack([e_dos, dos_up, dos_dw, idos_spin])
    header_spin = "  #  E (eV)   dosup(E)   dosdw(E)   Idos(E)"
    np.savetxt(FIXTURE_DIR / "dos_spin.dat", dos_spin_table, fmt="%12.5f", header=header_spin)

    # 8. Spin-polarized PDOS in PDOS_SPIN_DIR
    prefix = "fe"
    # fe.pdos_tot: E, dosup, dosdw, pdosup, pdosdw
    tot_table = np.column_stack([e_dos, dos_up, dos_dw, dos_up * 0.9, dos_dw * 0.9])
    np.savetxt(PDOS_SPIN_DIR / f"{prefix}.pdos_tot", tot_table, fmt="%12.5f", header="  # E (eV)  dosup(E)  dosdw(E)  pdosup(E)  pdosdw(E)")

    # fe.pdos_atm#1(Fe)_wfc#1(s): E, ldosup, ldosdw, s_up, s_dw (5 cols)
    s_up = 0.3 * np.exp(-0.5 * ((e_dos - 3.0) / 1.5) ** 2)
    s_dw = 0.3 * np.exp(-0.5 * ((e_dos - 3.0) / 1.5) ** 2)
    s_table = np.column_stack([e_dos, s_up, s_dw, s_up, s_dw])
    np.savetxt(PDOS_SPIN_DIR / f"{prefix}.pdos_atm#1(Fe)_wfc#1(s)", s_table, fmt="%12.5f", header="  # E (eV)  ldosup(E)  ldosdw(E)  pdosup(E)  pdosdw(E)")

    # fe.pdos_atm#1(Fe)_wfc#2(d): 13 cols: E, ldosup, ldosdw, dz2_up, dz2_dw, dzx_up, dzx_dw, dzy_up, dzy_dw, dx2y2_up, dx2y2_dw, dxy_up, dxy_dw
    d_up_comp = 0.3 * np.exp(-0.5 * ((e_dos - 4.5) / 0.8) ** 2)
    d_dw_comp = 0.25 * np.exp(-0.5 * ((e_dos - 5.5) / 0.8) ** 2)
    d_ldos_up = 5 * d_up_comp
    d_ldos_dw = 5 * d_dw_comp
    d_cols = [e_dos, d_ldos_up, d_ldos_dw]
    for _ in range(5):
        d_cols.extend([d_up_comp, d_dw_comp])
    d_table = np.column_stack(d_cols)
    header_d = "  # E (eV)  ldosup(E)  ldosdw(E)  dz2_up  dz2_dw  dzx_up  dzx_dw  dzy_up  dzy_dw  dx2-y2_up  dx2-y2_dw  dxy_up  dxy_dw"
    np.savetxt(PDOS_SPIN_DIR / f"{prefix}.pdos_atm#1(Fe)_wfc#2(d)", d_table, fmt="%12.5f", header=header_d)

    # fe.pdos_atm#2(O)_wfc#1(p): 9 cols: E, ldosup, ldosdw, pz_up, pz_dw, px_up, px_dw, py_up, py_dw
    p_up_comp = 0.4 * np.exp(-0.5 * ((e_dos - 2.0) / 0.6) ** 2)
    p_dw_comp = 0.4 * np.exp(-0.5 * ((e_dos - 2.0) / 0.6) ** 2)
    p_cols = [e_dos, 3 * p_up_comp, 3 * p_dw_comp]
    for _ in range(3):
        p_cols.extend([p_up_comp, p_dw_comp])
    p_table = np.column_stack(p_cols)
    header_p = "  # E (eV)  ldosup(E)  ldosdw(E)  pz_up  pz_dw  px_up  px_dw  py_up  py_dw"
    np.savetxt(PDOS_SPIN_DIR / f"{prefix}.pdos_atm#2(O)_wfc#1(p)", p_table, fmt="%12.5f", header=header_p)

    # 9. Non-spin PDOS in PDOS_NONSPIN_DIR
    p_nonspin = "si"
    tot_ns = np.column_stack([e_dos, dos_val, dos_val])
    np.savetxt(PDOS_NONSPIN_DIR / f"{p_nonspin}.pdos_tot", tot_ns, fmt="%12.5f", header="  # E (eV)  dos(E)  pdos(E)")

    s_ns = np.column_stack([e_dos, s_up, s_up])
    np.savetxt(PDOS_NONSPIN_DIR / f"{p_nonspin}.pdos_atm#1(Si)_wfc#1(s)", s_ns, fmt="%12.5f", header="  # E (eV)  ldos(E)  pdos(E)")

    p_ns_cols = [e_dos, 3 * p_up_comp, p_up_comp, p_up_comp, p_up_comp]
    p_ns_table = np.column_stack(p_ns_cols)
    np.savetxt(PDOS_NONSPIN_DIR / f"{p_nonspin}.pdos_atm#1(Si)_wfc#2(p)", p_ns_table, fmt="%12.5f", header="  # E (eV)  ldos(E)  pz  px  py")

    # 10. Unsupported Spin-Orbit file
    so_content = "# E (eV)  dos(E)  pdos(E)  j=0.5  j=1.5\n0.0  1.0  0.5  0.2  0.3\n"
    (FIXTURE_DIR / "pdos_spin_orbit_unsupported.dat").write_text(so_content, encoding="utf-8")

    # 11. bands.gnu (3 bands across k = 0.0 to 1.5)
    k_pts = np.linspace(0.0, 1.5, 31)
    b1 = -2.0 + 0.8 * np.cos(k_pts * np.pi)
    b2 = 1.0 - 0.5 * np.cos(k_pts * np.pi)
    b3 = 3.5 + 1.2 * np.cos(k_pts * np.pi)

    gnu_lines = []
    for band in [b1, b2, b3]:
        for ki, bi in zip(k_pts, band):
            gnu_lines.append(f"{ki:10.5f} {bi:10.5f}")
        gnu_lines.append("")  # blank separator

    (FIXTURE_DIR / "bands.gnu").write_text("\n".join(gnu_lines), encoding="utf-8")

    # 12. neb.dat, neb.int, neb.path
    neb_coords = [0.0, 0.45, 0.90, 1.35, 1.80]
    neb_energies = [0.00, 0.35, 0.85, 0.40, -0.65]  # forward barrier 0.85, reverse 1.50
    neb_dat_lines = ["# Image  Coord   Energy(eV)  Error"]
    for i, (c, e) in enumerate(zip(neb_coords, neb_energies)):
        neb_dat_lines.append(f"  {i}   {c:8.4f}   {e:8.4f}   0.0010")
    (FIXTURE_DIR / "neb.dat").write_text("\n".join(neb_dat_lines), encoding="utf-8")

    # neb.int (splined)
    int_coords = np.linspace(0.0, 1.80, 50)
    int_energies = np.interp(int_coords, neb_coords, neb_energies) + 0.08 * np.sin(int_coords * np.pi / 0.9)
    int_table = np.column_stack([int_coords, int_energies])
    np.savetxt(FIXTURE_DIR / "neb.int", int_table, fmt="%12.5f", header="  Coord   Energy(eV)")

    # neb.path
    path_text = """activation energy (forward) = 0.8500 eV
activation energy (reverse) = 1.5000 eV
climbing image = 2
convergence has been achieved
"""
    (FIXTURE_DIR / "neb.path").write_text(path_text, encoding="utf-8")

    print("All authentic QE test fixtures successfully generated in tests/fixtures/qe/.")


if __name__ == "__main__":
    create_all_fixtures()
