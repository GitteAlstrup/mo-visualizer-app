import streamlit as st
from pyscf import gto, scf, tools
import numpy as np
import py3Dmol
import glob
import os

def create_orbital_viewer(xyz_geometry, cube_file_path):
    """
    Genererer en interaktiv py3Dmol-viewer for et givent molekyle og en cube-fil.
    """
    with open(cube_file_path, 'r') as f:
        cube_data = f.read()

    view = py3Dmol.view(width=600, height=400)
    
    # Tilføj molekylemodellen og sæt stilen
    view.addModel(xyz_geometry, 'xyz')
    # Vis molekylet med atomer som kugler og bindinger som en wireframe for bedre klarhed
    view.setStyle({'sphere': {'radius': 0.3}, 'line': {}})
    
    # Tilføj de positive og negative lapper af orbitalen
    view.addVolumetricData(cube_data, 'cube', {'isoval': 0.03, 'color': "blue", 'opacity': 0.75})
    view.addVolumetricData(cube_data, 'cube', {'isoval': -0.03, 'color': "red", 'opacity': 0.75})
    
    view.zoomTo()
    return view

def run_pyscf_calculation(geometry_string, orbitals_to_plot):
    """
    Udfører en Hartree-Fock-beregning med PySCF og genererer cube-filer.
    """
    # Ryd op i gamle .cube filer
    for f in glob.glob("*.cube"):
        os.remove(f)

    # Konverter Psi4-geometri til PySCF-format
    lines = geometry_string.strip().split('\n')[1:]
    pyscf_geom = [f"{line.split()[0]} {line.split()[1]} {line.split()[2]} {line.split()[3]}" for line in lines]

    # Byg molekylet i PySCF
    mol = gto.Mole()
    mol.atom = "\n".join(pyscf_geom)
    mol.basis = 'sto-3g'
    mol.build()

    # Kør SCF-beregningen
    mf = scf.RHF(mol).run()
    scf_e = mf.e_tot

    # Generer cube-filer for de ønskede orbitaler
    for i in orbitals_to_plot:
        # PySCF's orbital-indeksering starter ved 0, så vi trækker 1 fra.
        orbital_index_pyscf = i - 1
        # Bevar det gamle filnavnsmønster for at sikre kompatibilitet
        # (Symmetri-delen er nu en pladsholder, da PySCF ikke giver den let)
        filename = f"Psi_a_{i}_A1.cube"
        tools.cubegen.orbital(mol, filename, mf.mo_coeff[:, orbital_index_pyscf])
    
    return scf_e

# --- Data for alle molekyler ---
molecules = {
    "water": {
        "geometry": """
            0 1
            O   0.000000    0.000000    0.117300
            H   0.000000    0.757200   -0.469200
            H   0.000000   -0.757200   -0.469200
        """,
        # 10 elektroner = 5 besatte MOs. HOMO=5 (lone pair), LUMO=6.
        "orbitals": [2, 3, 4, 5, 6],
        "orbital_labels": {
            4: "(σ binding)",
            5: "(HOMO, n-orbital)",
            6: "(LUMO, σ* antibindende)"
        },
        "description": """
        Vand er et godt eksempel på, hvordan MO-teori beskriver både **sigma (σ) bindinger** og **ikke-bindende orbitaler (lone pairs)**.
        *   **Orbital #4** viser en af de bindende σ-orbitaler.
        *   **Orbital #5 (HOMO)** er primært lokaliseret på iltatomet og repræsenterer et af de to 'lone pairs'.
        """
    },
    "ammonia": {
        "geometry": """
            0 1
            N   0.000000    0.000000    0.115329
            H   0.000000    0.939563   -0.269099
            H   0.813654   -0.469781   -0.269099
            H  -0.813654   -0.469781   -0.269099
        """,
        # 10 elektroner = 5 besatte MOs. HOMO=5 (lone pair), LUMO=6.
        "orbitals": [2, 3, 4, 5, 6],
        "orbital_labels": {
            5: "(HOMO, n-orbital)",
            6: "(LUMO)"
        }
    },
    "methane": {
        "geometry": """
            0 1
            C   0.000000   0.000000   0.000000
            H   0.629118   0.629118   0.629118
            H  -0.629118  -0.629118   0.629118
            H  -0.629118   0.629118  -0.629118
            H   0.629118  -0.629118  -0.629118
        """,
        # 10 elektroner = 5 besatte MOs. HOMO=5, LUMO=6.
        "orbitals": [2, 3, 4, 5, 6],
        "orbital_labels": {
            5: "(HOMO)",
            6: "(LUMO)"
        },
        "description": """
        Methan er det klassiske eksempel på **sp³-hybridisering**. I molekylorbitalteori er de fire C-H bindinger ikke repræsenteret af fire identiske, lokaliserede orbitaler.
        
        I stedet er de beskrevet af en kombination af orbitalerne #2, #3, #4 og #5, som er degenererede (har samme energi) og tilsammen skaber den velkendte tetraediske elektronfordeling.
        """
    },
    "ethane": {
        "geometry": """
            0 1
            C   0.000000    0.000000    0.765000
            C   0.000000    0.000000   -0.765000
            H   1.026621    0.000000    1.155000
            H  -0.513311   -0.889165    1.155000
            H  -0.513311    0.889165    1.155000
            H  -1.026621    0.000000   -1.155000
            H   0.513311    0.889165   -1.155000
            H   0.513311   -0.889165   -1.155000
        """,
        # 18 elektroner = 9 besatte MOs. HOMO=9, LUMO=10.
        "orbitals": [7, 8, 9, 10],
        "orbital_labels": {
            7: "(C-C σ-binding)",
            9: "(HOMO)",
            10: "(LUMO)"
        }
    },
    "ethene": {
        "geometry": """
            0 1
            C   0.000000    0.000000    0.670000
            C   0.000000    0.000000   -0.670000
            H   0.000000    0.920000    1.230000
            H   0.000000   -0.920000    1.230000
            H   0.000000   -0.920000   -1.230000
            H   0.000000    0.920000   -1.230000
        """,
        # 16 elektroner = 8 besatte MOs. HOMO=8 (pi), LUMO=9 (pi*).
        "orbitals": [7, 8, 9],
        "orbital_labels": {
            7: "(C-C σ-binding)",
            8: "(HOMO, π-binding)",
            9: "(LUMO, π*-antibindende)"
        },
        "description": """
        Ethen er det perfekte eksempel på **sp²-hybridisering**. Her kan man tydeligt adskille σ- og π-systemerne.
        *   **Orbital #7** viser den stærke **sigma (σ) binding** mellem de to carbonatomer.
        *   **Orbital #8 (HOMO)** er den berømte **pi (π) binding**, som dannes af p-orbitalerne og er afgørende for molekylets reaktivitet.
        *   **Orbital #9 (LUMO)** er den antibindende **pi-stjerne (π*) orbital**, som er tom og vigtig for at forstå molekylets reaktioner og UV/Vis-spektroskopi."""
    },
    "ethyne": {
        "geometry": """
            0 1
            C   0.000000    0.000000    0.602500
            C   0.000000    0.000000   -0.602500
            H   0.000000    0.000000    1.666500
            H   0.000000    0.000000   -1.666500
        """,
        # 14 elektroner = 7 besatte MOs. HOMO=7.
        "orbitals": [6, 7, 8, 9],
        "orbital_labels": {
            6: "(π-binding #1)",
            7: "(HOMO, π-binding #2)",
            8: "(LUMO, π*-antibindende)"
        }
    },
}

# --- Streamlit Web Application Interface ---

st.set_page_config(page_title="MO Visualisering", layout="wide")
st.title("Interaktiv Molekylær Orbital Visualisering")
st.write("Dette værktøj lader dig køre en kvantekemisk beregning og visualisere de resulterende molekylære orbitaler direkte i browseren.")

col1, col2 = st.columns([1, 2])

with col1:
    st.header("1. Vælg et molekyle")
    molecule_name = st.selectbox(
        "Vælg fra listen:",
        list(molecules.keys()),
        index=0 # Start med det første molekyle i listen
    )
    
    selected_molecule_data = molecules[molecule_name]
    
    st.info(f"**Geometri for {molecule_name.capitalize()}:**")
    st.code(selected_molecule_data["geometry"].strip(), language='text')

    if "description" in selected_molecule_data:
        st.markdown("---")
        st.markdown(selected_molecule_data["description"])

    st.header("2. Kør beregningen")
    run_button = st.button(f"Beregn orbitaler for {molecule_name.capitalize()}")

with col2:
    st.header("3. Visualiser resultater")
    
    if run_button:
        with st.spinner(f"Kører Hartree-Fock beregning for {molecule_name}... Dette kan tage et øjeblik."):
            try:
                energy = run_pyscf_calculation(
                    selected_molecule_data["geometry"],
                    selected_molecule_data["orbitals"]
                )
                st.success(f"Beregning færdig! Hartree-Fock Energi: {energy:.6f} Hartrees")
            except Exception as e:
                st.error(f"En fejl opstod under beregningen: {e}")
                st.stop()

        for orbital_index in selected_molecule_data["orbitals"]:
            search_pattern = f"Psi_a_{orbital_index}_*.cube"
            found_files = glob.glob(search_pattern)
            
            if found_files:
                cube_file = found_files[0]
                
                geo_lines = selected_molecule_data["geometry"].strip().split('\n')[1:]
                num_atoms = len(geo_lines)
                xyz_geometry = f"{num_atoms}\n{molecule_name}\n" + "\n".join(geo_lines)
                
                # Hent label for orbitalen, hvis den findes
                orbital_labels = selected_molecule_data.get("orbital_labels", {})
                label = orbital_labels.get(orbital_index, "")
                
                st.subheader(f"Orbital #{orbital_index} {label}")
                viewer = create_orbital_viewer(xyz_geometry, cube_file)
                st.components.v1.html(viewer._make_html(), height=420, width=620)
                
            else:
                st.warning(f"Kunne ikke finde .cube fil for orbital {orbital_index} (Søgte efter: {search_pattern})")
    else:
        st.info("Klik på knappen til venstre for at starte en beregning og se resultaterne her.")
