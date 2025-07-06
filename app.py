import streamlit as st
import psi4
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

def run_psi4_calculation(geometry_string, orbitals_to_plot):
    """
    Udfører en Hartree-Fock-beregning og genererer cube-filer.
    """
    psi4.core.clean()  # Ryd op i gamle filer
    psi4.set_memory('1 GB')
    psi4.set_num_threads(2)

    molecule = psi4.geometry(geometry_string)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'pk',
    })

    # Kør SCF-energiberegningen
    scf_e, scf_wfn = psi4.energy('scf', return_wfn=True)

    # Generer cube-filer for de specificerede orbitaler
    psi4.set_options({
        'cubeprop_tasks': ['orbitals'],
        'cubeprop_orbitals': orbitals_to_plot
    })
    psi4.cubeprop(scf_wfn)
    
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
        "orbitals": [2, 3, 4, 5, 6]
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
        "orbitals": [2, 3, 4, 5, 6]
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
        "orbitals": [2, 3, 4, 5, 6]
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
        "orbitals": [7, 8, 9, 10]
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
        "orbitals": [7, 8, 9]
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
        "orbitals": [6, 7, 8, 9]
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

    st.header("2. Kør beregningen")
    run_button = st.button(f"Beregn orbitaler for {molecule_name.capitalize()}")

with col2:
    st.header("3. Visualiser resultater")
    
    if run_button:
        with st.spinner(f"Kører Hartree-Fock beregning for {molecule_name}... Dette kan tage et øjeblik."):
            try:
                energy = run_psi4_calculation(
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
                
                geo_lines = selected_molecule_data["geometry"].strip().split('\n')[2:]
                num_atoms = len(geo_lines)
                xyz_geometry = f"{num_atoms}\n{molecule_name}\n" + "\n".join(geo_lines)
                
                st.subheader(f"Orbital #{orbital_index}")
                viewer = create_orbital_viewer(xyz_geometry, cube_file)
                st.components.v1.html(viewer._make_html(), height=420, width=620)
                
            else:
                st.warning(f"Kunne ikke finde .cube fil for orbital {orbital_index} (Søgte efter: {search_pattern})")
    else:
        st.info("Klik på knappen til venstre for at starte en beregning og se resultaterne her.")
