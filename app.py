import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol          # Keep py3Dmol import
import stmol            # Keep stmol import
import re
import numpy
from PIL import Image

# --- Helper Functions ---
def normalize_name(name):
    """Cleans molecule name"""
    name = name.strip()
    persian_nums = "۰۱۲۳۴۵۶۷۸۹"; arabic_nums = "٠١٢٣٤٥٦٧٨٩"; english_nums = "0123456789"
    name = name.translate(str.maketrans(persian_nums, english_nums)).translate(str.maketrans(arabic_nums, english_nums))
    name = re.sub(r'\s+', ' ', name).strip().lower()
    return name

def get_molecule_data(name):
    """Fetches data from PubChem"""
    try:
        results = pcp.get_compounds(name, 'name')
        return results[0] if results else None
    except pcp.PubChemHTTPError as e:
        st.error(f"PubChem Error: {'Not found' if 'NotFound' in str(e) else e}")
        return None
    except Exception as e:
        st.error(f"An unexpected error occurred during PubChem search: {e}")
        return None

# --- Streamlit User Interface ---
st.set_page_config(layout="wide")
st.title("🧪 Molecule Information Viewer ⌬")
st.markdown("Enter the English name of a molecule.")

# --- Get User Input ---
raw_molecule_name = st.text_input("Molecule Name (English):", placeholder="e.g., Water, Aspirin")

if raw_molecule_name:
    normalized_name = normalize_name(raw_molecule_name)
    st.write(f"Searching for: `{normalized_name}`")

    # --- Search and Display Information ---
    compound = None
    with st.spinner(f"Searching PubChem..."):
        compound = get_molecule_data(normalized_name)

    if compound:
        display_name = compound.iupac_name or raw_molecule_name
        st.success(f"Found: '{display_name}'!")

        # Basic Info
        st.subheader("Basic Information:")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Formula", compound.molecular_formula or "N/A")
        with col2:
            mol_weight_display = "N/A"
            if hasattr(compound, 'molecular_weight') and compound.molecular_weight is not None:
                try:
                    mol_weight_display = f"{float(compound.molecular_weight):.2f} g/mol"
                except (ValueError, TypeError):
                    mol_weight_display = f"{compound.molecular_weight} g/mol"
            st.metric("Mol. Weight", mol_weight_display)

        if compound.cid:
             st.markdown(f"**CID:** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")

        # Structures
        mol = None # RDKit mol object
        if hasattr(compound, 'isomeric_smiles') and compound.isomeric_smiles:
            try:
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
                if mol is None: st.warning("RDKit: MolFromSmiles failed.")
            except Exception as e:
                st.error(f"RDKit Error processing SMILES: {e}"); mol = None

        if mol:
            st.subheader("Molecule Structure:")
            col_2d, col_3d = st.columns(2)

            # 2D Structure
            with col_2d:
                st.markdown("**2D Structure:**")
                try:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                except Exception as e: st.error(f"Error generating 2D image: {e}")

            # --- 3D Structure using py3Dmol + stmol ---
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                if compound.cid:
                    try:
                        temp_sdf_file=f'c{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        with open(temp_sdf_file, 'r') as f: sdf_content = f.read()
                        if not sdf_content:
                            st.warning("Downloaded 3D SDF file was empty.")
                            sdf_content = None
                    except pcp.NotFoundError:
                        st.warning(f"3D SDF not found on PubChem (CID {compound.cid}).")
                    except Exception as e:
                        st.error(f"SDF Download Error: {e}")
                        sdf_content = None
                else:
                    st.warning("No PubChem CID available to download 3D structure.")

                if sdf_content:
                    try:
                        # 1. Create and configure the py3Dmol view object
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick':{}}) # Set style using py3Dmol
                        viewer.setBackgroundColor('0xeeeeee')
                        viewer.zoomTo()

                        # 2. Render the py3Dmol viewer object using stmol.showmol
                        stmol.showmol(viewer, height=400, width=400) # Pass the viewer object

                    except Exception as e:
                        st.error(f"Error rendering 3D view with py3Dmol/stmol: {e}")
                        st.error(f"Exception type: {type(e)}")
                # else: Handled by SDF download warnings
            # --- End of 3D Structure block ---

        # Synonyms
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names:")
             st.json(compound.synonyms[:5])

    elif raw_molecule_name:
        st.info("Search finished. See messages above for results or errors.")
else:
    st.info("Enter an English molecule name above.")
