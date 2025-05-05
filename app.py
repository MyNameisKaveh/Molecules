import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
# Removed py3Dmol, components, uuid
import stmol          # Added stmol import
import re
import numpy          # Keep numpy import
from PIL import Image

# --- Helper Functions (Keep as before) ---
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
                    mol_weight_display = f"{compound.molecular_weight} g/mol" # Display as is if conversion fails
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

            # --- 3D Structure using stmol ---
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                # Download SDF content (same logic as before)
                if compound.cid:
                    try:
                        st.write(f"[Debug] Downloading SDF CID: {compound.cid}") # Keep debug msg
                        temp_sdf_file=f'c{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        with open(temp_sdf_file, 'r') as f: sdf_content = f.read()
                        if sdf_content:
                            st.write(f"[Debug] SDF OK (len:{len(sdf_content)})") # Keep debug msg
                        else:
                            st.warning("[Debug] SDF empty") # Keep debug msg
                            sdf_content = None
                    except pcp.NotFoundError:
                        st.warning(f"3D SDF not found on PubChem (CID {compound.cid}).")
                    except Exception as e:
                        st.error(f"SDF Download Error: {e}")
                else:
                    st.warning("No PubChem CID available to download 3D structure.")

                # Render using stmol if SDF content exists
                if sdf_content:
                    try:
                        st.write("[Debug] Rendering with stmol...") # Keep debug msg
                        # --- Use stmol.showmol ---
                        stmol.showmol(sdf_content, style="stick", height=400, width=400)
                        # --- End of stmol usage ---
                        st.write("[Debug] stmol rendering attempted.") # Keep debug msg
                    except Exception as e:
                        st.error(f"Error rendering 3D view with stmol: {e}")
                        st.error(f"Exception type: {type(e)}")
                # else: No need for else, handled by SDF download warnings

            # --- End of 3D Structure block ---

        # Synonyms
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names:")
             # Display first 5 synonyms
             st.json(compound.synonyms[:5])

    elif raw_molecule_name: # Only show if search was attempted but failed
        st.info("Search finished. See messages above for results or errors.")

else: # Initial state
    st.info("Enter an English molecule name above.")
