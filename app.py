import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
import streamlit.components.v1 as components
import re
import numpy # Ensure numpy is imported if needed by dependencies
from PIL import Image
import uuid # To generate unique IDs for divs

# --- Helper Functions (Keep as before) ---
def normalize_name(name):
    name = name.strip()
    persian_nums = "۰۱۲۳۴۵۶۷۸۹"; arabic_nums = "٠١٢٣٤٥٦٧٨٩"; english_nums = "0123456789"
    name = name.translate(str.maketrans(persian_nums, english_nums)).translate(str.maketrans(arabic_nums, english_nums))
    name = re.sub(r'\s+', ' ', name).strip().lower()
    return name

def get_molecule_data(name):
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
        mol = None
        if hasattr(compound, 'isomeric_smiles') and compound.isomeric_smiles:
            try:
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
                if mol is None: st.warning("RDKit: MolFromSmiles failed.")
            except Exception as e:
                st.error(f"RDKit Error: {e}"); mol = None

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

            # 3D Structure
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                if compound.cid:
                    try:
                        st.write(f"[Debug] Downloading SDF CID: {compound.cid}")
                        temp_sdf_file=f'c{compound.cid}_3d.sdf'; pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        with open(temp_sdf_file, 'r') as f: sdf_content = f.read()
                        if sdf_content: st.write(f"[Debug] SDF OK (len:{len(sdf_content)})")
                        else: st.warning("[Debug] SDF empty"); sdf_content = None
                    except pcp.NotFoundError: st.warning(f"3D SDF not found (CID {compound.cid}).")
                    except Exception as e: st.error(f"SDF Download Error: {e}")
                else: st.warning("No CID for 3D SDF download.")

                if sdf_content:
                    try:
                        st.write("[Debug] Configuring py3Dmol...")
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}}); viewer.setBackgroundColor('0xeeeeee'); viewer.zoomTo()
                        st.write("[Debug] Configured. Generating manual HTML...")

                        # --- MANUAL HTML WORKAROUND V2 (Call viewer.js()) ---
                        viewer_js_content = None
                        try:
                            st.write(f"[Debug] Checking viewer.js type: {type(viewer.js)}")
                            if callable(viewer.js):
                                st.write("[Debug] Calling viewer.js()...")
                                viewer_js_content = viewer.js() # *** CALL THE METHOD ***
                            else:
                                st.write("[Debug] Accessing viewer.js attribute...")
                                viewer_js_content = viewer.js
                            
                            if not viewer_js_content or not isinstance(viewer_js_content, str):
                                raise TypeError(f"viewer.js did not return a valid string. Got: {type(viewer_js_content)}")
                            st.write(f"[Debug] viewer.js content OK (len:{len(viewer_js_content)})")

                        except Exception as js_err:
                            st.error(f"[Debug] Error getting viewer.js content: {js_err}")
                            raise # Stop if we can't get the JS

                        viewer_div_id = f"molviewer-{uuid.uuid4()}"
                        py3dmol_js_url = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"
                        jquery_url = "https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"

                        manual_html = f"""
                            <script src="{jquery_url}"></script>
                            <script src="{py3dmol_js_url}"></script>
                            <div id="{viewer_div_id}" style="height: 400px; width: 400px; position: relative;"></div>
                            <script>
                                (function() {{
                                    try {{
                                        var element = $('#{viewer_div_id}'); var config = {{}};
                                        var viewer = $3Dmol.createViewer(element, config);
                                        {viewer_js_content} // Inject JS commands
                                        viewer.render();
                                    }} catch (e) {{ console.error("3Dmol error:", e); $('#{viewer_div_id}').text("Error loading 3D view."); }}
                                }})();
                            </script>
                        """
                        st.write(f"[Debug] Manual HTML OK (len:{len(manual_html)}). Embedding...")
                        components.html(manual_html, height=410, width=410)
                        st.write("[Debug] Embedding attempted.")
                        # --- END WORKAROUND ---

                    except Exception as e:
                        st.error(f"Error in 3D View Setup/Generation: {e}")
                        st.error(f"Exception type: {type(e)}")

        # Synonyms
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names:")
             st.json(compound.synonyms[:5])

    elif raw_molecule_name:
        st.info("Search finished. See messages above.")
else:
    st.info("Enter an English molecule name.")
