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
    persian_nums = "۰۱۲۳۴۵۶۷۸۹"
    arabic_nums = "٠١٢٣٤٥٦٧٨٩"
    english_nums = "0123456789"
    translation_table_persian = str.maketrans(persian_nums, english_nums)
    translation_table_arabic = str.maketrans(arabic_nums, english_nums)
    name = name.translate(translation_table_persian)
    name = name.translate(translation_table_arabic)
    name = re.sub(r'\s+', ' ', name).strip()
    name = name.lower()
    return name

def get_molecule_data(name):
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0]
        else:
            return None
    except pcp.PubChemHTTPError as e:
        if "PUGREST.NotFound" in str(e):
             st.error(f"Molecule '{name}' not found in PubChem.")
        else:
             st.error(f"PubChem Error: {e}")
        return None
    except Exception as e:
        st.error(f"An unexpected error occurred during PubChem search: {e}")
        return None

# --- Streamlit User Interface ---
st.set_page_config(layout="wide")
st.title("🧪 Molecule Information Viewer ⌬")
st.markdown("""
Enter the English name of a molecule. The app will try to display its information and structure.
""")

# --- Get User Input ---
raw_molecule_name = st.text_input("Molecule Name (English):", placeholder="e.g., Water, Aspirin, 1-chlorobutane")

if raw_molecule_name:
    normalized_name = normalize_name(raw_molecule_name)
    st.write(f"Searching for: `{normalized_name}`")

    # --- Search and Display Information ---
    compound = None
    with st.spinner(f"Searching PubChem for '{normalized_name}'..."):
        compound = get_molecule_data(normalized_name)

    if compound:
        display_name = compound.iupac_name if compound.iupac_name else raw_molecule_name
        st.success(f"Found information for '{display_name}'!")

        # Basic Info
        st.subheader("Basic Information:")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Molecular Formula", compound.molecular_formula or "N/A")
        with col2:
            mol_weight_display = "N/A"
            if hasattr(compound, 'molecular_weight') and compound.molecular_weight is not None:
                try:
                    weight_float = float(compound.molecular_weight)
                    mol_weight_display = f"{weight_float:.2f} g/mol"
                except (ValueError, TypeError):
                    mol_weight_display = str(compound.molecular_weight)
                    if "g/mol" not in mol_weight_display.lower():
                         mol_weight_display += " g/mol"
            st.metric("Molecular Weight", mol_weight_display)

        if compound.cid:
             st.markdown(f"**PubChem CID:** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")

        # Structures
        mol = None
        if hasattr(compound, 'isomeric_smiles') and compound.isomeric_smiles:
            try:
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
                if mol is None:
                    st.warning("Could not generate molecule object from SMILES string.")
            except Exception as e:
                st.error(f"Error processing SMILES with RDKit: {e}")
                mol = None

        if mol:
            st.subheader("Molecule Structure:")
            col_2d, col_3d = st.columns(2)

            # 2D Structure
            with col_2d:
                st.markdown("**2D Structure:**")
                try:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                except Exception as e:
                    st.error(f"Error generating 2D image: {e}")

            # 3D Structure
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                if compound.cid:
                    try:
                        st.write(f"[Debug] Attempting to download 3D SDF for CID: {compound.cid}")
                        temp_sdf_file = f'cid_{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        with open(temp_sdf_file, 'r') as f:
                            sdf_content = f.read()
                        if sdf_content:
                             st.write(f"[Debug] SDF content downloaded (length: {len(sdf_content)})")
                        else:
                             st.warning("[Debug] Downloaded 3D SDF file was empty.")
                             sdf_content = None
                    except pcp.NotFoundError:
                        st.warning(f"3D structure (SDF) not found on PubChem for CID {compound.cid}.")
                    except Exception as e:
                        st.error(f"Error downloading 3D SDF from PubChem: {e}")
                else:
                    st.warning("Cannot fetch 3D structure without a PubChem CID.")

                if sdf_content:
                    try:
                        st.write("[Debug] Configuring py3Dmol viewer...")
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}})
                        viewer.setBackgroundColor('0xeeeeee')
                        viewer.zoomTo()
                        st.write("[Debug] Viewer configured. Attempting MANUAL HTML generation...")

                        # --- MANUAL HTML GENERATION WORKAROUND ---
                        # 1. Get the necessary JS commands from the viewer object
                        viewer_js = viewer.js
                        if not viewer_js:
                             raise ValueError("viewer.js attribute is empty, cannot generate manual HTML.")
                        st.write(f"[Debug] viewer.js content obtained (length: {len(viewer_js)})")

                        # 2. Define a unique ID for the container div
                        viewer_div_id = f"molviewer-{uuid.uuid4()}"

                        # 3. Construct the HTML string
                        # Using CDN links for py3Dmol and jQuery (py3Dmol dependency)
                        py3dmol_js_url = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"
                        jquery_url = "https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js" # Using a recent jQuery

                        manual_html = f"""
                            <script src="{jquery_url}"></script>
                            <script src="{py3dmol_js_url}"></script>
                            <div id="{viewer_div_id}" style="height: 400px; width: 400px; position: relative;"></div>
                            <script>
                                (function() {{
                                    try {{
                                        let element = $('#{viewer_div_id}');
                                        let config = {{}}; // Default config
                                        console.log("Creating 3Dmol viewer in div:", "{viewer_div_id}");
                                        let viewer = $3Dmol.createViewer(element, config);
                                        console.log("Applying viewer commands...");
                                        {viewer_js} // Inject the specific commands for this viewer instance
                                        console.log("Rendering viewer...");
                                        viewer.render();
                                        console.log("Viewer rendering process initiated.");
                                    }} catch (error) {{
                                        console.error("Error during 3Dmol manual script execution:", error);
                                        // Optionally display error in the div
                                        $('#{viewer_div_id}').html('<p style="color:red;">Error initializing 3D viewer. Check browser console.</p>');
                                    }}
                                }})();
                            </script>
                        """
                        st.write(f"[Debug] Manual HTML string created (length: {len(manual_html)})")

                        # 4. Embed using Streamlit components
                        st.write("[Debug] Embedding manual HTML component...")
                        components.html(manual_html, height=410, width=410)
                        st.write("[Debug] Manual HTML component embedding attempted.")
                        # --- END OF MANUAL HTML GENERATION WORKAROUND ---

                    except Exception as e:
                        # Catch errors during viewer setup or manual HTML generation itself
                        st.error(f"Error during py3Dmol setup or manual HTML generation: {e}")
                        st.error(f"Exception type: {type(e)}")

        # Synonyms
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names (Synonyms):")
             st.json(compound.synonyms[:5])

    elif raw_molecule_name:
        st.info("Search finished. See messages above for results or errors.")

else:
    st.info("Please enter an English molecule name in the box above.")
