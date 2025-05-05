import streamlit as st
from streamlit import cache_data # Import cache_data decorator
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol          # Keep py3Dmol import
import stmol            # Keep stmol import
import translators as ts # Import translators library
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

# --- Cache PubChem search results ---
@st.cache_data(ttl=3600) # Cache results for 1 hour (3600 seconds)
def get_molecule_data(name):
    """Fetches data from PubChem"""
    st.write(f"[Cache] Searching PubChem for: {name}") # To see when cache is missed
    try:
        results = pcp.get_compounds(name, 'name')
        # Return only the first result, or None if no results
        return results[0] if results else None
    except pcp.PubChemHTTPError as e:
        # Don't cache errors generally, but handle not found specifically
        if "NotFound" in str(e):
            st.error(f"Molecule '{name}' not found in PubChem.")
            return None # Return None for not found
        else:
            st.error(f"PubChem Error: {e}")
            # Re-raise other HTTP errors so they aren't cached as None
            raise e
    except Exception as e:
        st.error(f"An unexpected error occurred during PubChem search: {e}")
        # Re-raise other exceptions
        raise e

# --- Function to safely get attribute ---
def safe_get_attr(obj, attr, default="N/A"):
    """Safely get an attribute from an object, return default if not found"""
    if obj and hasattr(obj, attr):
        value = getattr(obj, attr)
        return value if value is not None else default
    return default

# --- Streamlit User Interface ---
st.set_page_config(layout="wide")
st.title("🧪 Advanced Molecule Information Viewer ⌬")
st.markdown("Search for molecules by English name, view properties, structures, and use the translator.")

# --- Translator Section ---
st.divider()
with st.expander("🌍 Persian to English Name Translator"):
    persian_name = st.text_area("Enter Persian name here:")
    if st.button("Translate to English"):
        if persian_name:
            try:
                with st.spinner("Translating..."):
                    # Use translators library (e.g., Google engine)
                    english_translation = ts.translate_text(persian_name, translator='google', from_language='fa', to_language='en')
                st.success(f"Suggested English name: `{english_translation}`")
                st.caption("Copy this name and paste it into the main search box below.")
            except Exception as e:
                st.error(f"Translation failed: {e}")
                st.info("Translation services can sometimes be unreliable. Please try again later or use a different tool.")
        else:
            st.warning("Please enter a Persian name first.")
st.divider()


# --- Main Search Section ---
raw_molecule_name = st.text_input("Molecule Name (English):", placeholder="e.g., Water, Aspirin, caffeine")

if raw_molecule_name:
    normalized_name = normalize_name(raw_molecule_name)
    st.write(f"Normalized search term: `{normalized_name}`")

    # --- Search and Display Information (using cached function) ---
    compound = None
    error_occurred = False
    try:
        with st.spinner(f"Searching... (Using cache if available)"):
            compound = get_molecule_data(normalized_name)
    except Exception as e:
        # Catch errors raised by the cached function if they weren't handled inside
        # Error message is already displayed inside get_molecule_data
        error_occurred = True


    if compound:
        display_name = safe_get_attr(compound, 'iupac_name', raw_molecule_name)
        st.success(f"Found: '{display_name}'!")

        # --- Organize results into Tabs ---
        tab_basic, tab_structure, tab_ids = st.tabs(["📊 Basic Info", "🔬 Structures", "🔑 Identifiers & Names"])

        # --- Basic Info Tab ---
        with tab_basic:
            st.subheader("General Properties")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Formula", safe_get_attr(compound, 'molecular_formula'))

                # Safely get melting point
                mp = safe_get_attr(compound, 'melting_point')
                st.metric("Melting Point", f"{mp} °C" if mp != "N/A" else "N/A")

            with col2:
                # Safely get and format molecular weight
                mw_val = safe_get_attr(compound, 'molecular_weight')
                mol_weight_display = "N/A"
                if mw_val != "N/A":
                    try:
                        mol_weight_display = f"{float(mw_val):.2f} g/mol"
                    except (ValueError, TypeError):
                        mol_weight_display = f"{mw_val} g/mol"
                st.metric("Mol. Weight", mol_weight_display)

                # Safely get boiling point
                bp = safe_get_attr(compound, 'boiling_point')
                st.metric("Boiling Point", f"{bp} °C" if bp != "N/A" else "N/A")

            # Add density if available (example, more properties can be added)
            density = safe_get_attr(compound, 'density')
            if density != "N/A":
                st.write(f"**Density:** {density}")


        # --- Structures Tab ---
        with tab_structure:
            st.subheader("Molecular Structures")
            mol_rdkit = None # RDKit mol object for structure generation
            smiles = safe_get_attr(compound, 'isomeric_smiles')
            if smiles != "N/A":
                try:
                    mol_rdkit = Chem.MolFromSmiles(smiles)
                    if mol_rdkit is None: st.warning("RDKit: Could not create molecule from SMILES.")
                except Exception as e:
                    st.error(f"RDKit Error processing SMILES: {e}"); mol_rdkit = None

            if mol_rdkit:
                col_2d, col_3d_config = st.columns([2, 3]) # Adjust column widths

                # 2D Structure
                with col_2d:
                    st.markdown("**2D Structure:**")
                    try:
                        img = Draw.MolToImage(mol_rdkit, size=(300, 300))
                        st.image(img)
                    except Exception as e: st.error(f"Error generating 2D image: {e}")

                # 3D Structure Configuration and Display
                with col_3d_config:
                    st.markdown("**3D Structure (Interactive):**")
                    # Style selection
                    style_options = ['stick', 'sphere', 'line', 'cartoon']
                    selected_style = st.selectbox("Select 3D Style:", style_options)

                    sdf_content = None
                    cid = safe_get_attr(compound, 'cid')
                    if cid != "N/A":
                        try:
                            temp_sdf_file=f'c{cid}_3d.sdf'
                            pcp.download('SDF', temp_sdf_file, cid, 'cid', record_type='3d', overwrite=True)
                            with open(temp_sdf_file, 'r') as f: sdf_content = f.read()
                            if not sdf_content:
                                st.warning("Downloaded 3D SDF file was empty.")
                                sdf_content = None
                        except pcp.NotFoundError:
                            st.warning(f"3D SDF not found on PubChem (CID {cid}).")
                        except Exception as e:
                            st.error(f"SDF Download Error: {e}")
                            sdf_content = None
                    else:
                        st.warning("No PubChem CID available to download 3D structure.")

                    if sdf_content:
                        try:
                            # 1. Create py3Dmol viewer
                            viewer = py3Dmol.view(width=450, height=400) # Slightly wider view
                            viewer.addModel(sdf_content, 'sdf')

                            # 2. Apply selected style
                            style_dict = {selected_style: {}}
                            if selected_style == 'cartoon':
                                style_dict['cartoon']['color'] = 'spectrum' # Common setting for cartoon
                            viewer.setStyle(style_dict)

                            viewer.setBackgroundColor('0xeeeeee')
                            viewer.zoomTo()

                            # 3. Render using stmol
                            stmol.showmol(viewer, height=400, width=450)

                        except Exception as e:
                            st.error(f"Error rendering 3D view with py3Dmol/stmol: {e}")
                            st.error(f"Exception type: {type(e)}")
                    # else: Handled by SDF download warnings

            else: # If mol_rdkit couldn't be created
                st.warning("Cannot display structures because molecule data (SMILES) is missing or invalid.")


        # --- Identifiers & Names Tab ---
        with tab_ids:
            st.subheader("Identifiers")
            cid = safe_get_attr(compound, 'cid')
            if cid != "N/A":
                st.markdown(f"**PubChem CID:** [{cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})")
            else:
                st.markdown("**PubChem CID:** N/A")

            st.write(f"**CAS:** {safe_get_attr(compound, 'cas')}")
            st.write(f"**InChI:** `{safe_get_attr(compound, 'inchi')}`")
            st.write(f"**InChIKey:** `{safe_get_attr(compound, 'inchikey')}`")

            st.subheader("Other Names (Synonyms)")
            synonyms = safe_get_attr(compound, 'synonyms')
            if synonyms != "N/A" and isinstance(synonyms, list) and len(synonyms) > 0:
                # Show first 5 directly
                st.write(", ".join(synonyms[:5]))
                # Put the rest in an expander if there are more
                if len(synonyms) > 5:
                    with st.expander("See all synonyms"):
                        st.json(synonyms)
            else:
                st.write("N/A")


    elif not error_occurred: # If compound is None but no exception was raised (e.g., Not Found handled in get_molecule_data)
        st.info("Molecule not found or no data retrieved.")

else: # Initial state when no input is given
    st.info("Enter an English molecule name above to start.")
