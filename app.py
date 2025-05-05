import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
import streamlit.components.v1 as components # Import Streamlit components
import re # For text manipulation (regex)
from PIL import Image # Needed by RDKit Draw sometimes explicitly, and Pillow is in requirements
import numpy # Import numpy - good practice although implicitly used

# --- Helper Functions ---

def normalize_name(name):
    """
    Cleans the molecule name for better search results.
    Removes leading/trailing whitespace, converts to lowercase.
    Handles spaces (converts multiple spaces to one).
    (Kept numeral conversion in case numbers are typed differently)
    """
    name = name.strip() # Remove leading/trailing whitespace
    # Convert Persian and Arabic numerals to English (just in case)
    persian_nums = "۰۱۲۳۴۵۶۷۸۹"
    arabic_nums = "٠١٢٣٤٥٦٧٨٩"
    english_nums = "0123456789"
    translation_table_persian = str.maketrans(persian_nums, english_nums)
    translation_table_arabic = str.maketrans(arabic_nums, english_nums)
    name = name.translate(translation_table_persian)
    name = name.translate(translation_table_arabic)

    # Replace multiple whitespace characters with a single space
    name = re.sub(r'\s+', ' ', name).strip()
    # Convert to lowercase (PubChem often works better with lowercase)
    name = name.lower()
    return name

def get_molecule_data(name):
    """
    Fetches molecule data from PubChem using PubChemPy.
    Includes basic error handling for PubChem queries.
    """
    try:
        # Search by name, limit to 1 result for simplicity
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0] # Return the first result
        else:
            return None
    except pcp.PubChemHTTPError as e:
        # Handle specific PubChem errors like not found or server issues more gracefully
        if "PUGREST.NotFound" in str(e):
             st.error(f"Molecule '{name}' not found in PubChem (PubChemHTTPError).")
        else:
             st.error(f"PubChem Error: {e}")
        return None
    except Exception as e:
        # Handle other potential errors (network connectivity, etc.)
        st.error(f"An unexpected error occurred during PubChem search: {e}")
        return None

# --- Streamlit User Interface ---

st.set_page_config(layout="wide") # Use wide layout for more space
st.title("🧪 Molecule Information Viewer ⌬")
st.markdown("""
Enter the English name of a molecule. The app will try to display its information and structure.
The search is case-insensitive and handles extra spaces (e.g., '1 chlorobutane' works).
""")

# --- Get User Input ---
raw_molecule_name = st.text_input("Molecule Name (English):", placeholder="e.g., Water, Aspirin, 1-chlorobutane")

if raw_molecule_name:
    # Normalize the input name directly
    normalized_name = normalize_name(raw_molecule_name)

    st.write(f"Searching for: `{normalized_name}`")

    # --- Search and Display Information ---
    compound = None # Initialize compound to None
    with st.spinner(f"Searching PubChem for '{normalized_name}'..."):
        compound = get_molecule_data(normalized_name)

    # Proceed only if compound data was successfully retrieved
    if compound:
        # Use IUPAC name if available, otherwise the raw input
        display_name = compound.iupac_name if compound.iupac_name else raw_molecule_name
        st.success(f"Found information for '{display_name}'!")

        # Display basic info
        st.subheader("Basic Information:")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Molecular Formula", compound.molecular_formula or "N/A")
        with col2:
            # CORRECTED Molecular Weight Handling
            mol_weight_display = "N/A" # Default value
            if hasattr(compound, 'molecular_weight') and compound.molecular_weight is not None:
                try:
                    weight_float = float(compound.molecular_weight)
                    mol_weight_display = f"{weight_float:.2f} g/mol"
                except (ValueError, TypeError):
                    mol_weight_display = str(compound.molecular_weight)
                    if "g/mol" not in mol_weight_display.lower():
                         mol_weight_display += " g/mol"
            st.metric("Molecular Weight", mol_weight_display)

        # Display CID and link to PubChem
        if compound.cid:
             st.markdown(f"**PubChem CID:** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")

        # --- Display Structures ---
        mol = None # RDKit molecule object
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

            # Display 2D structure
            with col_2d:
                st.markdown("**2D Structure:**")
                try:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                except Exception as e:
                    st.error(f"Error generating 2D image: {e}")

            # Display 3D structure using components.html
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                if compound.cid:
                    try:
                        st.write(f"[Debug] Attempting to download 3D SDF for CID: {compound.cid}") # Debug
                        temp_sdf_file = f'cid_{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        with open(temp_sdf_file, 'r') as f:
                            sdf_content = f.read()
                        if sdf_content:
                             st.write(f"[Debug] SDF content downloaded (length: {len(sdf_content)})") # Debug
                        else:
                             st.warning("[Debug] Downloaded 3D SDF file was empty.") # Debug
                             sdf_content = None # Ensure it's None if empty

                    except pcp.NotFoundError:
                        st.warning(f"3D structure (SDF) not found on PubChem for CID {compound.cid}.")
                    except Exception as e:
                        st.error(f"Error downloading 3D SDF from PubChem: {e}")
                else:
                    st.warning("Cannot fetch 3D structure without a PubChem CID.")

                # Render 3D view only if SDF content was obtained
                if sdf_content:
                    viewer_html = None # Initialize
                    try:
                         # Configure py3Dmol viewer
                        st.write("[Debug] Configuring py3Dmol viewer...") # Debug
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}})
                        viewer.setBackgroundColor('0xeeeeee')
                        viewer.zoomTo()
                        st.write("[Debug] Viewer configured. Generating HTML...") # Debug

                        # --- Generate HTML with specific error checking ---
                        try:
                            viewer_html = viewer.to_html()
                            # Check if HTML generation was successful
                            if viewer_html:
                                st.write(f"[Debug] HTML generated (type: {type(viewer_html)}, length: {len(viewer_html)})") # Debug
                            else:
                                st.error("[Debug] viewer.to_html() returned None or empty string.")
                                viewer_html = None # Ensure it's None if empty/failed

                        except TypeError as type_err:
                            st.error(f"[Debug] TypeError occurred specifically during viewer.to_html(): {type_err}") # Debug specific error
                            viewer_html = None # Ensure html is None on error
                        except Exception as html_gen_error:
                            st.error(f"[Debug] Non-TypeError exception during viewer.to_html(): {html_gen_error}") # Debug other errors
                            st.error(f"[Debug] Exception type: {type(html_gen_error)}")
                            viewer_html = None # Ensure html is None on error


                        # Embed the HTML using Streamlit components only if HTML was generated successfully
                        if viewer_html:
                            st.write("[Debug] Embedding HTML component...") # Debug
                            components.html(viewer_html, height=410, width=410)
                            st.write("[Debug] HTML component embedding attempted.") # Debug
                        # else: # No need for else, error already shown if viewer_html is None
                        #    st.error("Failed to generate HTML for 3D view.")


                    except Exception as e:
                        # Catch potential errors during py3Dmol setup (before to_html) or rendering
                        st.error(f"Error during py3Dmol setup or rendering: {e}")
                        st.error(f"Exception type: {type(e)}")

        # Display Synonyms (Other names)
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names (Synonyms):")
             st.json(compound.synonyms[:5])

    elif raw_molecule_name:
        # Error/Not found message displayed within get_molecule_data
        st.info("Search finished. See messages above for results or errors.")

else:
    # Initial message when the input box is empty
    st.info("Please enter an English molecule name in the box above.")
