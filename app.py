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
            # --- CORRECTED Molecular Weight Handling ---
            mol_weight_display = "N/A" # Default value
            # Check if the attribute exists and is not None
            if hasattr(compound, 'molecular_weight') and compound.molecular_weight is not None:
                try:
                    # Attempt conversion to float BEFORE formatting
                    weight_float = float(compound.molecular_weight)
                    mol_weight_display = f"{weight_float:.2f} g/mol"
                except (ValueError, TypeError):
                    # If conversion fails, it might already be a string description
                    # Or it's a non-numeric value. Display it as is.
                    mol_weight_display = str(compound.molecular_weight)
                    # Optionally add units if it doesn't seem to have them
                    if "g/mol" not in mol_weight_display.lower():
                         mol_weight_display += " g/mol"


            # Display the result using st.metric
            st.metric("Molecular Weight", mol_weight_display)
            # --- End of Correction ---


        # Display CID and link to PubChem
        if compound.cid:
             st.markdown(f"**PubChem CID:** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")

        # --- Display Structures ---
        mol = None # RDKit molecule object
        # Check if SMILES data is available
        if hasattr(compound, 'isomeric_smiles') and compound.isomeric_smiles:
            try:
                # Create RDKit molecule object from SMILES string
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
                if mol is None: # Check if MolFromSmiles failed silently
                    st.warning("Could not generate molecule object from SMILES string.")

            except Exception as e:
                st.error(f"Error processing SMILES with RDKit: {e}")
                mol = None # Ensure mol is None if creation failed

        # Proceed only if the RDKit molecule object was created successfully
        if mol:
            st.subheader("Molecule Structure:")
            col_2d, col_3d = st.columns(2)

            # Display 2D structure
            with col_2d:
                st.markdown("**2D Structure:**")
                try:
                    # Generate 2D image using RDKit
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                except Exception as e:
                    st.error(f"Error generating 2D image: {e}")

            # Display 3D structure using components.html
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                # Attempt to download 3D SDF only if we have a CID
                if compound.cid:
                    try:
                        temp_sdf_file = f'cid_{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)
                        # Read the content of the downloaded SDF file
                        with open(temp_sdf_file, 'r') as f:
                            sdf_content = f.read()
                        if not sdf_content: # Check if file was empty
                             st.warning("Downloaded 3D SDF file was empty.")

                    except pcp.NotFoundError:
                        st.warning(f"3D structure (SDF) not found on PubChem for CID {compound.cid}.")
                    except Exception as e:
                        st.error(f"Error downloading 3D SDF from PubChem: {e}")
                else:
                    st.warning("Cannot fetch 3D structure without a PubChem CID.")

                # Render 3D view only if SDF content was obtained
                if sdf_content:
                    try:
                         # Configure py3Dmol viewer
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}}) # Display as stick model
                        viewer.setBackgroundColor('0xeeeeee') # Light gray background
                        viewer.zoomTo()

                        # Generate HTML representation using the correct method
                        viewer_html = viewer.to_html()

                        # Embed the HTML using Streamlit components
                        components.html(viewer_html, height=410, width=410)

                    except Exception as e:
                        # Catch potential errors during py3Dmol rendering
                        st.error(f"Error rendering 3D view with py3Dmol/components: {e}")
                        # Also print the type of exception for debugging
                        st.error(f"Exception type: {type(e)}")

        # Display Synonyms (Other names)
        # Check if compound and synonyms attribute exist and synonyms list is not empty
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names (Synonyms):")
             # Show only the first 5 synonyms to avoid clutter
             st.json(compound.synonyms[:5])

    # This message is shown only if a search was attempted (raw_molecule_name is not empty)
    # but get_molecule_data returned None (or failed)
    elif raw_molecule_name:
        # The error message should have been displayed within get_molecule_data
        # Add a generic message here if needed, but avoid redundancy
        st.info("Search finished. See messages above for results or errors.")


else:
    # Initial message when the input box is empty
    st.info("Please enter an English molecule name in the box above.")
