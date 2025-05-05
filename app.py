import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
import streamlit.components.v1 as components # Import Streamlit components
import re # For text manipulation (regex)
from PIL import Image # Needed by RDKit Draw sometimes explicitly, and Pillow is in requirements

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
    """
    try:
        # Search by name, limit to 1 result for simplicity
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0] # Return the first result
        else:
            return None
    except pcp.PubChemHTTPError as e:
        # Handle specific PubChem errors like not found or server issues
        st.error(f"PubChem Error: {e}")
        return None
    except Exception as e:
        # Handle other potential errors (network, etc.)
        st.error(f"An unexpected error occurred: {e}")
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
    with st.spinner(f"Searching for information on {raw_molecule_name}..."):
        compound = get_molecule_data(normalized_name)

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
            # Molecular Weight Handling
            mol_weight_display = "N/A" # Default value
            if compound.molecular_weight is not None: # Check if it exists
                try:
                    # Try converting to float before formatting
                    weight_float = float(compound.molecular_weight)
                    mol_weight_display = f"{weight_float:.2f} g/mol"
                except (ValueError, TypeError):
                    # If conversion fails, display as is (might be a string already)
                    mol_weight_display = f"{compound.molecular_weight} g/mol" # Fallback

            st.metric("Molecular Weight", mol_weight_display)

        # Display CID and link to PubChem
        if compound.cid:
             st.markdown(f"**PubChem CID:** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")

        # --- Display Structures ---
        mol = None
        if compound.isomeric_smiles:
            try:
                # Create RDKit molecule object from SMILES string
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
            except Exception as e:
                st.error(f"Error processing SMILES with RDKit: {e}")
                # Don't proceed if mol object creation failed
                mol = None

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
                try:
                    # Attempt to download 3D SDF format from PubChem
                    # Use CID for reliable downloading
                    if compound.cid:
                        temp_sdf_file = f'cid_{compound.cid}_3d.sdf'
                        pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)

                        # Read the content of the downloaded SDF file
                        with open(temp_sdf_file, 'r') as f:
                            sdf_content = f.read()
                    else:
                         st.warning("Cannot fetch 3D structure without a PubChem CID.")


                except pcp.NotFoundError:
                     st.warning("3D structure (SDF) not found on PubChem for this compound.")
                except Exception as e:
                    st.error(f"Error downloading 3D SDF from PubChem: {e}")

                if sdf_content:
                    try:
                         # Configure py3Dmol viewer
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}}) # Display as stick model
                        viewer.setBackgroundColor('0xeeeeee') # Light gray background
                        viewer.zoomTo()

                        # --- CORRECT WAY to get HTML for embedding ---
                        # Generate HTML representation using the correct method
                        viewer_html = viewer.to_html()
                        # --- End of correction ---

                        # Embed the HTML using Streamlit components
                        components.html(viewer_html, height=410, width=410)

                    except Exception as e:
                        st.error(f"Error rendering 3D view with py3Dmol/components: {e}")
                # Added a condition for when SDF couldn't be fetched/read but no specific error was raised before
                elif compound.cid: # Only show this if we attempted download (i.e., had CID)
                    st.info("Could not display 3D structure (SDF content missing or download failed).")

        # Display Synonyms (Other names)
        # Check if compound and synonyms attribute exist
        if compound and hasattr(compound, 'synonyms') and compound.synonyms:
             st.subheader("Other Names (Synonyms):")
             # Show only the first 5 synonyms to avoid clutter
             st.json(compound.synonyms[:5])

    # This message is shown only if a search was attempted (raw_molecule_name is not empty)
    # but get_molecule_data returned None
    elif raw_molecule_name:
        st.error(f"Sorry, no molecule found matching '{raw_molecule_name}' (or normalized name `{normalized_name}`).")
        st.info("Please check the spelling and use the English name.")

else:
    # Initial message when the input box is empty
    st.info("Please enter an English molecule name in the box above.")
