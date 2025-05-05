import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
# Removed: from streamlit_py3dmol import st_py3dmol
import streamlit.components.v1 as components # Import Streamlit components
import re # For text manipulation (regex)
from PIL import Image # Needed by RDKit Draw sometimes explicitly, and Pillow is in requirements

# --- Helper Functions ---

def normalize_name(name):
    """
    Cleans the molecule name for better search results.
    Removes leading/trailing whitespace, converts to lowercase.
    Converts Persian/Arabic numerals to English numerals.
    Handles spaces (converts multiple spaces to one).
    """
    name = name.strip() # Remove leading/trailing whitespace
    # Convert Persian and Arabic numerals to English
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
    except Exception as e:
        st.error(f"Error connecting to PubChem: {e}")
        return None

def detect_language(text):
    """
    Simple detection of Persian language based on character range.
    """
    # Check if any character falls within the Persian Unicode range
    if re.search(r'[\u0600-\u06FF]', text):
        return "Persian"
    else:
        # Assume English or other non-Persian if no Persian chars found
        return "English"

# --- Streamlit User Interface ---

st.set_page_config(layout="wide") # Use wide layout for more space
st.title("🧪 Molecule Information Viewer ⌬")
st.markdown("""
Enter the name of a molecule in English or Persian. The app will try to display its information and structure.
**Note:** Searching by the exact English name usually yields better results. Persian names might require translation (not implemented yet).
""")

# --- Get User Input ---
raw_molecule_name = st.text_input("Molecule Name:", placeholder="e.g., Water, Aspirin, 1-chlorobutane, آب")

if raw_molecule_name:
    # Detect input language
    lang = detect_language(raw_molecule_name)
    st.write(f"Detected language: {'Persian' if lang == 'Persian' else 'English/Other'}")

    if lang == "Persian":
        st.warning("⚠️ Note: Searching with Persian names may be inaccurate or return no results. Using the English name is recommended.")
        # Future improvement: Add a translation service here.
        # For now, just try the normalized Persian name (low chance of success on PubChem)
        normalized_name = normalize_name(raw_molecule_name)
    else:
        # Normalize the English name
        normalized_name = normalize_name(raw_molecule_name)

    st.write(f"Searching for: `{normalized_name}`")

    # --- Search and Display Information ---
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
            # Format weight to 2 decimal places if available
            mol_weight = f"{compound.molecular_weight:.2f} g/mol" if compound.molecular_weight else "N/A"
            st.metric("Molecular Weight", mol_weight)

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

            # Display 3D structure
            with col_3d:
                st.markdown("**3D Structure (Interactive):**")
                sdf_content = None
                try:
                    # Attempt to download 3D SDF format from PubChem
                    temp_sdf_file = f'cid_{compound.cid}_3d.sdf'
                    pcp.download('SDF', temp_sdf_file, compound.cid, 'cid', record_type='3d', overwrite=True)

                    # Read the content of the downloaded SDF file
                    with open(temp_sdf_file, 'r') as f:
                       sdf_content = f.read()

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

                        # --- NEW WAY to display using components.html ---
                        # Generate HTML representation from py3Dmol viewer
                        viewer_html = viewer._repr_html_()
                        # Embed the HTML using Streamlit components
                        components.html(viewer_html, height=410, width=410) # Slightly larger height/width can help avoid scrollbars
                        # --- End of new way ---

                    except Exception as e:
                        st.error(f"Error rendering 3D view with py3Dmol/components: {e}")

                # else: (Handled by the warning above if sdf_content is None due to NotFoundError)
                     # st.warning("Could not display 3D structure (SDF content missing or download failed).")


        else:
            # This message appears if SMILES string was missing or RDKit failed to parse it
            st.warning("Molecule structure data (SMILES) not available or couldn't be processed.")

        # Display Synonyms (Other names)
        if compound.synonyms:
             st.subheader("Other Names (Synonyms):")
             # Show only the first 5 synonyms to avoid clutter
             st.json(compound.synonyms[:5])


    elif raw_molecule_name: # Only show error if user actually typed something and nothing was found
        st.error(f"Sorry, no molecule found matching '{raw_molecule_name}' (or normalized name `{normalized_name}`).")
        st.info("Please check the spelling. Using the English name is recommended.")

else:
    # Initial message when the input box is empty
    st.info("Please enter a molecule name in the box above.")
