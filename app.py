import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
from st_py3dmol import st_py3dmol # برای نمایش سه بعدی در استریملیت
import re # برای کار با متن

# --- توابع کمکی ---

def normalize_name(name):
    """
    اسم مولکول رو تمیز می‌کنه تا جستجو بهتر انجام بشه.
    فاصله‌ها رو حذف می‌کنه، به حروف کوچک تبدیل می‌کنه.
    اعداد فارسی رو به انگلیسی تبدیل می‌کنه.
    """
    name = name.strip() # حذف فاصله‌های اضافی اول و آخر
    # تبدیل اعداد فارسی و عربی به انگلیسی
    persian_nums = "۰۱۲۳۴۵۶۷۸۹"
    arabic_nums = "٠١٢٣٤٥٦٧٨٩"
    english_nums = "0123456789"
    translation_table_persian = str.maketrans(persian_nums, english_nums)
    translation_table_arabic = str.maketrans(arabic_nums, english_nums)
    name = name.translate(translation_table_persian)
    name = name.translate(translation_table_arabic)

    # حذف همه فاصله‌ها (بحث‌برانگیز، شاید بهتر باشه فقط فاصله‌های اضافی حذف شه)
    # name = name.replace(" ", "") # روش ساده ولی ممکنه مشکل‌ساز باشه
    # روش بهتر: حذف فاصله‌های تکراری و تبدیل به یک فاصله
    name = re.sub(r'\s+', ' ', name).strip()
    # به حروف کوچک تبدیل کن (PubChem معمولا با حروف کوچک بهتر کار می‌کنه)
    name = name.lower()
    return name

def get_molecule_data(name):
    """
    با استفاده از PubChemPy اطلاعات مولکول رو پیدا می‌کنه.
    """
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0] # اولین نتیجه رو برمی‌گردونه
        else:
            return None
    except Exception as e:
        st.error(f"خطا در ارتباط با PubChem: {e}")
        return None

def detect_language(text):
    """
    تشخیص ساده زبان فارسی بر اساس وجود حروف فارسی.
    """
    if re.search(r'[\u0600-\u06FF]', text):
        return "Persian"
    else:
        # فرض می‌کنیم بقیه انگلیسی هستن (یا حداقل قابل جستجو در PubChem)
        return "English"

# --- رابط کاربری Streamlit ---

st.set_page_config(layout="wide") # صفحه رو عریض‌تر می‌کنه
st.title("🧪 اپلیکیشن نمایش اطلاعات مولکول ⌬")
st.markdown("""
اسم یک مولکول را به فارسی یا انگلیسی وارد کنید. اپلیکیشن سعی می‌کند اطلاعات و ساختار آن را نمایش دهد.
**نکته:** جستجو بر اساس نام انگلیسی دقیق‌تر است. برای نام‌های فارسی، ممکن است نیاز به ترجمه باشد.
""")

# --- دریافت ورودی از کاربر ---
raw_molecule_name = st.text_input("نام مولکول:", placeholder="مثلا: آب، Water, ۱-کلروبوتان, 1-chlorobutane")

if raw_molecule_name:
    # تشخیص زبان ورودی
    lang = detect_language(raw_molecule_name)
    st.write(f"زبان تشخیص داده شده: {'فارسی' if lang == 'Persian' else 'انگلیسی/غیره'}")

    if lang == "Persian":
        st.warning("⚠️ توجه: جستجو با نام فارسی ممکن است دقیق نباشد یا نتیجه‌ای نداشته باشد. بهتر است نام انگلیسی مولکول را وارد کنید.")
        # در آینده می‌توان اینجا یک سرویس ترجمه اضافه کرد
        # فعلا، فقط نام فارسی نرمال‌شده را امتحان می‌کنیم (احتمال موفقیت کم است)
        normalized_name = normalize_name(raw_molecule_name)
    else:
        # نرمال‌سازی نام انگلیسی
        normalized_name = normalize_name(raw_molecule_name)

    st.write(f"جستجو برای: `{normalized_name}`")

    # --- جستجو و نمایش اطلاعات ---
    with st.spinner(f"در حال جستجوی اطلاعات برای {raw_molecule_name}..."):
        compound = get_molecule_data(normalized_name)

    if compound:
        st.success(f"اطلاعات مولکول '{compound.iupac_name or raw_molecule_name}' پیدا شد!")

        # نمایش اطلاعات پایه
        st.subheader("اطلاعات پایه:")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("فرمول مولکولی", compound.molecular_formula or "N/A")
        with col2:
            st.metric("وزن مولکولی", f"{compound.molecular_weight or 0:.2f} g/mol")

        # نمایش CID و لینک به PubChem
        if compound.cid:
             st.markdown(f"**شناسه PubChem (CID):** [{compound.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid})")


        # --- نمایش ساختارها ---
        mol = None
        if compound.isomeric_smiles:
            try:
                mol = Chem.MolFromSmiles(compound.isomeric_smiles)
            except Exception as e:
                st.error(f"خطا در پردازش SMILES با RDKit: {e}")

        if mol:
            st.subheader("ساختار مولکول:")
            col_2d, col_3d = st.columns(2)

            # نمایش ساختار 2 بعدی
            with col_2d:
                st.markdown("**ساختار دو بعدی:**")
                try:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img)
                except Exception as e:
                    st.error(f"خطا در تولید تصویر دو بعدی: {e}")

            # نمایش ساختار 3 بعدی
            with col_3d:
                st.markdown("**ساختار سه بعدی (تعاملی):**")
                try:
                    # گرفتن فرمت SDF از PubChem برای اطلاعات سه بعدی
                    sdf_3d = pcp.download('SDF', f'cid_{compound.cid}_3d.sdf', compound.cid, 'cid', record_type='3d', overwrite=True)
                    # خواندن محتوای فایل SDF دانلود شده
                    with open(f'cid_{compound.cid}_3d.sdf', 'r') as f:
                       sdf_content = f.read()

                    if sdf_content:
                         # تنظیمات نمایشگر py3Dmol
                        viewer = py3Dmol.view(width=400, height=400)
                        viewer.addModel(sdf_content, 'sdf')
                        viewer.setStyle({'stick': {}}) # نمایش به صورت میله‌ای
                        viewer.zoomTo()
                        # نمایش در Streamlit با استفاده از کامپوننت سفارشی
                        st_py3dmol(viewer)
                    else:
                         st.warning("اطلاعات ساختار سه بعدی در دسترس نیست.")

                except Exception as e:
                    # اگر دانلود SDF 3D شکست خورد یا py3Dmol مشکل داشت
                    st.error(f"خطا در نمایش سه بعدی: {e}")
                    st.info("ممکن است ساختار سه بعدی برای این مولکول در PubChem موجود نباشد یا خطایی رخ داده باشد.")

        else:
            st.warning("ساختار مولکولی (SMILES) برای پردازش در دسترس نیست.")

        # نمایش اطلاعات بیشتر (Synonyms)
        if compound.synonyms:
             st.subheader("نام‌های دیگر:")
             # نمایش ۵ نام اول برای جلوگیری از شلوغی
             st.json(compound.synonyms[:5])


    elif raw_molecule_name: # فقط اگر کاربر چیزی وارد کرده بود و نتیجه‌ای نبود
        st.error(f"متاسفانه مولکولی با نام '{raw_molecule_name}' (یا معادل نرمال‌شده آن `{normalized_name}`) پیدا نشد.")
        st.info("لطفاً از صحت نام مطمئن شوید و ترجیحاً از نام انگلیسی استفاده کنید.")

else:
    st.info("لطفاً نام یک مولکول را در کادر بالا وارد کنید.")
