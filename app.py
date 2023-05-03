
import pandas as pd
import streamlit as st
import rdkit
from rdkit.Chem import Draw
from rdkit import Chem

with st.sidebar:
    default_smiles = "\n".join(["CN2C(=O)N(C)C(=O)C1=C2N=CN1C",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"])
    text = st.text_area("Enter SMILES here", height=100, key="smiles", value=default_smiles)
    lookup_text = st.text_area("Enter Query SMILES here", height=100, key="query_smiles", value=default_smiles)

    smiles = [x.strip() for x in text.split("\n")]
    query_smiles = [x.strip() for x in lookup_text.split("\n")]

    canon_smiles = [Chem.CanonSmiles(smile) for smile in smiles if smile]
    canon_query_smiles = [Chem.CanonSmiles(smile) for smile in query_smiles if smile]

mols = [Chem.MolFromSmiles(smile) for smile in canon_smiles]
imgs = [Draw.MolToImage(mol, size=(600,400)) for mol in mols]

for i, (img, smi, csmi) in enumerate(zip(imgs, smiles, canon_smiles)):
    if smi != csmi:
        caption = f"Original: {smi}\n\nCanonical: {csmi}"
    else:
        caption = smi

    with st.expander(f"{i+1}: {smi}"):
        st.markdown(f"```\n{caption}\n```")
        st.image(img, width=600)

match_positions = []
for x in canon_smiles:
    try:
        match_positions.append(canon_query_smiles.index(x))
    except ValueError:
        match_positions.append(None)


df = pd.DataFrame({
    'Input SMILES': smiles,
    'Canonical SMILES': canon_smiles,
    'is cannonical': [x == y for x, y in zip(smiles, canon_smiles)],
    'query position': match_positions,
})

st.dataframe(df)
# add a button to download the dataframe as a csv
st.download_button(
    label="Download data as CSV (you can open it in excel)",
    data=df.to_csv().encode("utf-8"),
    file_name="canonical_smiles.csv",
    mime="text/csv",
)