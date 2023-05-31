import numpy as np
import pandas as pd
import streamlit as st
import rdkit
from rdkit.Chem import Draw
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import altair as alt

with st.sidebar:
    default_smiles = "\n".join(
        ["CN2C(=O)N(C)C(=O)C1=C2N=CN1C", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"]
    )
    text = st.text_area(
        "Enter SMILES here", height=100, key="smiles", value=default_smiles
    )
    lookup_text = st.text_area(
        "Enter Query SMILES here", height=100, key="query_smiles", value=default_smiles
    )

    smiles = [x.strip() for x in text.split("\n")]
    smiles = [x for x in smiles if x]
    query_smiles = [x.strip() for x in lookup_text.split("\n")]
    query_smiles = [x for x in query_smiles if x]

    canon_smiles = [Chem.CanonSmiles(smile) for smile in smiles if smile]
    canon_query_smiles = [Chem.CanonSmiles(smile) for smile in query_smiles if smile]

mols = [Chem.MolFromSmiles(smile) for smile in canon_smiles]
query_mols = [Chem.MolFromSmiles(smile) for smile in canon_query_smiles]
imgs = [Draw.MolToImage(mol, size=(600, 400)) for mol in mols]

st.markdown("# Input SMILES (and structures)")
for i, (img, smi, csmi) in enumerate(zip(imgs, smiles, canon_smiles)):
    if smi != csmi:
        caption = f"Original: {smi}\n\nCanonical: {csmi}"
    else:
        caption = smi

    with st.expander(f"{i+1}: {smi}"):
        st.markdown(f"```\n{caption}\n```")
        st.image(img, width=600)

match_positions = [None for _ in range(len(canon_smiles))]
for i, x in enumerate(canon_smiles):
    try:
        match_positions[i] = canon_query_smiles.index(x)
    except ValueError:
        match_positions[i] = None


target_fps = [FingerprintMols.FingerprintMol(x) for x in mols]
mws = [ExactMolWt(x) for x in mols]
query_fps = [FingerprintMols.FingerprintMol(x) for x in query_mols]

similarity_matrix = [
    [DataStructs.FingerprintSimilarity(x, y) for x in query_fps] for y in target_fps
]
closest_match = [np.argmax(x) for x in similarity_matrix]

st.markdown("# Results Table")
df = pd.DataFrame(
    {
        "Input SMILES": smiles,
        "MW": mws,
        "Closest match": [query_smiles[i] for i in closest_match],
        "Canonical SMILES": canon_smiles,
        "is cannonical": [x == y for x, y in zip(smiles, canon_smiles)],
        "query position": match_positions,
        "closest query position": closest_match,
        "similarity": [x[i] for x, i in zip(similarity_matrix, closest_match)],
    }
)
st.dataframe(df, use_container_width=True)
# add a button to download the dataframe as a csv
st.download_button(
    label="Download data as CSV (you can open it in excel)",
    data=df.to_csv().encode("utf-8"),
    file_name="canonical_smiles.csv",
    mime="text/csv",
)

st.markdown("# Similarity Matrix")

x, y = np.meshgrid(range(1, len(query_fps) + 1), range(1, len(target_fps) + 1))
sim_array = np.array(similarity_matrix)
chart_df = pd.DataFrame(
    {"target": x.ravel(), "query": y.ravel(), "similarity": sim_array.ravel()}
)
chart = (
    alt.Chart(chart_df)
    .mark_rect()
    .encode(x="target:O", y="query:O", color="similarity:Q")
)

st.altair_chart(chart, use_container_width=True)
