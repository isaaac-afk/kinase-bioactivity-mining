"""
Kinase Bioactivity Mining

Queries the ChEMBL database for IC50 bioactivity records against a kinase
target, cleans the data, computes pIC50, and visualizes the potency
distribution. Outputs a sorted CSV and prints the top-10 most potent compounds.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from chembl_webresource_client.new_client import new_client


# ---- Configuration ----
TARGET_CHEMBL_ID = "CHEMBL203"
TARGET_NAME      = "EGFR"


def fetch_bioactivity(target_id: str) -> pd.DataFrame:
    print(f"Querying ChEMBL for IC50 records on {target_id} ...")
    activity = new_client.activity
    records = activity.filter(
        target_chembl_id=target_id,
        standard_type="IC50",
    ).only([
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_value",
        "standard_units",
        "standard_relation",
    ])
    df = pd.DataFrame(records)
    print(f"  Retrieved {len(df)} raw records")
    return df


def clean(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["standard_units"] == "nM"]
    df = df[df["standard_relation"] == "="]
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value", "canonical_smiles"])
    df = df[df["standard_value"] > 0]
    df["pIC50"] = -np.log10(df["standard_value"] * 1e-9)
    df = df.groupby("molecule_chembl_id", as_index=False).agg(
        canonical_smiles=("canonical_smiles", "first"),
        IC50_nM=("standard_value", "median"),
        pIC50=("pIC50", "median"),
    )
    print(f"  Cleaned: {len(df)} unique compounds with valid IC50 (nM)")
    return df


def plot_distribution(df: pd.DataFrame, target_name: str, out_path: str):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(df["pIC50"], bins=40, color="#1565C0",
            edgecolor="black", alpha=0.85)
    median = df["pIC50"].median()
    ax.axvline(median, color="red", ls="--",
               label=f"Median pIC50 = {median:.2f}")
    ax.set_xlabel("pIC50  =  -log10(IC50 in M)")
    ax.set_ylabel("Number of compounds")
    ax.set_title(f"{target_name} - inhibitor potency distribution "
                 f"({len(df):,} compounds from ChEMBL)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"  Saved plot: {out_path}")


def main():
    raw = fetch_bioactivity(TARGET_CHEMBL_ID)
    if raw.empty:
        print("No records returned. Check your target ID and internet connection.")
        return

    df = clean(raw)
    df = df.sort_values("pIC50", ascending=False).reset_index(drop=True)

    csv_path = f"{TARGET_NAME.lower()}_inhibitors.csv"
    plot_path = f"{TARGET_NAME.lower()}_pic50_distribution.png"

    df.to_csv(csv_path, index=False)
    print(f"  Saved table: {csv_path}")

    plot_distribution(df, TARGET_NAME, plot_path)

    print(f"\nTop 10 most potent {TARGET_NAME} inhibitors:")
    top10 = df.head(10)[["molecule_chembl_id", "IC50_nM", "pIC50"]]
    print(top10.to_string(index=False))


if __name__ == "__main__":
    main()