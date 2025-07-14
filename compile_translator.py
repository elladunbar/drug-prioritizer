#!/usr/bin/env uv run

import datetime
import os
from collections.abc import Sequence
from typing import Any

import pandas as pd


class TranslatorData:
    _drug_categories = [f"biolink:{biolink_cat}" for biolink_cat in ("Drug", "SmallMolecule")]
    _invalid_id_types = ("CHV:", "MESH:", "MONDO:", "UMLS:")
    _valid_id_to_drugbank = {
        "CHEBI": "ChEBI",
        "CHEMBL.COMPOUND": "ChEMBL",
        "DRUGBANK": "drugbank-id",
        "PUBCHEM.COMPOUND": "PubChem Compound",
        "PUBCHEM.SUBSTANCE": "PubChem Substance",
        "UNII": "unii",
    }

    def __init__(self, folder: str):
        files = tuple(filter(lambda s: s.endswith(".csv"), os.listdir(folder)))

        dfs = []
        for file in files:
            with open(folder + file, "r") as f:
                df = pd.read_csv(f, low_memory=False)

            df["search term"] = file.split("_")[0].lower()  # downloader names = query-name_other-stuff
            dfs.append(df)

        self.data = pd.concat(dfs, ignore_index=True)

    def filter_equal(self, column: str | Sequence[str], value: Any, *, in_place: bool = False):
        if isinstance(column, str):
            column = (column,)

        mask = pd.Series(False, index=self.data.index)
        for col in column:
            if isinstance(value, Sequence) and not isinstance(value, str):
                mask |= self.data[col].isin(value)
            else:
                mask |= self.data[col] == value

        df = self.data[mask]
        if in_place:
            self.data = df
        return df

    def _aggregate(self, df: pd.DataFrame):
        agg_funcs = {}
        grouped_df = df.groupby("result_id")

        for col in df.columns:
            if col == "result_id":
                continue

            unqique_vals_per_group = grouped_df[col].nunique()
            if (unqique_vals_per_group <= 1).all():
                agg_funcs[col] = lambda x: x.iloc[0]
            else:
                agg_funcs[col] = lambda x: set(x)

        return grouped_df.agg(agg_funcs)

    def get_drug_list(self):
        df = self.filter_equal(
            ("result_subjectNode_cat", "result_objectNode_cat"),
            self._drug_categories,
        ).copy(deep=True)

        # combine subject and object to get chemicals only
        subject_is_drug = df["result_subjectNode_cat"].isin(self._drug_categories)
        for col_ending in ("name", "id", "cat"):
            df[f"result_{col_ending}"] = df[f"result_subjectNode_{col_ending}"].where(
                subject_is_drug, df[f"result_objectNode_{col_ending}"]
            )

        # combine repeat drugs
        df = self._aggregate(df).reset_index()

        # only valid drugs
        df = df[~df["result_id"].str.startswith(self._invalid_id_types)]

        # annotate id source more plainly
        df["id_type"] = df["result_id"].str.split(":").str[0].map(self._valid_id_to_drugbank)
        df["id_type"] = df["id_type"].fillna("InChIKey")

        # chop off id type
        df["result_id"] = df["result_id"].str.replace(r".*:(.*)$", r"\1", regex=True)

        # add current date and time
        df["date_time"] = datetime.datetime.now()

        # get relevant columns
        cols = ["result_name", "result_id", "id_type", "search term", "id_type", "date_time"]
        df = df[cols]

        return df


if __name__ == "__main__":
    translator_folder = "data/translator/"
    translator_results = TranslatorData(translator_folder)
    drug_list = translator_results.get_drug_list()
    print(drug_list)

    # db = database.Database()
    # db.import_dataframe()
