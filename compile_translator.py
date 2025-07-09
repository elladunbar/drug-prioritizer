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

    def filter_equal(self, column: str | Sequence[str], value: Any, *, in_place: bool = False) -> pd.DataFrame:
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

    def get_drug_list(self) -> pd.DataFrame:
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
        df = df.groupby("result_id").agg(lambda x: set(x)).reset_index()

        # only valid drugs
        df = df[~df["result_id"].str.startswith(self._invalid_id_types)]

        # annotate id source more plainly
        df["id_type"] = df["result_id"].str.split(":").str[0].map(self._valid_id_to_drugbank)
        df["id_type"] = df["id_type"].fillna("InChIKey")

        # chop off id type
        needs_id_removal = df["result_id"].str.contains(":", regex=False)
        df["result_id"][needs_id_removal] = df["result_id"][needs_id_removal].str.split(":").str[1]

        # add current date and time
        df["date_time"] = datetime.datetime.now()

        return df


if __name__ == "__main__":
    translator_folder = "data/translator/"
    translator_results = TranslatorData(translator_folder)
    drug_list = translator_results.get_drug_list()
    with open("data/translator_input.json", "w") as f:
        drug_list.to_json(f, orient="records")
