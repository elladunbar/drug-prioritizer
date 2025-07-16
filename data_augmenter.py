import xml.etree.ElementTree as ET
from collections.abc import Iterable
from typing import Callable

import catboost as cb
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from tqdm import tqdm

import maplight_gnn


class DrugBank:
    namespaces = {"": "http://www.drugbank.ca"}

    def _execute_if_found(self, element: ET.Element, name: str, on_found: Callable):
        result = element.find(name, self.namespaces)
        self._if_not_none(result, on_found)

    def _if_not_none(self, value, on_present: Callable):
        if value is not None:
            on_present(value)

    def __init__(self, filename: str, ids, id_types, names):
        self.filename = filename
        self.ids = ids
        self.id_types = id_types
        self.names = names.str.lower()
        self.get_ids = {
            "cas-number": self.cas_number,
            "ChEBI": self.chebi,
            "ChEMBL": self.chembl,
            "drugbank-id": self.drugbank,
            "InChIKey": self.inchikey,
            "PubChem Compound": self.pubchem_compound,
            "PubChem Substance": self.pubchem_substance,
            "unii": self.unii,
        }

    def get_matches(self):
        for _, element in tqdm(ET.iterparse(self.filename, ["end"])):
            # don't care about non-drug entries
            if element.tag[24:None:None] != "drug":
                continue

            matches = self.check_match(element)

            # make sure there are matches before doing more work
            if not matches.any():
                continue

            yield (matches, element)

    def check_match(self, element):
        matches = pd.Series(False, index=self.ids.index)
        for id_type, id_func in self.get_ids.items():
            id_val = id_func(element)
            if id_val is None:
                continue
            id_matches = (self.id_types == id_type) & (self.ids == id_val)
            matches = matches | id_matches

        # names can't use same logic as other id types
        generic_names, brand_names = self.all_names(element)
        matches = matches | self.names.isin(generic_names)
        matches = matches | self.names.isin(brand_names)
        return matches

    def all_names(self, element):
        generic_names = set()
        brand_names = set()
        main_name = self.name(element)
        generic_names.add(main_name.lower()) if main_name is not None else None
        it = element.find("synonyms", self.namespaces)
        if it is not None:
            for synonym in it.iter():
                generic_names.add(synonym.text.lower()) if synonym is not None and synonym.text is not None else None

        it = element.find("products", self.namespaces)
        if it is not None:
            for product in it.iter():
                brand_name = product.find("name", self.namespaces)
                brand_names.add(brand_name.text.lower()) if brand_name is not None else None

        generic_names = tuple(filter(lambda s: "\n" not in s, generic_names))
        brand_names = tuple(filter(lambda s: "\n" not in s, brand_names))
        return (generic_names, brand_names)

    def cas_number(self, element):
        it = element.find("cas-number", self.namespaces)
        return it.text if it is not None else None

    def chebi(self, element):
        return self.from_external_identifiers(element, "ChEBI")

    def chembl(self, element):
        return self.from_external_identifiers(element, "ChEMBL")

    def drugbank(self, element):
        it = element.find("drugbank-id", self.namespaces)
        return it.text if it is not None else None

    def fda_approval(self, element):
        it = element.find("groups", self.namespaces)
        return "approved" in tuple(it.itertext()) if it is not None else None

    def inchikey(self, element):
        return self.from_calculated_properties(element, "InChIKey")

    def indication(self, element):
        it = element.find("indication", self.namespaces)
        return it.text if it is not None else None

    def mechanism(self, element):
        it = element.find("mechanism-of-action", self.namespaces)
        return it.text if it is not None else None

    def name(self, element):
        it = element.find("name", self.namespaces)
        return it.text if it is not None else None

    def prices(self, element):
        it = element.find("prices", self.namespaces)
        if it is not None:
            prices = list()
            for price_element in it.iterfind("price", self.namespaces):
                price = price_element.find("cost", self.namespaces)
                prices.append(price.text + price.attrib.get("currency")) if price is not None else None
            return prices
        else:
            return None

    def pubchem_compound(self, element):
        return self.from_external_identifiers(element, "PubChem Compound")

    def pubchem_substance(self, element):
        return self.from_external_identifiers(element, "PubChem Substance")

    def smiles(self, element):
        return self.from_calculated_properties(element, "SMILES")

    def unii(self, element):
        it = element.find("unii", self.namespaces)
        return it.text if it is not None else None

    def from_external_identifiers(self, element, resource_type):
        it = element.find("external-identifiers", self.namespaces)
        if it is not None:
            for external_identifier in it.iterfind("external-identifier", self.namespaces):
                if external_identifier.findtext("resource", namespaces=self.namespaces) == resource_type:
                    return external_identifier.findtext("identifier", namespaces=self.namespaces)
        else:
            return None

    def from_calculated_properties(self, element, kind_type):
        it = element.find("calculated-properties", self.namespaces)
        if it is not None:
            for property in it.iterfind("property", self.namespaces):
                if property.findtext("kind", namespaces=self.namespaces) == kind_type:
                    return property.findtext("value", namespaces=self.namespaces)
        else:
            return None


class DataAugmenter:
    def __init__(self, filename):
        self.filename = filename
        self.drug_list = None
        self.admet_models = None

    def load_drug_queries(self):
        if self.filename.endswith(".csv"):
            with open(self.filename, "r") as f:
                self.drug_list = pd.read_csv(f)
        elif self.filename.endswith(".json"):
            with open(self.filename, "r") as f:
                self.drug_list = pd.read_json(f, orient="records")
        else:
            raise ValueError("Data file must be .csv or .json")
        return self

    def load_admet_models(self, models):
        self.admet_models = dict()
        for name, path in models.items():
            model = cb.CatBoostClassifier()
            model.load_model(path)
            self.admet_models[name] = model
        return self

    def save_drug_info(self, filename):
        if self.drug_list is None:
            raise ValueError("drug-list must be loaded first.")
        with open(filename, "w") as f:
            return self.drug_list.to_json(f, orient="records")

    def _unwrap_list(self, x):
        return x[0] if isinstance(x, list) else x

    def _is_collection(self, x):
        return isinstance(x, Iterable) and not isinstance(x, (str, bytes))

    def _flatten(self, coll):
        if not self._is_collection(coll):
            return [coll]
        result = []
        for element in coll:
            result.extend(self._flatten(element))
        return result

    def match_drugbank(self, filename, id_col_name, id_type_col_name, name_col_name):
        if self.drug_list is None:
            raise ValueError("drug-list is not defined. Call load-drug-queries before match-drugbank.")

        # make sure the cols are strings and not lists of strings
        id_col = self.drug_list[id_col_name].apply(self._unwrap_list)
        id_type_col = self.drug_list[id_type_col_name].apply(self._unwrap_list)
        name_col = self.drug_list[name_col_name].apply(self._unwrap_list)

        # tedius column making
        all_names_column = "All Names"
        self.drug_list[all_names_column] = self.drug_list.apply(lambda _: list(), axis=1)
        cas_column = "CAS Registry Number"
        self.drug_list[cas_column] = None
        fda_column = "FDA Approved"
        self.drug_list[fda_column] = None
        indication_column = "Indication"
        self.drug_list[indication_column] = None
        mechanism_column = "Mechanism"
        self.drug_list[mechanism_column] = None
        name_column = "DrugBank Name"
        self.drug_list[name_column] = None
        price_column = "Prices"
        self.drug_list[price_column] = self.drug_list.apply(lambda _: list(), axis=1)
        smiles_column = "SMILES"
        self.drug_list[smiles_column] = None
        unii_column = "UNII"
        self.drug_list[unii_column] = None
        drugbank = DrugBank(filename, id_col, id_type_col, name_col)
        for matches, element in drugbank.get_matches():
            self.drug_list.loc[matches, all_names_column] = self.drug_list.loc[matches, all_names_column].apply(
                lambda _: self._flatten(drugbank.all_names(element))
            )
            self.drug_list.loc[matches, cas_column] = drugbank.cas_number(element)
            self.drug_list.loc[matches, fda_column] = drugbank.fda_approval(element)
            self.drug_list.loc[matches, indication_column] = drugbank.indication(element)
            self.drug_list.loc[matches, mechanism_column] = drugbank.mechanism(element)
            self.drug_list.loc[matches, name_column] = drugbank.name(element)
            self.drug_list.loc[matches, price_column] = self.drug_list.loc[matches, price_column].apply(
                lambda _: drugbank.prices(element)
            )
            self.drug_list.loc[matches, smiles_column] = drugbank.smiles(element)
            self.drug_list.loc[matches, unii_column] = drugbank.unii(element)

    def deduplicate(self):
        if self.drug_list is None:
            raise ValueError("drug-list is not defined. Call load-drug-queries before deduplicate.")
        if "DrugBank Name" not in self.drug_list.columns:
            raise ValueError("ID data does not exist yet. Run match-drugbank to create it.")

        def deduplication_aggregator(x):
            y = []
            for item in x:
                y.extend(item) if isinstance(item, list) else y.append(item)
            z = set(y)
            z.discard(None)
            return None if len(z) == 0 else z.pop() if len(z) == 1 else z if True else None

        self.drug_list = self.drug_list.groupby("DrugBank Name").agg(deduplication_aggregator).reset_index()

    def predict_admet(self):
        if self.drug_list is None:
            raise ValueError("drug-list is not defined. Call load-drug-queries before predict-admet.")
        if self.admet_models is None:
            raise ValueError("admet-models is not defined. Call load-admet-models before predict-admet.")
        if "SMILES" not in self.drug_list.columns:
            raise ValueError("SMILES data does not exist yet. Run match-drugbank to create it.")

        RDLogger.DisableLog("rdApp.*")

        smiles_mask = self.drug_list["SMILES"].notna()
        smiles = self.drug_list.loc[smiles_mask, "SMILES"]
        molecules = smiles.apply(Chem.MolFromSmiles)
        molecules_mask = molecules.notna()
        fingerprints = self.get_fingerprints(molecules[molecules_mask])
        combined_mask = pd.Series(False, index=self.drug_list.index)
        combined_mask.loc[smiles[molecules_mask].index] = True
        for name, model in self.admet_models.items():
            predictions = model.predict_proba(fingerprints)
            self.drug_list.loc[combined_mask, name] = predictions[slice(None, None), 1]

    def get_fingerprints(self, molecules):
        fingerprints = list()
        fingerprints.append(maplight_gnn.get_morgan_fingerprints(molecules))
        fingerprints.append(maplight_gnn.get_avalon_fingerprints(molecules))
        fingerprints.append(maplight_gnn.get_erg_fingerprints(molecules))
        fingerprints.append(maplight_gnn.get_rdkit_features(molecules))
        fingerprints.append(maplight_gnn.get_gin_supervised_masking(molecules))
        return np.concatenate(fingerprints, axis=1)


if __name__ == "__main__":
    augmenter = (
        DataAugmenter("data/translator_drugs.json")
        .load_drug_queries()
        .load_admet_models(
            {
                "Blood Brain Barrier": "data/admet/bbb_martins-0.916-0.002.dump",
                "Bioavailability": "data/admet/bioavailability_ma-0.74-0.01.dump",
                "Human Intestinal Absorption": "data/admet/hia_hou-0.989-0.001.dump",
            }
        )
    )
    augmenter.drug_list["id_type"] = "cas-number"
    _hy_gensym_f_1 = augmenter
    _hy_gensym_f_1.match_drugbank("data/src/drugbank.xml", "result_id", "id_type", "result_name")
    _hy_gensym_f_1.deduplicate()
    _hy_gensym_f_1.predict_admet()
    _hy_gensym_f_1.save_drug_info("data/translator_drug_list.json")
