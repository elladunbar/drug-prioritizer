#!/usr/bin/env pypy3

import xml.etree.ElementTree as ET

TERMINAL_TYPES = {
    "xs:anyURI": str,
    "xs:boolean": bool,
    "xs:date": str,
    "xs:float": float,
    "xs:integer": int,
    "xs:string": str,
}


class SchemaType:
    def __init__(self, type_name: str, singular: bool):
        self.singular = singular
        self.terminal = type_name in TERMINAL_TYPES


def parse_tree():
    for _, element in ET.iterparse("data/drugbank.xsd", ["end"]):
        print(element.tag)


parse_tree()
