#!/usr/bin/env pypy3

import xml.etree.ElementTree as ET
from collections import deque

TERMINAL_TYPES = {
    "{http://www.w3.org/2001/XMLSchema}anyURI": str,
    "{http://www.w3.org/2001/XMLSchema}boolean": bool,
    "{http://www.w3.org/2001/XMLSchema}date": str,
    "{http://www.w3.org/2001/XMLSchema}float": float,
    "{http://www.w3.org/2001/XMLSchema}integer": int,
    "{http://www.w3.org/2001/XMLSchema}string": str,
}
TYPE_TAGS = ("element", "simpleType", "complexType")


class SchemaType:
    def __init__(self, type_name: str):
        self.name = type_name
        self.singular = False
        self.list_like = False
        self.terminal = False

        if "list-type" in self.name:
            self.list_like = True
        elif type_name in TERMINAL_TYPES:
            self.singular = True
            self.terminal = True


schema_namespace = "{http://www.w3.org/2001/XMLSchema}"
schema = ET.parse("data/drugbank.xsd")

root_type_element = schema.find(f"./{schema_namespace}element[@name='drugbank']")
type_mapping = {root_type_element.attrib["name"]: root_type_element.attrib["type"]}
unfinished_types = deque([root_type_element.attrib["type"]])
while unfinished_types:
    current_type = unfinished_types.popleft()
    element = schema.find(f".//*[@name='{current_type}']")
    if element.tag == schema_namespace + "element":
        type_mapping[current_type] = element.attrib["type"]
    elif element.tag == schema_namespace + "simpleType":
        restriction = element.find(f"./{schema_namespace}restriction")
        type_mapping[current_type] = restriction.attrib["base"]
    elif element.tag == schema_namespace + "complexType":
        continue
    else:
        continue
