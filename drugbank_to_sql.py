#!/usr/bin/env uv run --script --python pypy
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "xmlschema",
# ]
# ///

import xmlschema

schema = xmlschema.XMLSchema("data/drugbank.xsd")
root = schema.elements["drugbank"]


def walk_element(name, elem_type, depth=0):
    indent = "  " * depth
    if elem_type.is_complex():
        print(f"{indent}{name} (complex)")
        content = elem_type.content
        if hasattr(content, "particles"):
            for particle in content.iter_elements():
                child_name = particle.name
                if child_name:
                    walk_element(child_name, particle.type, depth + 1)
    else:
        print(f"{indent}{name} (simple)")


# Start recursive walk
walk_element(root.name, root.type)
