# Organizing Thoughts

## DrugBank XML to SQL Transition

### Currently

- drugbank.xml is hierarchical text with schema stored in drugbank.xsd
- want to normalize and turn into tables with foreign keys in drugs.db

### To Transition

- flatten drugbank.xml into something that resembles main table and child tables
- make sure types match up according to schema
    - dictionary per drug?

### Algorithms

Create SQL tables/schema from XML schema


```
1 parse whole document at once
2 for element in document
    3 if no attribute "name"; continue
    4 if type is in terminal types, add to mapping as-is
    5 
```

Create SQL from XML

```
1 create tables according to parsed schema
2 parse document -> element
    3 if element is "drug", then create new dictionary -> drug
    4 for event, tag in drug:
        5 if tag is singular: drug[tag] = text
        6 else if list-like: drug[tag] = list(sub_text)
        7 else: create new dictionary to push into
        8 pop up dictionary on event-end
    9 send dictionary to database using schema
```
