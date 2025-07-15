#!/usr/bin/env uv run

import json
import re
import sys
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib.backends.backend_pdf import PdfPages

##############################
# GET PARTS OF FILE NAME
##############################
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Convert string to datetime object
datetime_object = datetime.strptime(timestamp, "%Y%m%d_%H%M%S")

year = datetime_object.year
month = datetime_object.month
day = datetime_object.day
hour = datetime_object.hour
minute = datetime_object.minute
second = datetime_object.second

print(f"Year: {year}, Month: {month}, Day: {day}, Hour: {hour}, Minute: {minute}, Second: {second}")

##############################
# SET GLOBAL VARIABLES
##############################
# Initialize the list to hold all rows
all_rows = []
edges = {}
nodes = {}
results = []
qualifiers = []
auxiliary_graphs = {}


##############################
# PROCESS SOP
##############################
def get_edge(row_result_data, edge_id):
    # edge ID = support0-UBERON:0002107-has_part-PUBCHEM.COMPOUND:1102-via_subclass
    # row_result_data = {'pk': 'cd108944-3281-4aea-adf5-90200503bbdb', 'ara': 'infores:aragorn', 'result_subjectNode_name': 'carboacyl group', 'result_subjectNode_id': 'CHEBI:37838', 'result_subjectNode_cat': 'biolink:PhysicalEssence', 'result_objectNode_name': 'Bethlem myopathy', 'result_objectNode_id': 'MONDO:0008029', 'result_objectNode_cat': 'biolink:Disease', 'rank': 720, 'sugeno_score': 0.0, 'comp_confidence_score': 0.0, 'comp_novelty_score': 0, 'comp_clinical_evidence_score': 0.0, 'weighted_mean_score': 0.0, 'normalized_score': 49.227, 'ARAX': False, 'ARAX_score': -0.0, 'unsecret': False, 'unsecret_score': -0.0, 'improving_agent': False, 'improving_agent_score': -0.0, 'biothings_explorer': False, 'biothings_explorer_score': -0.0, 'aragorn': True, 'aragorn_score': 0.0, 'ARA_list': ['infores:aragorn'], 'ARA_count': 1, 'result_counter': 3}
    # support_edge = ca63d5667cf0
    # edge_data = row_result_data
    global edges, auxiliary_graphs, qualifiers

    try:
        edge = edges[edge_id]
    except KeyError:
        print(f"KeyError: The key {edge_id} was not found in the 'edges' dictionary.")
        print(f"KeyError: The row_result_data {row_result_data} was not fedge we were on when the error occured.")
        return
        # You can add additional logic here if needed
        # For example, you might want to set a default value or raise a different exception

    edge_objectNode_id = edge["object"]
    edge_subjectNode_id = edge["subject"]

    # CHECK TO SEE IF THE NODE EXISTS IN NODES AND THEN GET THE DATA
    if edge_subjectNode_id in nodes:
        edge_subjectNode_name = nodes[edge_subjectNode_id].get("name", "not provided")
        edge_subjectNode_cats = nodes[edge_subjectNode_id].get("categories", "not provided")
    else:
        edge_subjectNode_name = "not provided"
        edge_subjectNode_cats = "not provided"
        print(f"Node id {edge_subjectNode_id} not found in nodes")

    if edge_objectNode_id in nodes:
        edge_objectNode_name = nodes[edge_objectNode_id].get("name", "not provided")
        edge_objectNode_cats = nodes[edge_objectNode_id].get("categories", "not provided")
    else:
        edge_objectNode_name = "not provided"
        edge_objectNode_cats = "not provided"
        print(f"Node id {edge_objectNode_name} not found in nodes")

    # edge_objectNode_name = nodes[edge_objectNode_id].get('name', 'not provided')
    # edge_objectNode_cats = nodes[edge_objectNode_id].get('categories', 'not provided')
    # edge_subjectNode_name = nodes[edge_subjectNode_id].get('name', 'not provided')
    # edge_subjectNode_cats = nodes[edge_subjectNode_id].get('categories', 'not provided')
    qualifiers = edge.get("qualifiers", [])

    # INITIALIZE QUALIFIED PREDICATE DATA
    qualified_predicate = ""
    causal_mechanism_qualifier = ""
    object_direction_qualifier = ""
    subject_direction_qualifier = ""
    subject_form_or_variant_qualifier = ""
    object_form_or_variant_qualifier = ""
    object_aspect_qualifier = ""
    subject_aspect_qualifier = ""

    if qualifiers is not None and qualifiers != []:
        for qualifier in qualifiers:
            if qualifier["qualifier_type_id"].split(":")[-1] == "qualified_predicate":
                qualified_predicate = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "causal_mechanism_qualifier":
                causal_mechanism_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "object_direction_qualifier":
                object_direction_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "subject_direction_qualifier":
                subject_direction_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "subject_form_or_variant_qualifier":
                subject_form_or_variant_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "object_form_or_variant_qualifier":
                object_form_or_variant_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "object_aspect_qualifier":
                object_aspect_qualifier = qualifier["qualifier_value"].split(":")[-1]

            elif qualifier["qualifier_type_id"].split(":")[-1] == "subject_aspect_qualifier":
                subject_aspect_qualifier = qualifier["qualifier_value"].split(":")[-1]

            # elif qualifier['qualifier_type_id'].split(':')[-1] == 'causal_mechanism_qualifier':
            #   causal_mechanism_qualifier = qualifier['qualifier_value'].split(':')[-1]

            # elif qualifier['qualifier_type_id'].split(':')[-1] == 'causal_mechanism_qualifier':
            #   causal_mechanism_qualifier = qualifier['qualifier_value'].split(':')[-1]

    edge_data = {
        # 'objectNode_name' : objectNode_name,
        # 'subjectNode_name' : subjectNode_name,
        "edge_id": edge_id,
        "edge_object": edge["object"],
        "edge_objectNode_name": edge_objectNode_name,
        "edge_objectNode_cats": edge_objectNode_cats,
        "edge_objectNode_cat": edge_objectNode_cats[0],
        "edgeg_subject": edge["subject"],
        "edge_subjectNode_name": edge_subjectNode_name,
        "edge_subjectNode_cats": edge_subjectNode_cats,
        "edge_subjectNode_cat": edge_subjectNode_cats[0],
        "predicate": edge["predicate"],
        "qualifiers": qualifiers,
        "edge_type": "one-hop",
        "qualified_predicate": qualified_predicate,
        "causal_mechanism_qualifier": causal_mechanism_qualifier,
        "subject_direction_qualifier": subject_direction_qualifier,
        "subject_aspect_qualifier": subject_aspect_qualifier,
        "subject_form_or_variant_qualifier": subject_form_or_variant_qualifier,
        "object_direction_qualifier": object_direction_qualifier,
        "object_aspect_qualifier": object_aspect_qualifier,
        "object_form_or_variant_qualifier": object_form_or_variant_qualifier,
        "publications": [],
    }

    #  GET PHRASE
    edge_data["phrase"] = recombobulation(edge_data)

    # GET SOURCE INFORMATION
    agg_counter = 1
    sources = edge["sources"]

    # CHECK FOR PRIMARY KS AND AGGREAGATOR
    for source in sources:
        role = source["resource_role"]
        resource_id = source["resource_id"]

        if role == "primary_knowledge_source":
            edge_data["primary_source"] = resource_id
        elif role == "aggregator_knowledge_source" and agg_counter <= 2:
            edge_data[f"agg{agg_counter}"] = resource_id
            agg_counter += 1

    # CHECK FOR SUPPORT GRAPHS AND SET
    try:
        attributes = edge.get("attributes", [])
    except Exception as e:
        print(f"edge_id = {edge_id}")  # json.dumps(edge, indent=4)
        print(f"edge = {json.dumps(edge, indent=4)}")  # json.dumps(edge, indent=4)
        print(f"An unexpected error occurred: {e}")
        exit()
    has_supportgraphs = False
    support_graphs_ids = []

    if attributes is not None and attributes != []:
        try:
            for attribute in attributes:
                if attribute["attribute_type_id"] == "biolink:support_graphs":
                    edge_data["edge_type"] = "creative"
                    has_supportgraphs = True
                    support_graphs_ids = attribute["value"]
                if attribute["attribute_type_id"] == "biolink:publications":
                    edge_data["publications"] = attribute["value"]

        except Exception as e:  # This will catch any other type of exception
            print(f"edge_id = {edge_id}")
            print(f"attributes = {attributes}")
            print(f"An unexpected error occurred: {e}")
            exit()

    edge_data["publications_count"] = len(edge_data["publications"])

    # ADD THE NEW EDGE DATA TO THE RESULT DATA AND APPEND TO ALL ROWS
    result_edge_data = {**row_result_data, **edge_data}

    all_rows.append(result_edge_data)

    ##############################
    # LOOP AGAIN IF THERE ARE SUPPORT GRAPHS - auxiliary_graphs
    ##############################

    # GET EDGES FROM SUPPORT GRAPH
    if has_supportgraphs == True:
        for support_graph_id in support_graphs_ids:
            try:
                aux_graph = auxiliary_graphs[support_graph_id]
                support_edges = aux_graph["edges"]
                # print(f"support_edges = {support_edges}")

                for support_edge in support_edges:
                    # print(f"support_edge = {support_edge}")
                    get_edge(row_result_data, support_edge)
            except Exception as e:
                print("aux_graph error")
                print(f"edge ID = {support_graph_id}")
                print(f"row_result_data = {row_result_data}")
                print(f"support_edge = {support_edge}")
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)


##############################
# REASSEMBLE THE SPO STATEMENT
##############################
def recombobulation(edge_data):
    objectNode_name = edge_data["edge_objectNode_name"]
    subjectNode_name = edge_data["edge_subjectNode_name"]
    # predicate = edge_data['predicate']
    direction = ""
    direction = ""
    object_aspect = ""
    subject_aspect = ""
    object_of = ""

    ##############################
    # GET PARTS OF THE PHRASE
    # DOWNREGULATES AND ABUNDANCE HAVE CHANGES THAT SEEM TO BE CORRECTING FOR A MISTAKE
    ##############################

    # SET VARIABLE TO HOLD THE PREDICATE TO BE USED IN THE PHRASE
    predicate_used = edge_data["predicate"].replace("biolink:", "", 1)

    # IF THERE IS A object_direction_qualifier THEN THERE SHOULD BE A
    # CREATE VARIABLE TO SHOW IF I SHOULD MODIFY THE QUALIFIED PREDICATE
    chang_qualified_predicate = False
    # GET COMPONENTS FOR PHRASE
    # ## ORDER MATTER FOR THE PROCESS BELOW

    #  GET REPLACE THE PREDICATE WITH A QUALIFIED IF EXISTS
    if edge_data["qualified_predicate"] != "":
        predicate_used = edge_data["qualified_predicate"].replace("biolink:", "", 1)

    # GET THE ASPECT THAT IS REFERED TO FOR THE OBJECT
    if edge_data["object_aspect_qualifier"] != "":
        object_aspect = edge_data["object_aspect_qualifier"]
        if object_aspect == "abundance":
            object_aspect = "the abundance"

    # GET THE ASPECT THAT IS REFERED TO FOR THE SUBJECT
    if edge_data["subject_aspect_qualifier"] != "":
        subject_aspect = edge_data["subject_aspect_qualifier"]
        if subject_aspect == "abundance":
            subject_aspect = "the abundance"

    # IF THERE IS AN object_direction_qualifier THEN THE PREDICATE READS BETTER IF THERE IS A QUALIFIED PREDICATE CAUSES
    if edge_data["object_direction_qualifier"] != "":
        direction = edge_data["object_direction_qualifier"]
        predicate_used = "causes"
        if edge_data["object_direction_qualifier"] == "downregulated":
            predicate_used = "downregulated"

    # ADD THE 'OF' TO MAKE THE PHRASE BETTER
    if edge_data["object_aspect_qualifier"] != "":
        object_of = "of"

    # CHANGE THE TENSE OF qualified_predicate
    if edge_data["qualified_predicate"] == "" and edge_data["qualified_predicate"] != "":
        direction = direction[:-1] + "s"

    # DEFAULT PHRASE
    infered_phrase = "DEFAULT PHRASE"
    # PUT THE PHRASE TOGETHER
    if edge_data["qualified_predicate"] == "causes":
        infered_phrase = (
            f"{subjectNode_name} {predicate_used} {direction} {object_aspect} {object_of} {objectNode_name}"
        )

    if edge_data["qualified_predicate"] == "caused_by":
        # infered_phrase = (f"{subjectNode_name} {predicate_used} {direction} {object_aspect} {object_of} {objectNode_name}")
        infered_phrase = f"{edge_data['subject_direction_qualifier']} {edge_data['subject_aspect_qualifier']} of {objectNode_name} is {edge_data['qualified_predicate']} {subjectNode_name}"
    # subject_direction_qualifier
    # THE PROCESS ABOVE CAN RESULT IN DOUBLE SPACES - THIS REMOVES THEM
    infered_phrase = re.sub(" +", " ", infered_phrase)
    infered_phrase = re.sub("_", " ", infered_phrase)

    return infered_phrase


##############################
# Function to run when button is clicked
##############################
def run_on_click(b):
    global file_count
    global edges, nodes, resuls, auxiliary_graphs
    # result_count = 0
    # try:
    print("STARTED")
    # Get the base URL for the selected environment
    base_url = url_dict.get(env_input.value)
    pk_for_file = pk_input.value

    # Fetch data from the ARS API
    print("link with trace")
    print(f"{base_url}/ars/api/messages/{pk_input.value}?trace=y")
    response = requests.get(f"{base_url}/ars/api/messages/{pk_input.value}?trace=y")
    data = response.json()

    print("GETTING MERGED PK")
    merged_version = data["merged_version"]
    print(f"{base_url}/ars/api/messages/{merged_version}")

    # GET CHILD DATA
    mergeResponse = requests.get(f"{base_url}/ars/api/messages/{merged_version}")
    merged_json_data = mergeResponse.json()
    print("API response returned")
    # print(json.dumps(data, indent=2))  # Pretty print the API response

    # SAVE EARLY FOR TROUBLE SHOOTING
    if attribute_check.value == "Do Download JSON":
        save_json(merged_json_data, "nonCreativeMerged")

    # Initialize the result counter
    result_counter = 0

    # GET ALL OF THE PARTS OF THE RESULT NEEDED BROKEN OUT FOR USABILITY
    results = merged_json_data["fields"]["data"]["message"]["results"]
    print(f"results length = {len(results)}")
    nodes = merged_json_data["fields"]["data"]["message"]["knowledge_graph"]["nodes"]
    print(f"nodes length = {len(nodes)}")

    edges = merged_json_data["fields"]["data"]["message"]["knowledge_graph"]["edges"]
    auxiliary_graphs = merged_json_data["fields"]["data"]["message"]["auxiliary_graphs"]

    # LOOP THROUGH RESULTS
    print("looping through results")
    for result in results:
        # RESULT COUNT
        result_counter += 1

        # GET RESULT DATA
        comp_novelty = "N/A"
        comp_confidence = "N/A"
        comp_clinical_evidence = "N/A"

        rank = result.get("rank", "N/A")
        sugeno = result.get("sugeno", "N/A")
        weighted_mean = result.get("weighted_mean", "N/A")
        normalized_score = result.get("normalized_score", "N/A")
        ordering_components = result.get("ordering_components", [])
        if ordering_components != []:
            comp_novelty = result["ordering_components"]["novelty"]
            comp_confidence = result["ordering_components"]["confidence"]
            comp_clinical_evidence = result["ordering_components"]["clinical_evidence"]

        # rank = result['rank']
        # sugeno = result['sugeno']
        # weighted_mean = result['weighted_mean']
        # normalized_score = result['normalized_score']
        # comp_novelty = result['ordering_components']['novelty']
        # comp_confidence = result['ordering_components']['confidence']
        # comp_clinical_evidence = result['ordering_components']['clinical_evidence']

        #  GET RESULT NODES INFORMATION
        # objectNode_id = result['node_bindings']['on'][0]['id']
        # subjectNode_id = result['node_bindings']['sn'][0]['id']

        # GET KEYS - PREVIOUSLY ON AND SN
        node_bindings = result["node_bindings"]
        node_binding_keys = node_bindings.keys()
        node_group_one = []
        node_group_two = []
        node_group_one_names = []
        node_group_two_names = []
        node_group_one_cat = []
        node_group_two_cat = []
        node_group_counter = 1

        for key in node_binding_keys:
            node_group_array = node_bindings[key]
            try:
                for node_group in node_group_array:
                    if node_group_counter == 1:
                        node_group_one = node_group["id"]
                        # node_id = node_group['id']
                        node_group_one_names = nodes[node_group_one].get("name", "N/A")
                        node_group_one_cat = nodes[node_group_one].get("categories", ["N/A"])
                    if node_group_counter == 2:
                        node_group_two = node_group["id"]
                        # node_id = node_group['id']
                        node_group_two_names = nodes[node_group_two].get("name", "N/A")
                        node_group_two_cat = nodes[node_group_two].get("categories", ["N/A"])
                node_group_counter += 1
            except Exception as e:
                print(f"error = {e}")
                print(f"node_group = {node_group}")

        try:
            result_subjectNode_name = node_group_two_names
            result_subjectNode_id = node_group_two
            result_subjectNode_cat = node_group_two_cat[0]
            result_objectNode_name = node_group_one_names
            result_objectNode_id = node_group_one
            result_objectNode_cat = node_group_one_cat[0]
        except Exception as e:
            print(f"Exception = {e}")
            print(f"nodes = {nodes}")
            # print(f"KeyError: The key {edge_id} was not found in the 'edges' dictionary.")
            # print(f"KeyError: The row_result_data {row_result_data} was not fedge we were on when the error occured.")
            print(f"node_group_one_cat: {node_group_one_cat}")
            print(f"node_group_two_cat: {node_group_two_cat}")
            # return

        # GET THE ARAS THAT CONTRIBUTED TO
        ## SET VALUES TO DEFAULT STATES
        improving_agent = False
        improving_agent_score = -0.0001

        ARAX = False
        ARAX_score = -0.0001

        unsecret = False
        unsecret_score = -0.0001

        improving_agent = False
        improving_agent_score = -0.0001

        biothings_explorer = False
        biothings_explorer_score = -0.0001

        aragorn = False
        aragorn_score = -0.0001

        ## CHECK FOR DATA
        analyses = result["analyses"]
        # print(analyses)
        # ara = analyses['resource_id']
        aras = []

        for analysis in analyses:
            # print("analysis['resource_id']")
            # print(analysis['resource_id'])
            ara = analysis["resource_id"]
            if analysis["resource_id"] == "infores:improving-agent":
                aras.append("infores:improving-agent")
                improving_agent = True
                improving_agent_score = analysis["score"]

            elif analysis["resource_id"] == "infores:rtx-kg2":
                aras.append("infores:ARAX_rtx-kg2")
                ARAX = True
                ARAX_score = analysis["score"]

            elif analysis["resource_id"] == "infores:improving-agent":
                aras.append("infores:improving-agent")
                improving_agent = True
                improving_agent_score = analysis["score"]

            elif analysis["resource_id"] == "infores:biothings-explorer":
                aras.append("infores:biothings-explorer")
                biothings_explorer = True
                biothings_explorer_score = analysis["score"]

            elif analysis["resource_id"] == "infores:unsecret-agent":
                aras.append("infores:unsecret-agent")
                biotunsecrethings_explorer = True
                unsecret_score = analysis["score"]

            elif analysis["resource_id"] == "infores:aragorn":
                aras.append("infores:aragorn")
                aragorn = True
                aragorn_score = analysis["score"]

            else:
                aras.append(analysis["resource_id"])

        # CREATE ROW
        row_result_data = {
            "pk": pk_input.value,
            "ara": ara,
            "result_subjectNode_name": result_subjectNode_name,
            "result_subjectNode_id": result_subjectNode_id,
            "result_subjectNode_cat": result_subjectNode_cat,
            "result_objectNode_name": result_objectNode_name,
            "result_objectNode_id": result_objectNode_id,
            "result_objectNode_cat": result_objectNode_cat,
            "rank": rank,
            "sugeno_score": sugeno,
            "comp_confidence_score": comp_confidence,
            "comp_novelty_score": comp_novelty,
            "comp_clinical_evidence_score": comp_clinical_evidence,
            "weighted_mean_score": weighted_mean,
            "normalized_score": normalized_score,
            "ARAX": ARAX,
            "ARAX_score": ARAX_score,
            "unsecret": unsecret,
            "unsecret_score": unsecret_score,
            "improving_agent": improving_agent,
            "improving_agent_score": improving_agent_score,
            "biothings_explorer": biothings_explorer,
            "biothings_explorer_score": biothings_explorer_score,
            "aragorn": aragorn,
            "aragorn_score": aragorn_score,
            "ARA_list": aras,
            "ARA_count": len(aras),
            "result_counter": result_counter,
        }

        for key in row_result_data.keys():
            if key.endswith("_score"):
                # row_result_data[key] = round(row_result_data[key], 3)
                # ALLOW FOR N/A TO BE THERE BUT ROUND IF NUMBER
                if isinstance(row_result_data[key], (int, float)):
                    row_result_data[key] = round(row_result_data[key], 3)

        # GO BACK THROUGH ALL OF THE ANALYSES GET THE T_EDGES
        for analysis in analyses:
            edge_bindings = analysis["edge_bindings"]
            edge_bindings_keys = edge_bindings.keys()
            # CHECK FOR EDGE BINDINGS  - LIKLEY ONLY 1
            for edge_bindings_key in edge_bindings_keys:
                edge_objects = edge_bindings[edge_bindings_key]
                # GET EDGE IDS
                for edge_object in edge_objects:
                    edge_id = edge_object["id"]
                    get_edge(row_result_data, edge_id)

    print("PRINTING ALL ROWS")
    print(f"# Results = {len(results)}")
    print(f"# Rows found = {len(all_rows)}")
    # print(all_rows)
    save_file(all_rows, merged_json_data)
    # make_histo(all_rows)


##################
# CREATE HISTOGRAM OF UNIQUE SCORES
# THIS IS RUN FROM THE SAVE CSV FUNCTION SO THAT THE NAME CAN BE REUSED
##################
# Function to safely divide two numbers
def safe_divide(numerator, denominator):
    if denominator == 0:
        return 0
    else:
        return numerator / denominator


def make_histo(all_rows, calc_pdf_filename):
    # Convert all_rows to a pandas DataFrame
    df = pd.DataFrame(all_rows)

    # Deduplicate the DataFrame based on 'result_counter'
    df_unique = df.drop_duplicates(subset="result_counter").copy()

    # Normalize 'normalized_score' by dividing by 100
    df_unique["normalized_score"] = df_unique["normalized_score"].apply(
        lambda x: x / 100 if isinstance(x, (int, float)) else x
    )

    # Name of the PDF file
    pdf_filename = "data/translator/" + calc_pdf_filename + ".pdf"

    # Create a PdfPages object
    with PdfPages(pdf_filename) as pdf:
        # Define the order of columns for each page
        page1_columns = [
            "aragorn_score",
            "ARAX_score",
            "biothings_explorer_score",
            "improving_agent_score",
            "unsecret_score",
        ]
        page2_columns = ["sugeno_score", "comp_confidence_score", "weighted_mean_score", "normalized_score"]

        # CREATE PAGE 1
        fig, axes = plt.subplots(3, 2, figsize=(15, 10))  # 3x2 grid for 5 plots
        axes = axes.flatten()
        for i, column in enumerate(page1_columns):
            ax = axes[i]
            # CREATE PAGE 1
            # ... [code to create histogram for each 'column' in page1_columns]
            if column in df_unique.columns:
                # Drop rows with N/A scores and filter out zero values
                hist_data = df_unique[df_unique[column] != "N/A"]
                hist_data = hist_data[hist_data[column] > 0]

                # Convert column to numeric type
                hist_data[column] = pd.to_numeric(hist_data[column])

                # Calculate the percentage for each bin
                bin_edges = np.arange(0.05, 1.05, 0.05)  # Bins from 0.05 to 1.0 in increments of 0.05
                hist, _ = np.histogram(hist_data[column], bins=bin_edges)
                total = hist.sum()
                hist_percentage = (hist / total) * 100 if total > 0 else np.zeros_like(hist)

                # Plot histogram on the selected axis
                patches = ax.bar(bin_edges[:-1], hist_percentage, width=0.05, align="edge", alpha=0.7)

                # Add value above each bar
                for patch, value in zip(patches, hist_percentage):
                    if value > 0:
                        ax.text(patch.get_x() + patch.get_width() / 2, value, f"{value:.0f}%", ha="center", va="bottom")

                # Set x-axis and y-axis limits
                ax.set_xlim(0.05, 1.0)
                ax.set_ylim(0, 100)

                # Add titles and labels
                ax.set_title(f"Histogram of {column}")
                ax.set_xlabel("Score")
                ax.set_ylabel("Percentage of Non-Zero Scores")

                # Calculate and display the percentage of unique scores
                unique_scores = len(hist_data[column].unique())
                non_zero_scores = len(hist_data)
                percent_unique = safe_divide(unique_scores, non_zero_scores) * 100
                ax.title.set_size(10)
                ax.set_title(f"{ax.get_title()}, % Unique: {percent_unique:.2f}", fontsize=10)
            else:
                ax.axis("off")  # Turn off axis if the column is not in the dataframe

        plt.suptitle(calc_pdf_filename + " page 2")  # Set page title
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout
        pdf.savefig(fig)  # Save the first page
        plt.close(fig)  # Close the figure after saving

        # Page 2
        # ... [earlier parts of your code]

        # Page 2
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))  # 2x2 grid for 4 plots
        axes = axes.flatten()
        for i, column in enumerate(page2_columns):
            ax = axes[i]
            if column in df_unique.columns:
                # Drop rows with N/A scores and filter out zero values
                hist_data = df_unique[df_unique[column] != "N/A"]
                hist_data = hist_data[hist_data[column] > 0]

                # Convert column to numeric type
                hist_data[column] = pd.to_numeric(hist_data[column])

                # Calculate the percentage for each bin
                bin_edges = np.arange(0.05, 1.05, 0.05)  # Bins from 0.05 to 1.0 in increments of 0.05
                hist, _ = np.histogram(hist_data[column], bins=bin_edges)
                total = hist.sum()
                hist_percentage = (hist / total) * 100 if total > 0 else np.zeros_like(hist)

                # Plot histogram on the selected axis
                patches = ax.bar(bin_edges[:-1], hist_percentage, width=0.05, align="edge", alpha=0.7)

                # Add value above each bar
                for patch, value in zip(patches, hist_percentage):
                    if value > 0:
                        ax.text(patch.get_x() + patch.get_width() / 2, value, f"{value:.0f}%", ha="center", va="bottom")

                # Set x-axis and y-axis limits
                ax.set_xlim(0.05, 1.0)
                ax.set_ylim(0, 100)

                # Add titles and labels
                ax.set_title(f"Histogram of {column}")
                ax.set_xlabel("Score")
                ax.set_ylabel("Percentage of Non-Zero Scores")

                # Calculate and display the percentage of unique scores
                unique_scores = len(hist_data[column].unique())
                non_zero_scores = len(hist_data)
                percent_unique = safe_divide(unique_scores, non_zero_scores) * 100
                ax.title.set_size(10)
                ax.set_title(f"{ax.get_title()}, % Unique: {percent_unique:.2f}", fontsize=10)
            else:
                ax.axis("off")  # Turn off axis if the column is not in the dataframe

        # Adjust layout and save the figure
        # plt.tight_layout()
        # pdf.savefig(fig)  # Save the second page
        # plt.close(fig)  # Close the figure after saving

        # ... [code for downloading the PDF]

        plt.suptitle(calc_pdf_filename + " page 2")  # Set page title
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout
        pdf.savefig(fig)  # Save the second page
        plt.close(fig)  # Close the figure after saving


##################
# SAVE DATA AS CSV
# AND RUN THE PDF OF HISTOGRAM FUNCTION
##################
def save_file(all_rows, merged_json_data):
    global qualifiers
    # Download the CSV file
    print("ABOUT TO SAVE CSV")
    pk_for_file = pk_input.value

    # MAKE PANDA DATAFRAME
    df = pd.DataFrame(all_rows)
    print("MAKING DF")
    df.head()

    # GET GENE LIST FROM objectNode_name
    # gene_name_list = df.explode('objectNode_name')
    # gene_id_list = df.explode('objectNode_id')
    # qualifiers_list = df.explode('qualifiers')
    # print(f"gene_list = {gene_list}")
    # unique_gene_names = gene_name_list['objectNode_name'].unique()
    # unique_gene_id = gene_id_list['objectNode_id'].unique()
    # unique_qualifiers_list = qualifiers_list['qualifiers']
    # print(f"unique_genes = {unique_gene_names}")
    # print(f"unique_gene_id = {unique_gene_id}")
    # print(f"qualifiers = {qualifiers}")

    # GET NAME OF CREATIVE TARGET
    most_common_object_name = ""
    most_common_object_name = df["result_objectNode_name"].mode()[0]

    # CREATE NAME FROM INFO
    filename = "{}_{}_{}_{}_{}_{}_{}".format(most_common_object_name, pk_for_file, year, month, day, hour, minute)
    print(filename)  # Outputs: 'Year_Month_Day_Hour_Minute_Second'

    # RUN THE HISTOGRAM FUNCTION
    make_histo(all_rows, filename)
    # Save the DataFrame to a CSV file with your filename
    df.to_csv(f"data/translator/{filename}.csv", index=False)

    del df

    # check to see if you want to download the json also
    print(f"save json value = {attribute_check.value}")
    if attribute_check.value == "Do Download JSON":
        save_json(merged_json_data, filename)


##################
# SAVE DATA AS JSON
##################
def save_json(data, name):
    print("ABOUT TO SAVE JSON")
    # Convert the Python dictionary to a JSON string
    json_string = json.dumps(data, indent=4)  # indent=4 for pretty formatting

    # Save the JSON string to a file
    with open(f"data/translator/{name}.json", "w") as json_file:
        json_file.write(json_string)


# Dictionary mapping environment names to base URLs
url_dict = {
    "test": "https://ars.test.transltr.io",
    "CI": "https://ars.ci.transltr.io",
    "dev": "https://ars-dev.transltr.io",
    "prod": "https://ars-prod.transltr.io",
}


# Duck typing
class ValueWrapper:
    def __init__(self, arg):
        self.value = arg


if __name__ == "__main__":
    pk_input = ValueWrapper("297cc1a8-79c7-40a7-8546-05f5757836b8")
    env_input = ValueWrapper("prod")
    attribute_check = ValueWrapper("DO NOT Download JSON")

    run_on_click(None)
