import json

import pandas
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pandas import core
import numpy as np
import openpyxl
import random
import requests
import time
from operator import getitem
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.api as sm
from scipy import stats
from scipy.stats import f_oneway, kstest, mannwhitneyu, false_discovery_control

def main():

    p_list = []

    name_list = ["BRAF", "CFTR", "DDX3Y", "PDE5A", "SLC5A2", "VKORC1",
                 "HBB", "ACE", " BRCA2", "ADRB2", "CYP2D6", "SMARCE1"]

    p_list.append(take_query('', name_list[3]))
    p_list.append(take_query('', name_list[3]))
    p_list.append(take_query('', name_list[3]))
    p_list.append(take_query('', name_list[3]))
    p_list.append(take_query('02489e70-d19e-42dd-9cf7-822fa5dfdede', name_list[4]))
    p_list.append(take_query('caed9eff-da30-42cd-b351-669210461384', name_list[5]))

    p_list.append(take_query('ebe022f4-652c-4380-b9f6-da7a3eda0887', name_list[6]))
    p_list.append(take_query('ce531d16-f25c-4a22-a047-5909c915acc0', name_list[7]))
    p_list.append(take_query('fc98b929-346d-479e-86c8-f056f5cca3af', name_list[8]))
    p_list.append(take_query('854676de-c32f-4475-95fa-eac71c6e9454', name_list[9]))
    p_list.append(take_query('58028aa9-1ba9-4022-a55d-86b7811fce3d', name_list[10]))
    p_list.append(take_query('b47a79bd-83b5-4688-a142-e3b9642450d4', name_list[11]))

    rejected, p_vals_corrected, x, y = multipletests(p_list, alpha=0.05, method="fdr_bh")


    p_df = pd.DataFrame({"Genes": name_list, "Values":p_list, "Adjusted P Values": p_vals_corrected, "Outcome" : rejected})
    print()

def take_query(api, filePath):
    url = "https://ars.test.transltr.io/ars/api/messages/" + api + "?trace=y"
    r = requests.get(url, allow_redirects=True)
    time.sleep(360)

    json_data = r.json()
    merged_version_pk = json_data['merged_version']
    url2 = "https://ars.test.transltr.io/ars/api/messages/" + merged_version_pk
    r2 = requests.get(url2, allow_redirects=True)
    open('results.json', 'wb').write(r2.content)
    gene = json.load(open('results.json', 'rb'))

    # Data with only results including atc_classifications
    gene_data = run(gene)
    gene_def_data = run2(gene)
    all_gene_scores = displayScores(gene)
    all_gene_scores.sort()

    # Data with sugeno averages for items with same atc_name, also has count
    gene_stats = run_stats(gene_data)


    #Graph of scores with atc classifications
    plt.hist((gene_data["sugeno"]).sort_values(), bins=50)
    plt.title("Atc Scores")
    plt.xlabel(filePath)
    plt.show()

    #Graph without atc classifications
    plt.hist((gene_def_data["sugeno"]).sort_values(), bins = 50)
    plt.title("Non-Atc Scores")
    plt.xlabel(filePath)
    plt.show()

    #Statistical Test

    atc_scores = np.array(gene_data["sugeno"])
    non_atc_scores = np.array(gene_def_data["sugeno"])
    stat, p_value = mannwhitneyu(atc_scores, non_atc_scores)
    print(filePath)
    print("Statistic: ", stat)
    print("P-Value: ", p_value)
    print()


    # Each atc_code in the set, with its list of respective sugeno scores
    atc_codes = atc_dict(gene_data)

    atc_code_index = {}
    index = 0
    for code in atc_codes:
        atc_code_index[index] = (code, atc_codes[code]["names"], atc_codes[code]["scores"])
        index += 1

    atc_code_df = pd.DataFrame.from_dict(atc_code_index, orient='index', columns=['code', 'names', 'scores'])
    return p_value






def atc_dict(df):
    temp = {}
    for row in df.iterrows():
        sugeno = row[1]["sugeno"]
        atc_code = row[1]["atc_code"]
        atc_name = row[1]["atc_name"]
        if atc_code in temp.keys():
            temp[atc_code]["scores"].append(sugeno)
            temp[atc_code]["names"].append(atc_name)
        else:
            temp[atc_code] = {
                "scores": [sugeno],
                "names": [atc_name]
            }

    return temp
def run_stats(df):
    stats = {}
    temp = {}
    i = 0
    for row in df.iterrows():
        atc_name = row[1]["atc_name"]
        sugeno = row[1]["sugeno"]
        atc_code = row[1]["atc_code"]

        if atc_name in stats.keys():
            stats[atc_name]["count"] = stats[atc_name]["count"] + 1
            stats[atc_name]["sugeno_average"] = (sugeno * (1/stats[atc_name]["count"])) + ((stats[atc_name]["sugeno_average"]) * (stats[atc_name]["count"] - 1) / (stats[atc_name]["count"]))
        else:
            stats[atc_name] = {
                "count": 1,
                "sugeno_average": sugeno,
            }

    df2 = pd.DataFrame.from_dict(stats, orient='index')
    df2.sort_values(by=['count'], ascending=False, inplace=True)
    return df2


def loadPandas(scoreSet):
    i = 1
    pandasData = []

    for item in scoreSet:
        score = (scoreSet[item]['sugeno'])
        name = (scoreSet[item]['node_bindings']['sn'][0]['id'])
        for atc_instance in scoreSet[item]['atc']:
            atcName = (atc_instance['level4']['name'])
            atcCode = (atc_instance['level4']['code'])
            data =  {
                'name' : name,
                'sugeno' : score,
                'atc_name' : atcName,
                'atc_code' : atcCode
            }
            pandasData.append(data)

    df = pd.DataFrame(pandasData)
    return df

def loadPandas2(scoreSet):
    i = 1
    pandasData = []

    for item in scoreSet:
        score = (scoreSet[item]['sugeno'])
        name = (scoreSet[item]['node_bindings']['sn'][0]['id'])
        data = {
            'name' : name,
            'sugeno' : score
        }
        pandasData.append(data)

    df = pd.DataFrame(pandasData)
    return df
def run(data):

    # Checks to see the the data is not empty/null
    results = getResults(data)

    # Loads the knowledgeGraph from the results dictionary
    knowledgeGraph = getKnowledgeGraph(data)

    # Loads the list of result names
    idlist = getresultIDs(results)
    nodeMap = {}

    # Searches the knowledge graph for the result names and then loads the nodes with that name
    for zids in idlist:
        for zid in zids:

            nodeMap[zid] = getNode(knowledgeGraph, zid)

    #Map of results with atc classifications
    atcMap = getMaps(nodeMap)

    #Returns a set of
    pr = processResults(results, atcMap)

    #Retrieves the scores with atc attributes from the processed results
    atcResultsSet = showAtcData(pr)

    #Loads score onto dataframe to return
    dataframe = loadPandas(atcResultsSet)


    return dataframe

def run2(data):

    results = getResults(data)


    knowledgeGraph = getKnowledgeGraph(data)


    idlist = getresultIDs(results)


    nodeMap = {}
    for zids in idlist:
        for zid in zids:

            nodeMap[zid] = getNode(knowledgeGraph, zid)


    atcMap = getMaps(nodeMap)

    pr = processResults(results, atcMap)
    atcResultsSet = showData(pr)
    dataframe = loadPandas2(atcResultsSet)


    return dataframe
def displayScores(data):
    results = getResults(data)
    scoreSet = getAllScores(results)
    return scoreSet

def getAtcByLevel(result, level):
    atcNames = []
    stringLevel = ""
    match level:
        case 1:
            stringLevel = "level1"
        case 2:
            stringLevel = "level2"
        case 3:
            stringLevel = "level3"
        case 4:
            stringLevel = "level4"
        case 5:
            stringLevel = "level5"

    if "atc" in result.keys():
        for atc_instance in result["atc"]:

            atcString = atc_instance[stringLevel]["name"]
            atcNames.append(atcString)
        return atcNames
    else:
        return None

def showData(pr):
    scoreSet = {}
    i = 0
    for att in pr:
        if "atc" not in att:
            scoreSet[i] = pr[i]
        i += 1

    return scoreSet
def showAtcData(pr):
    scoreSet = {}
    i = 0
    for att in pr:
        if "atc" in att:
            scoreSet[i] = pr[i]
        i += 1

    return scoreSet


def processResults(results, atcMap):
    for result in results:
        node_bindings = get_safe(result, "node_bindings")
        # print()
        for key,value in node_bindings.items():
            for item in value:
               # print()
                if "id" in item.keys() and item["id"] is not None:
                    if item["id"] in atcMap.keys():
                        result["atc"] = atcMap[item["id"]]
    return results

def getAllScores(results):
    scores = []
    for item in results:
        scores.append(item["sugeno"])
    return scores
def getMaps(nodeMap):
    atcMap = {}

    for key,value in nodeMap.items():
        attributes = get_safe(value, "attributes")
        if attributes is not None:

            for attribute in attributes:
                if "value" in attribute.keys():

                    if attribute["value"] is not None:
                        attribute_values = attribute["value"]

                        if isinstance(attribute_values, list):
                            for attribute_value in attribute_values:

                                if isinstance(attribute_value, dict):

                                    if "atc_classifications" in attribute_value.keys():

                                        if attribute_value["atc_classifications"] is not None:
                                            atcMap[key] = attribute_value["atc_classifications"]

    return atcMap
def getNode(knowledgeGraph, nodeID):
    nodes = get_safe(knowledgeGraph, "nodes")
    for key,value in nodes.items():
        if nodeID == key:
            return value


def getresultIDs(results):
    id_lists_lists = []
    if results is not None:
        for result in results:
            nb = get_safe(result, "node_bindings")
            id_lists = []

            for key, value in nb.items():
                #print(key)
               # print(value)

                for obj in value:

                    for key, value in obj.items():

                        if key == "id":
                            id_lists.append(value)
                       # print(str(id_lists))
                        id_lists_lists.append(id_lists)
    return id_lists_lists


def getResults(data):
    results = get_safe(data, "fields","data","message","results")
    return results

def getKnowledgeGraph(data):
    knowledgeGraph = get_safe(data, "fields","data","message","knowledge_graph")
    return knowledgeGraph

def get_safe(element,*keys):
    '''
    :param element: JSON to be processed
    :param keys: list of keys in order to be traversed. e.g. "fields","data","message","results
    :return: the value of the terminal key if present or None if not
    '''
    if element is None:
        return None
    _element = element
    for key in keys:
        try:
            _element = _element[key]
            if _element is None:
                return None
            if key == keys[-1]:
                return _element
        except KeyError:
            return None
    return None


if __name__ == '__main__':
    main()








