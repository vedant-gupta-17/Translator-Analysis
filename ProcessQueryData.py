import json

import pandas
import pandas as pd
from pandas import core
import numpy as np
import openpyxl
import random
import requests
import time
from operator import getitem
import matplotlib.pyplot as plt
from scipy.stats import kstest
from scipy.stats import f_oneway
import statsmodels.api as sm


def main():

    processData('/Users/vedant/HBB.json')
    processData('/Users/vedant/ACE.json')
    processData('/Users/vedant/Desktop/Translator NIH/JSON data/BRCA2.json')
    processData('/Users/vedant/ADRB2.json')
    processData('/Users/vedant/CYP2D6.json')
    processData('/Users/Vedant/SMARCE1.json')
def processData(filePath):

    gene = json.loads(open(filePath).read())
    # Data with only results including atc_classifications
    gene_data = run(gene)
    # Data with sugeno averages for items with same atc_name, also has count
    gene_stats = run_stats(gene_data)
    print("Count and average for names:")
    print(gene_stats)
    # Each atc_code in the set, with its list of respective sugeno scores
    atc_codes = atc_dict(gene_data)

    # KS Test of Significance
    ks_test = run_ks(gene, atc_codes)
    ks_df = pd.DataFrame.from_dict(ks_test, orient='index')
    print("Dataframe Results of KS Test:")
    print(ks_df)

    atc_code_index = {}
    index = 0
    for code in atc_codes:
        atc_code_index[index] = (code, atc_codes[code]["names"], atc_codes[code]["scores"])
        index += 1

    #Analysis of Variance(Anova) Significance Test
    atc_code_df = pd.DataFrame.from_dict(atc_code_index, orient='index', columns=['code', 'names', 'scores'])
    anova_test = anova(atc_code_df)
    print("Anova Results:")
    print(anova_test)
    print()

def take_query(api):
    url = "https://ars.test.transltr.io/ars/api/messages/" + api + "?trace=y"
    r = requests.get(url, allow_redirects=True)
    time.sleep(600)

    json_data = r.json()
    #print(json_data['merged_version'])
    merged_version_pk = json_data['merged_version']
    url2 = "https://ars.test.transltr.io/ars/api/messages/" + merged_version_pk
    r2 = requests.get(url2, allow_redirects=True)
    open('results.json', 'wb').write(r2.content)
    json_final = r2.json
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


def anova(atc):
    score_list = []
    for item in atc["scores"]:
        score_list.append(item)
    output = f_oneway(*score_list)
    return output

def run_ks(message, atc_res):
    results = get_safe(message, "fields","data","message","results")
    N = len(results)
    result_scores = []
    for result in results:
        if "sugeno" in result.keys():
            result_scores.append(result["sugeno"])
        else:
            print("no sugeno")

    ks_list = {}
    for code,content in atc_res.items():
        #print(content["names"])
        results_ks = kstest(result_scores, content["scores"])

        ks_list[code] = {
            "results": results_ks,
            "names": content["names"]
        }
    return ks_list
def run_stats(df):
    stats = {}
    temp = {}
    i = 0
    for row in df.iterrows():
        print()
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
def run(data):
    results = getResults(data)
    knowledgeGraph = getKnowledgeGraph(data)
    idlist = getresultIDs(results)
    nodeMap = {}
    for zids in idlist:
        for zid in zids:

            nodeMap[zid] = getNode(knowledgeGraph, zid)


    atcMap = getMaps(nodeMap)
    pr = processResults(results, atcMap)
    scoreSet = showData(pr)
    dataframe = loadPandas(scoreSet)
    return dataframe




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
