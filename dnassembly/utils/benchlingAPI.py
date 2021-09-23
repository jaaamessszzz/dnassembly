import requests
import json
import os

#baseURL = "https://outpacebiotest.benchling.com/api/" #testURL
#baseURL = "https://outpacebio.benchling.com/api/" #realURL
baseURL = os.environ.get('BENCHLING_URL')

apiVersion = os.environ.get('BENCHLING_API_VERSION')

#key = #testkey
key =  os.environ.get('BENCHLING_API_KEY') #realkey

#Define the common parameters for new DNA plasmid sequences
newSeq = {}
newSeq["isCircular"] = True
newSeq["namingStrategy"] = os.environ.get('SEQ_NAMING_STRAT')

#Real Registry
newSeq["registryId"] = os.environ.get('SEQ_REGISTRY_ID') #Outpace Registry
newSeq["schemaId"] = os.environ.get('SEQ_SCHEMA_ID') #ID code for Plasmid
project = {"value": [os.environ.get('BENCHLING_PROJ_ID')]} #MoClo Project
carb = os.environ.get('CARB') #Real
kan = os.environ.get('KAN') #Real

#Test Registry
#newSeq["folderId"] = "lib_1A6kl8LI" #Test, does this value change for each folder within a project?
#newSeq["registryId"] = "src_D2ebrtJZ" #Test
#newSeq["schemaId"] = "ts_aPuxhTTQ" #Test
#project = {"value": ["sfso_gVlMamQK"]} #Test

#Define the common parameters for new DNA Parts
newPart = {}
newPart["isCircular"] = False
newPart["namingStrategy"] = os.environ.get('PART_NAMING_STRAT')

#Real Registry
newPart["registryId"] = os.environ.get('PART_REGISTRY_ID') #Outpace Registry
newPart["schemaId"] = os.environ.get('PART_SCHEMA_ID') #ID code for DNA Part

#TODO - populate with test registry
#newPart["registryId"] = "" #Outpace Registry
#newPart["schemaId"] = "" #ID code for DNA Part

class BadRequestException(Exception):
    pass

def getBenchling(path, query):
    #request = f"{baseURL}{apiVersion}/{path}{query}"
    request = baseURL+'v2/'+path+query
    r = requests.get(request, auth=(key,"")) #real version

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            )
        )

    return r.json()

def postSeqBenchling(bases, name, assembly_type, assembledFrom='', assembledFromID=''):
    #request = f"{baseURL}{apiVersion}/dna-sequences"
    request = baseURL+'v2/dna-sequences'

    newSeq["bases"] = bases
    newSeq["name"] = name

    #Test
    #antibioticId =  "sfso_0CQFeeLN" #test
    #resistance = {"value": antibioticId}
    #newSeq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    #Real
    if assembly_type == 'cassette': #Set the antibiotic to Carb for Stage 2
        resistance = {"value": carb}

    else: #Set the antibiotic to Kan for both Stage 1 and Stage 3
        resistance = {"value": kan}

    MoClo_Assembled_From = {"value": assembledFrom}
    MoClo_Assembled_From_seqID = {"value": assembledFromID}

    #Set the fields
    newSeq["fields"] = {"Antibiotic Resistance": resistance, "Project": project, "MoClo_Assembled_From": MoClo_Assembled_From, "MoClo_Assembled_From_seqID": MoClo_Assembled_From_seqID}

    if assembly_type == 'cassette': #Set location to the parent Stage 2 folder
        newSeq["folderId"] = 'lib_lIpZ86uz'

    elif assembly_type == 'MC': #Set location to the parent Stage 3 folder
        newSeq["folderId"] = 'lib_hWRdTDLG'

    else: #Set location to somewhere in Stage 1 folder
        partType = assembly_type.strip()

        if partType in ['1']:
            newSeq["folderId"] = 'lib_ZaMHBULI'

        elif partType in ['2a','2b']:
            newSeq["folderId"] = 'lib_QBZvgB3q'

        elif partType in ['3a','3b','3c','3d','3e']:
            newSeq["folderId"] = 'lib_tzkMsLGC'

        elif partType in ['4a','4b']:
            newSeq["folderId"] = 'lib_Y275V89m'

        elif partType in ['5']:
            newSeq["folderId"] = 'lib_RcvOpOZC'

        elif partType in ['6']:
            newSeq["folderId"] = 'lib_QOR8IQ4Y'

        elif partType in ['7']:
            newSeq["folderId"] = 'lib_QmgULI3v'

        else:
            newSeq["folderId"] = 'lib_kCnFLwBS' #Send to the parent Stage 1 folder

    r = requests.post(request, json=newSeq, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            )
        )

    return r.json()

def postPartBenchling(bases, name, partType):
    #request = f"{baseURL}{apiVersion}/dna-sequences"
    request = baseURL+'v2/dna-sequences'

    newPart["bases"] = bases
    newPart["name"] = name

    newPart["folderId"] = 'lib_R5pTYwpd' #Sends part to DNA Part Folder

    if partType.strip() in ['1']:
        DNAtype = {"value": 'sfso_ZIYN45Gl'}

    elif partType.strip() in ['2a','2b']:
        DNAtype = {"value": 'sfso_xvdSPEIX'}

    elif partType.strip() in ['3a','3b','3c','3d','3e']:
        DNAtype = {"value": 'sfso_a5vEtqm5'}

    elif partType.strip() in ['4a','4b']:
        DNAtype = {"value": "sfso_U3jQSIjd"}

    elif partType.strip() in ['5']:
        DNAtype = {"value": 'sfso_gyH8kxAF'}

    elif partType.strip() in ['6','7']:
        DNAtype = {"value": 'sfso_eDKuXLe5'}

    newPart["fields"] = {"Project": project, "DNA Type": DNAtype}
    r = requests.post(request, json=newPart, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            )
        )

    return r.json()

def searchSeqBenchling(bases):
    #request = f"{baseURL}{apiVersion}-beta/dna-sequences:search-bases"
    request = baseURL+'v2-beta/dna-sequences:search-bases'
    query = {"bases": bases,
        "registryId": "src_oh4knSj0", #real
        "schemaId": "ts_itQ9daT4"}

        #Test
        #"registryId": "src_D2ebrtJZ",
        #"schemaId": "ts_aPuxhTTQ"}

    r = requests.post(request, json=query, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            )
        )

    templates = r.json()
    return templates['dnaSequences']

def annotatePartBenchling(seqIDs):
    #request = f"{baseURL}{apiVersion}/dna-sequences:autofill-parts"
    request = baseURL+'v2/dna-sequences:autofill-parts'

    annotateList = {"dnaSequenceIds": seqIDs}
    r = requests.post(request, json=annotateList, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            )
        )
    return
