import requests
import json

baseURL = "https://outpacebiotest.benchling.com/api/" #testURL
#baseURL = "https://outpacebio.benchling.com/api/" #realURL

key = #testkey
#key = #realkey

# These are the backend locations of things in benchlingAPI
#part3a = "lib_iH2yDHof" #Define the locations of each folder for each part
#registry_id = "src_oh4knSj0"

# For Test Only!
newseq = {}
#newseq["bases"] = "ATCG"
newseq["folderId"] = "lib_1A6kl8LI" #Test, does this value change for each folder within a project?
newseq["isCircular"] = True
newseq["name"] = "Andrew Test"
newseq["registryId"] = "src_D2ebrtJZ" #Test
newseq["namingStrategy"] = "NEW_IDS"
newseq["schemaId"] = "ts_aPuxhTTQ" #Test

#resistance = {
#          "value": "sfso_0CQFeeLN" #use the variable name for the field
#          }

project = {"value": ["sfso_gVlMamQK"]} #use the variable name for the field

#newseq["fields"] = {"Antibiotic Resistance": resistance, "Project": project} #Need to define each field as a dict

class BadRequestException(Exception):
    def __init__(self, message, rv):
        super(BadRequestException, self).__init__(message)
        self.rv = rv

def getBenchling(path, query):
    request = baseURL+'v2/'+path+query
    r = requests.get(request, auth=(key,"")) #real version
    return r.json()

def postSeqBenchling(bases, name, resistance):
    request = baseURL+'v2/dna-sequences'

    newseq["bases"] = bases
    newseq["name"] = name

    if resistance == 'Kanamycin':
        antibioticId = "sfso_0CQFeeLN" #Double check this value in real registry
        resistance = {"value": antibioticId}
        newseq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    elif resistance == 'Carbenicillin':
        antibioticId = "sfso_0CQFeeLN" #Double check this value in real registry
        resistance = {"value": antibioticId}
        newseq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    else:
        return "Your Antibiotic Resistance Is Non-Standard!"

    r = requests.post(request, json=newseq, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            ),
            r,
        )

    return r.json()

def searchSeqBenchling(bases):
    request = baseURL+'v2-experimental/dna-sequences:search-bases'
    query = {"bases": bases}
    r = requests.post(request, json=query, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            ),
            r,
        )

    templates = r.json()
    return templates['dnaSequences']

#def main():
    """
    newseq = {}
    newseq["bases"] =
    newseq["folderId"] =
    newseq["isCircular"] = True
    newseq["name"] =
    newseq["registryId"] = registry_id
    newseq["namingStrategy"] = "NEW_IDS"

    BenchlingAPI("dna-sequences")
    """

#if __name__ == "__main__":
    #main()
