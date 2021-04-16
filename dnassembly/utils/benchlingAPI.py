import requests

baseURL = 'https://outpacebio.benchling.com/api/v2/'
key = 'sk_w9JBBlWWljMkjBE8znl5KSbOesw0S'
registry_id = 'src_oh4knSj0'

def BenchlingAPI(path, query):
    request = baseURL+path+query
    r = requests.get(request, auth=(key,''))

    return r.json()
