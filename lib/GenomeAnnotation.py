try:
    import json
except ImportError:
    import sys
    sys.path.append('simplejson-2.3.3')
    import simplejson as json
    
import urllib



class GenomeAnnotation:

    def __init__(self, url):
        if url != None:
            self.url = url

    def genomeTO_to_reconstructionTO(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.genomeTO_to_reconstructionTO',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomeTO_to_feature_data(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.genomeTO_to_feature_data',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def reconstructionTO_to_roles(self, reconstructionTO):

        arg_hash = { 'method': 'GenomeAnnotation.reconstructionTO_to_roles',
                     'params': [reconstructionTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def reconstructionTO_to_subsystems(self, reconstructionTO):

        arg_hash = { 'method': 'GenomeAnnotation.reconstructionTO_to_subsystems',
                     'params': [reconstructionTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def annotate_genome(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.annotate_genome',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def call_RNAs(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.call_RNAs',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def call_CDSs(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.call_CDSs',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def find_close_neighbors(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.find_close_neighbors',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def assign_functions_to_CDSs(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.assign_functions_to_CDSs',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def annotate_proteins(self, genomeTO):

        arg_hash = { 'method': 'GenomeAnnotation.annotate_proteins',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None




        
