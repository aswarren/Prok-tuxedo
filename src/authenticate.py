import os
import sys
import re
import requests
import urllib
import json
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

LOG = sys.stderr

PatricUser = None

def createQueryGet(api_url=None):
    if api_url == None:
        api_url="https://www.patricbrc.org/api/"
    Session = requests.Session()
    #Session.headers.update({ 'accept': "text/tsv" })
    #Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })
    if not authenticateByEnv(Session):
        authenticateByFile(None, Session)
    return Session


def authenticateByFile(tokenFile=None, Session=None):
    if not tokenFile:
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
    if os.path.exists(tokenFile):
        LOG.write("reading auth key from file %s\n"%tokenFile)
        with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            authenticateByString(tokenString, Session)
        return True
    return False


def authenticateByEnv(Session):
    if os.environ.has_key("KB_AUTH_TOKEN"):
        LOG.write("reading auth key from environment\n")
        authenticateByString(os.environ.get('KB_AUTH_TOKEN'), Session)
        return True
    else:
        return authenticateByFile(None, Session)

def authenticateByString(tokenString, Session):
    Session.headers.update({ 'Authorization' : tokenString })
    if "Authorization" in Session.headers:
        global PatricUser
        PatricUser = Session.headers["Authorization"].split(r"|")[3].split("=")[1]
        LOG.write("Patric user = %s\n"%PatricUser)
