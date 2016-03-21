import os
import drmaa
import pyproteins.services.utils
import re


def hhSgeDump(queryFilePath, templateFilePath, scriptFile):
    string  ='#!/bin/bash\nprintenv\nhhalign '
    string += '-i ' + queryFilePath + ' -t ' + templateFilePath + ' '
    string += '-o hhAlign.out -M first -nocons -nopred'
    with open(scriptFile, "w") as f:
        f.write(string)
    #/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh

def hhAlignParse(resultFilePath):
    results = { 'Q' : [], 'T' : [] }

    with open (resultFilePath, 'r') as f:
        templateCurrPos = 0
        queryCurrPos = 0
        for l in f:
            m=re.match('^(Q|T)[\s]+[\S]+[\s]+([\d]+)[\s]+([\S]+)[\s]+[\d]+[\s]+\([\d]+\)[\s]*$', l)

            if m:
                bQuery = False
                if m.group(1) == 'Q':
                    queryCurrPos = int(m.group(2))
                    bQuery = True
                if m.group(1) == 'T':
                    templateCurrPos = int(m.group(2))
                for c in m.group(3):
                    bUp = False if c == '-' else True
                    if bQuery:
                        results[m.group(1)].append({'aa' : c, 'pos' : queryCurrPos})
                        queryCurrPos = queryCurrPos + 1 if bUp else queryCurrPos
                    else:
                        results[m.group(1)].append({'aa' : c, 'pos' : templateCurrPos})
                        templateCurrPos = templateCurrPos + 1 if bUp else templateCurrPos
    return results


#def hhAlign(query=None, template=None, bSge=False, workDir=os.getcwd(), hhBinDir=None):
#    query.hhDump(filePath=workDir + '/' + 'query.hhFasta')
#    for i, tObj in enumerate(template):
#        sgeFolder = workDir + '/hhAlign_' + str(i)
#        pyproteins.services.utils.mkdir(sgeFolder)
#        tObj.hhDump(filePath=sgeFolder + '/' + 'template.hhFasta')


def hhAlign(query=None, template=None, bSge=False, workDir=os.getcwd(), hhBinDir=None):

    if isinstance(template, list) and bSge:
        queryFilePath = workDir + '/' + 'query.hhFasta'
        query.hhDump(filePath=queryFilePath)

        with drmaa.Session() as s:
            print('hhAlign drmaa session was started successfully')
            jobListId = []
            jt = s.createJobTemplate()
            for i, tObj in enumerate(template):

                sgeFolder = workDir + '/hhAlign_' + str(i)
                pyproteins.services.utils.mkdir(sgeFolder)

                templateFilePath = sgeFolder + '/' + 'template.hhFasta'
                tObj.hhDump(filePath=templateFilePath)

                sgeScript = sgeFolder + '/hhSgeRun.sh'
                hhSgeDump(queryFilePath, templateFilePath, sgeScript)
                pyproteins.services.utils.chmodX(sgeScript)

                jt.workingDirectory=sgeFolder
                jt.joinFiles = True
                jt.nativeSpecification= "-q short.q";
                jt.outputPath = ':' + sgeFolder
                jt.errorPath = ':' + sgeFolder
                jt.jobEnvironment = {'PATH': os.environ['PATH'], 'HHLIB': os.environ['HHLIB']}
                jt.remoteCommand = sgeScript
                jobListId.append(s.runJob(jt))
                print 'HHalign ---> template ' + tObj.id + ', sgeID : ' + jobListId[-1]

            for i,curjob in enumerate(jobListId):
                print('Collecting job ' + curjob)
                retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                print('Job: {0} finished with status {1}'.format(retval.jobId,
                                                    retval.hasExited))
                sgeFolder = workDir + '/hhAlign_' + str(i)
                datum = hhAlignParse(sgeFolder + '/hhAlign.out')
                hhMask(query, template[i], datum)
                #print datum


    # Loop over query sequence
    # if it is found in data we
    # We load from:to starting position and end position of template PDB structure

def hhMask(query, template, hhDatum):
    qStrand = ''
    tStrand = ''

    bIn = False
    for i, aa in enumerate(query.fasta):
        index = i + 1
        #print ' CHK ' + aa + ' ' + str(i) + " ...\n"
        if str(hhDatum['Q'][0]['pos']) == str(index):
            if hhDatum['Q'][0]['aa'] != aa:
                raise ValueError, " Oups to begin was expecting " + aa + ' at pos ' + str(index) + ' Got ' + str(hhDatum['Q'][0]) + ' instead'
            else:
                bIn = True
                break
        qStrand += aa
        tStrand += '-'

    if not bIn:
        raise ValueError, 'Unable to find starting hh align position in query sequence'

    trail = None
    for i, aa in enumerate(hhDatum['Q']):
        if hhDatum['T'][i]['aa'] != '-' and not template.isPdbDefined(int(hhDatum['T'][i]['pos'])):
            if hhDatum['Q'][i]['aa'] == '-':
                continue
            else:
                tStrand += '-'
        else:
            tStrand += hhDatum['T'][i]['aa']

        qStrand += hhDatum['Q'][i]['aa']

        trail = i

    for i in range(trail + 1, len(query.fasta), 1):
        qStrand += query.fasta[i]
        tStrand += '-'

    qStrand += '*'
    tStrand += '*'

    print '>example of query strand\n' + pyproteins.services.utils.lFormat(qStrand) + '\n'
    print '>example of template strand ' + template.id + '\n' + pyproteins.services.utils.lFormat(tStrand) + '\n'

def pirDump(chainID='A'):
    pass





class HomologyModel(object):
    def __init__(self, aliFile=None, sequence=None, templateArray=None):
        self.templates = []
        pass
    def addTemplate(self, pdbFile):
        self.templates.append(Template(pdbFile))

    def addQuery(self, **kwargs):
        if kwargs is not None:
        #    for key, value in kwargs.iteritems():
        #        print "%s == %s" %(key,value)
            if "sequence" in kwargs:
                self.query = Query(**kwargs)

    def link():
        if not self.msa:
            print "no msa found for query have to blast it, plz wait..."
            self.msa = self.peptide.blast()
        for template in self.templates:
            msa.map(self.msa, template.msa)

