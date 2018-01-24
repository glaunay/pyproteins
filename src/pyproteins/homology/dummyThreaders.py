# For reading JSON Config file 
import json
from pprint import pprint
import tempfile
import xml.etree.ElementTree as ET
#########################

import os
import sys
sys.path.append('/Users/mbachir/dev/pyproteins/src/')
#import drmaa
import pyproteins.services.utils
import pyproteinsExt.structure.operations as op
import pyproteinsExt.structure.coordinates as PDB
import subprocess
import re

classesPotentials = {
    'LL' : -2.0814,
    'LV' : -1.9197,
    'LI' : -1.9500,
    'LM' : -2.0480,
    'LC' : -2.0827,
    'LA' : 0.0163,
    'LS' : 0.2074,
    'LT' : 0.1433,
    'LP' : 0.1394,
    'LG' : -0.0938,
    'LF' : -1.7221,
    'LY' : -1.7252,
    'LW' : -1.6908,
    'LE' : 0.7525,
    'LD' : 0.7401,
    'LN' : 0.8242,
    'LQ' : 0.7431,
    'LK' : 0.5364,
    'LR' : 0.5910,
    'LH' : 0.4589,
    'VV' : -2.1920,
    'VI' : -1.9238,
    'VM' : -1.9835,
    'VC' : -2.0472,
    'VA' : 0.1146,
    'VS' : 0.2116,
    'VT' : 0.2513,
    'VP' : 0.0121,
    'VG' : -0.1449,
    'VF' : -1.5688,
    'VY' : -1.6748,
    'VW' : -1.7111,
    'VE' : 0.9260,
    'VD' : 0.8314,
    'VN' : 0.7816,
    'VQ' : 0.7946,
    'VK' : 0.5887,
    'VR' : 0.5305,
    'VH' : 0.5580,
    'II' : -2.1643,
    'IM' : -2.0520,
    'IC' : -2.0224,
    'IA' : 0.0864,
    'IS' : 0.1458,
    'IT' : 0.1271,
    'IP' : 0.0603,
    'IG' : -0.1068,
    'IF' : -1.7212,
    'IY' : -1.6617,
    'IW' : -1.7238,
    'IE' : 0.7838,
    'ID' : 0.8180,
    'IN' : 0.8392,
    'IQ' : 0.7753,
    'IK' : 0.5740,
    'IR' : 0.5250,
    'IH' : 0.5086,
    'MM' : -2.0442,
    'MC' : -2.0654,
    'MA' : 0.1399,
    'MS' : 0.1473,
    'MT' : 0.1117,
    'MP' : 0.1207,
    'MG' : -0.1215,
    'MF' : -1.7564,
    'MY' : -1.6629,
    'MW' : -1.7555,
    'ME' : 0.8048,
    'MD' : 0.8040,
    'MN' : 0.7711,
    'MQ' : 0.7623,
    'MK' : 0.5935,
    'MR' : 0.5359,
    'MH' : 0.4874,
    'CC' : -2.0816,
    'CA' : 0.0954,
    'CS' : 0.1533,
    'CT' : 0.1308,
    'CP' : 0.1121,
    'CG' : -0.1710,
    'CF' : -1.7221,
    'CY' : -1.6991,
    'CW' : -1.7483,
    'CE' : 0.7560,
    'CD' : 0.7756,
    'CN' : 0.8008,
    'CQ' : 0.7861,
    'CK' : 0.5313,
    'CR' : 0.5145,
    'CH' : 0.5053,
    'AA' : -0.0485,
    'AS' : -0.0623,
    'AT' : 0.0646,
    'AP' : -0.0113,
    'AG' : 0.0913,
    'AF' : -0.4573,
    'AY' : -0.3710,
    'AW' : -0.4687,
    'AE' : 0.1622,
    'AD' : 0.1807,
    'AN' : 0.1291,
    'AQ' : 0.1605,
    'AK' : -0.1504,
    'AR' : -0.2445,
    'AH' : -0.1913,
    'SS' : -0.0111,
    'ST' : -0.0019,
    'SP' : 0.0359,
    'SG' : 0.1351,
    'SF' : -0.3676,
    'SY' : -0.4782,
    'SW' : -0.4277,
    'SE' : 0.0829,
    'SD' : 0.2486,
    'SN' : 0.0939,
    'SQ' : 0.2113,
    'SK' : -0.2430,
    'SR' : -0.2503,
    'SH' : -0.2954,
    'TT' : 0.0610,
    'TP' : 0.0669,
    'TG' : 0.1655,
    'TF' : -0.4872,
    'TY' : -0.3248,
    'TW' : -0.4496,
    'TE' : 0.2252,
    'TD' : 0.2100,
    'TN' : 0.0822,
    'TQ' : 0.1486,
    'TK' : -0.1127,
    'TR' : -0.2423,
    'TH' : -0.2279,
    'PP' : -0.0482,
    'PG' : 0.0547,
    'PF' : -0.4287,
    'PY' : -0.4705,
    'PW' : -0.4167,
    'PE' : 0.1327,
    'PD' : 0.1466,
    'PN' : 0.2046,
    'PQ' : 0.1412,
    'PK' : -0.2173,
    'PR' : -0.2262,
    'PH' : -0.2584,
    'GG' : 0.1006,
    'GF' : -0.6404,
    'GY' : -0.5133,
    'GW' : -0.5142,
    'GE' : 0.1110,
    'GD' : 0.0053,
    'GN' : 0.0174,
    'GQ' : -0.0019,
    'GK' : 0.4217,
    'GR' : 0.1033,
    'GH' : 0.1919,
    'FF' : -0.6684,
    'FY' : -0.6699,
    'FW' : -0.6378,
    'FE' : 0.3916,
    'FD' : 0.4041,
    'FN' : 0.2861,
    'FQ' : 0.3182,
    'FK' : -0.1673,
    'FR' : -0.3047,
    'FH' : -0.2188,
    'YY' : -0.6301,
    'YW' : -0.6352,
    'YE' : 0.3785,
    'YD' : 0.4749,
    'YN' : 0.3803,
    'YQ' : 0.3841,
    'YK' : -0.1128,
    'YR' : -0.1948,
    'YH': -0.2804,
    'WW' : -0.6509,
    'WE' : 0.3910,
    'WD' : 0.3548,
    'WN' : 0.3760,
    'WQ' : 0.3461,
    'WK' : -0.2514,
    'WR' : -0.2470,
    'WH' : -0.2427,
    'EE' : 0.5848,
    'ED' : 0.6123,
    'EN' : 0.6167,
    'EQ' : 0.5703,
    'EK' : -0.2430,
    'ER' : -0.4305,
    'EH' : -0.2890,
    'DD' : 0.5279,
    'DN' : 0.4841,
    'DQ' : 0.5382,
    'DK' : -0.3138,
    'DR' : -0.2624,
    'DH' : -0.2705,
    'NN' : 0.4959,
    'NQ' : 0.5521,
    'NK': -0.2523,
    'NR' : -0.2412,
    'NH' : -0.2494,
    'QQ' : 0.4399,
    'QK' : -0.1168,
    'QR' : -0.2640,
    'QH' : -0.3106,
    'KK' : 0.6717,
    'KR' : 0.6508,
    'KH' : 0.5253,
    'RR' : 0.4516,
    'RH' : 0.5019,
    'HH' : 0.5265
    }

    # Create a batch file to execute HHAlign in order to align the template full SEQRES Sequence
    # with the Query fasta sequence#

def hhSgeDump(queryFilePath, templateFilePath, scriptFile, hhBinDir, outDir, glob = False):

    # Modifications
    os.environ['HHLIB']='~/Desktop/pipeline/hhsuite/lib/hh/'
    string  = '#!/bin/bash\nexport HHLIB='+os.environ['HHLIB']+'\n'+hhBinDir+'hhalign '
    string += '-i ' + queryFilePath + ' -t ' + templateFilePath + ' '
    string += '-o '+outDir+'/hhAlign.out -M first -nocons -nopred'
    if glob:
        string += ' -glob'
    
    with open(scriptFile, "w") as f:
        f.write(string)
    rq = subprocess.call('sh '+scriptFile ,shell=True, stdout=subprocess.PIPE)
    #print rq.communicate('n\n')[0]
    # End
    #/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh

def psiblast(queryFilePath, dbFilePath, iteration, hhBinDir, scriptFile ,outDir):
    os.environ['HHLIB']='~/Desktop/pipeline/hhsuite/lib/hh/'

    string  = '#!/bin/bash\nexport HHLIB='+os.environ['HHLIB']+'\n'+hhBinDir+'blastpgp '
    string += '-i ' + queryFilePath + ' -d ' + dbFilePath + ' '
    string += '-o '+outDir+'/blast.out -m 7 -j '+iteration
    
    with open(scriptFile+'_blast', "w") as f:
        f.write(string)
    rq = subprocess.call('sh '+scriptFile+'_blast' ,shell=True, stdout=subprocess.PIPE)


    # Parse the Query Sequence et the Template Sequence from the HHAlign output
def hhAlignParse(resultFilePath):
    results = { 'Q' : [], 'T' : [] ,'Coverage' : None}

    with open (resultFilePath, 'r') as f:
        templateCurrPos = 0
        queryCurrPos = 0
        for l in f:
            m=re.match('^(Q|T)[\s]+[\S]+[\s]+([\d]+)[\s]+([\S]+)[\s]+[\d]+[\s]+\([\d]+\)[\s]*$', l)

            if 'Identities' in l:
                results['Coverage'] = l.split('  ')[4].split('=')[1]

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

def hhAlign(query=None, template=None, bSge=False, workDir=os.getcwd(), hhBinDir=None, nModel=0):
    data = []

    if isinstance(template, list) and bSge:

        # Get the first key of my template's seqres attribute corresponding of the first chain ID 

        queryFilePath = query.filePath
        #queryFilePath = workDir + '/' + 'query.hhFasta'
        #query.hhDump(filePath=queryFilePath)

        print('hhAlign drmaa session was started successfully')
        #print template[0]
        #jobListId = []
        #jt = s.createJobTemplate()
        for i, tObj in enumerate(template):
            
            if not len(tObj._seqres) < 1:
                fKey = sorted(list(tObj._seqres))[0]
                fkey = tObj._seqres[fKey]
            else:
                fKey = tObj.chainList[0]

            sgeFolder = workDir + '/' + query.id.split('|')[1] + '_' + fKey
            pyproteins.services.utils.mkdir(sgeFolder)

            data.append({'workDir' : sgeFolder})

            # Create a temporary FASTA file of SEQRES sequence form template (Structure) Object
            templateFilePath = tempfile.NamedTemporaryFile(suffix = '.fasta', dir = workDir , delete=False)
            fasta_file = '>'+tObj.name+'\n'
            # We took the chain A by default
            fasta_file += tObj._seqres[fKey] if tObj._seqres else tObj.chain(fKey).fasta
            templateFilePath.write(fasta_file)
            templateFilePath.close()

            #templateFilePath = '/Users/mbachir/Desktop/dummyTest/template/ADH1A.a3m'
            #templateFilePath = sgeFolder + '/' + 'template.hhFasta'
            #tObj.hhDump(filePath=templateFilePath)

            sgeScript = sgeFolder + '/hhSgeRun.sh'
            hhSgeDump(queryFilePath, templateFilePath.name, sgeScript, hhBinDir, sgeFolder)
            pyproteins.services.utils.chmodX(sgeScript)


            #jt.workingDirectory=sgeFolder
            #jt.joinFiles = True
            #jt.nativeSpecification= "-q short.q";
            #jt.outputPath = ':' + sgeFolder
            #jt.errorPath = ':' + sgeFolder
            #jt.jobEnvironment = {'PATH': os.environ['PATH'], 'HHLIB': os.environ['HHLIB']}
            #jt.remoteCommand = sgeScript
            #jobListId.append(s.runJob(jt))
            #print 'HHalign ---> template ' + tObj.id + ', sgeID : ' + jobListId[-1]

            #for i,curjob in enumerate(jobListId):
            #print('Collecting job ' + curjob)
            #retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            #print('Job: {0} finished with status {1}'.format(retval.jobId, retval.hasExited))

            data[i]['chainStructure'] =  fKey
            data[i]['hhAlignFilePath'] =  data[i]['workDir'] + '/hhAlign.out'
            data[i]['hhAlignResults'] = hhAlignParse( data[i]['hhAlignFilePath'] )

            #data[i]['contactMap'] = checkContactMap(data[i]['hhAlignResults'], query, tObj)

            data[i]['hhAlignStrands'] = hhMask(query, template[i],  data[i]['hhAlignResults'], data[i]['chainStructure'])
            data[i]['pirFilePath'] = data[i]['workDir'] + '/default.pir'
            pirDump(data[i], template=template[i], query=query, chainID=data[i]['chainStructure'],filePath=data[i]['pirFilePath'], )


        """
        with drmaa.Session() as s:
            print('hhAlign drmaa session was started successfully')
            print template
            jobListId = []
            jt = s.createJobTemplate()
            for i, tObj in enumerate(template):

                sgeFolder = workDir + '/hhAlign_' + str(i)
                pyproteins.services.utils.mkdir(sgeFolder)

                data.append({'workDir' : sgeFolder})

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

                data[i]['hhAlignFilePath'] =  data[i]['workDir'] + '/hhAlign.out'
                data[i]['hhAlignResults'] = hhAlignParse( data[i]['hhAlignFilePath'] )
                data[i]['hhAlignStrands'] = hhMask(query, template[i],  data[i]['hhAlignResults'])
                data[i]['pirFilePath'] = data[i]['workDir'] + '/default.pir'
                pirDump(data[i], template=template[i], query=query, chainID='A',filePath=data[i]['pirFilePath'])
    """
        if nModel > 0:
            modellerAll(query, template, data, nModel=nModel)

                #print datum

def checkContactMap(query = list, templatePDB = list, hhAlignFiles = list , hhBinDir=None):

    hhDatum = []

    if len(query) == len(templatePDB):
        '''
        for i , (q, t) in enumerate(zip(query, templatePDB)):

            if t._seqres:
                fKey = sorted(list(t._seqres))[0]
                fkey = t._seqres[fKey]
            else:
                fKey = t.chainList[0]

            sgeFolder = workDir + '/hhAlign_' + str(i) +'_' + q.id.split('|')[1]
            pyproteins.services.utils.mkdir(sgeFolder)

             # Create a temporary FASTA file of SEQRES sequence form template (Structure) Object
            templateFilePath = tempfile.NamedTemporaryFile(suffix = '.fasta', dir = sgeFolder , delete=False)
            fasta_file = '>'+t.name+'\n'
            # We took the chain A by default
            fasta_file += t._seqres[fKey] if t._seqres else t.chain(fKey).fasta
            templateFilePath.write(fasta_file)
            templateFilePath.close()

            sgeScript = sgeFolder + '/hhSgeRun.sh'
            hhSgeDump(q.filePath, templateFilePath.name, sgeScript, hhBinDir, sgeFolder, True)
            pyproteins.services.utils.chmodX(sgeScript)

            hhDatum.append(hhAlignParse( sgeFolder + '/hhAlign.out'))
        '''

        #parser = PDB.Parser()
        templateSeq1 = templatePDB[0]
        templateSeq2 = templatePDB[1]

        parseSeq1 = hhAlignParse(hhAlignFiles[0])
        parseSeq2 = hhAlignParse(hhAlignFiles[1])

        querySeq1Fasta = [{aa['pos'] : aa['aa'] for aa in parseSeq1['Q']}]
        querySeq2Fasta = [{aa['pos'] : aa['aa'] for aa in parseSeq2['Q']}]

        templateSeq1Fasta = [{aa['pos'] : aa['aa'] for aa in parseSeq1['T']}]
        templateSeq2Fasta = [{aa['pos'] : aa['aa'] for aa in parseSeq2['T']}]

        matrixTemplate = op.ContactMap(templateSeq1, templateSeq2)
        positionMatches = matrixTemplate.positionOfContact(matrixTemplate.residuesInterfacialBool(9, False))

        energyT = 0
        energyQ = 0
        nbrGapQa = 0
        nbrGapQb = 0

        for p in positionMatches:
            
            #Count Nbr Gap
 
            if p[0] >= min(querySeq1Fasta[0].keys()) and p[0] <= max(querySeq1Fasta[0].keys()) and p[1] >= min(querySeq2Fasta[0].keys()) and p[1] <= max(querySeq2Fasta[0].keys()):
                if querySeq1Fasta[0][p[0]] == '-':
                    nbrGapQa += 1
                elif querySeq2Fasta[0][p[1]] == '-':
                    nbrGapQb += 1

                if str(querySeq1Fasta[0][p[0]]+querySeq2Fasta[0][p[1]]) in classesPotentials:
                    if str(querySeq1Fasta[0][p[0]]) != "X" and str(querySeq2Fasta[0][p[1]]) != "X":
                        energyQ += classesPotentials[str(querySeq1Fasta[0][p[0]]+querySeq2Fasta[0][p[1]])]
                else:
                    if str(querySeq1Fasta[0][p[0]]) != "X" and str(querySeq2Fasta[0][p[1]]) != "X":
                        energyQ += classesPotentials[str(querySeq2Fasta[0][p[1]])+querySeq1Fasta[0][p[0]]]

            if p[0] >= min(templateSeq1Fasta[0].keys()) and p[0] <= max(templateSeq1Fasta[0].keys()) and p[1] >= min(templateSeq2Fasta[0].keys()) and p[1] <= max(templateSeq2Fasta[0].keys()):
                if templateSeq1Fasta[0][p[0]] == '-':
                    nbrGapQa += 1
                elif templateSeq2Fasta[0][p[1]] == '-':
                    nbrGapQb += 1

                if str(templateSeq1Fasta[0][p[0]]+templateSeq2Fasta[0][p[1]]) in classesPotentials:
                    if str(templateSeq1Fasta[0][p[0]]) != "X" and str(templateSeq2Fasta[0][p[1]]) != "X":
                        energyT += classesPotentials[str(templateSeq1Fasta[0][p[0]]+templateSeq2Fasta[0][p[1]])]
                    
                else:
                    if str(templateSeq1Fasta[0][p[0]]) != "X" and str(templateSeq2Fasta[0][p[1]]) != "X":
                        energyT += classesPotentials[str(templateSeq2Fasta[0][p[1]]+templateSeq1Fasta[0][p[0]])]
        
        print 'Nbr Contact : '+str(len(positionMatches))+'\n'+'Energy Template : ' +str(energyT) + '\nEnergy Query : ' +str(energyQ) + '\nNumber Gap Q/T A : ' +str(nbrGapQa) + '\nNumber Gap Q/T B : ' +str(nbrGapQb)


def hhMask(query, template, hhDatum, chain):
    
    # Create 2 Strands:
    # - Query strand (qStrand): Sequence from the HHAlign output
    # - Template strand (tStrand) : Resultat of comparision position per position of amino acid nature bewteen
    #                               template sequence from HHAlign alignement and the fasta sequence of the
    #                               ATOM Record fields for a given chain in the PDB Sructure.
    #                               If we have similarity for a given position we store the amino acid 
    #                               If not we store a gap '-'
    #       In order to pirDump function

    qStrand = ''
    tStrand = ''
    bIn = False

    for i, aa in enumerate(query.fasta):
        index = i + 1
        #print ' CHK ' + aa + ' ' + str(i) + " ...\n"
        if str(hhDatum['Q'][0]['pos']) == str(index):
            print str(hhDatum['Q'][0]['pos']), str(index)
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
        # isPdbDefined() get index position of an amino acid (in sequence ?)

        if hhDatum['T'][i]['aa'] != '-' and not template.hasCoordinates(hhDatum['T'][i], chain):
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

    return { 'Q' : qStrand, 'T' : tStrand }

def pirDump(datum, template=None, query=None, chainID='A',filePath=None):
    
    # Generate a PIR alignement file format giving the "structure" sequence --> tStrand
    #                                              the "query" Sequence     --> qStrand

    strands = datum['hhAlignStrands']
    hhDatum = datum['hhAlignResults']

    #Replace all X caracters in sequence with a dot "." --> X corresponding HETATM field.
    # Modeller read a HETATM record field by a dot in both Q|T Strands
    # Need to replace by a dot at the same position as the template to the query strand
    indexDots = [index for index in template._hetatm]
    
    if int(hhDatum['Q'][0]["pos"]) >= int(hhDatum['T'][0]["pos"]):
        shift =  int(hhDatum['Q'][0]["pos"]) - int(hhDatum['T'][0]["pos"])
    else:
        shift = (int(hhDatum['T'][0]["pos"]) - ( int(hhDatum['T'][0]["pos"]) - int(hhDatum['Q'][0]["pos"]))) + 1

    if len(indexDots) > 0:
        strands['T'] = list(strands['T']) 
        strands['Q'] = list(strands['Q'])

        for index in indexDots:
        # Replacing in Q/T Strand the HETATM
            
            for i, aa in enumerate(strands['T']):
                
                limit = None

                if index-1 <= hhDatum['T'][-1]["pos"]-1:
                    limit = index - 1
                else:
                    limit = hhDatum['T'][-1]["pos"]-1
                    
                if i >= hhDatum['Q'][0]["pos"] and i <= limit :
                    if aa == "-":
                        shift += 1

            print "Limit T : "+str(hhDatum['T'][-1]["pos"]-1) + " Shift : "+str(shift - 1) + " Index : " +str(index-1) + " Q min : " + str(hhDatum['Q'][0]["pos"]) + " T min : " + str(hhDatum['T'][0]["pos"])

            #if index > hhdatum[T] si T > Q en terme de pos[0]

            if index-1 < hhDatum['T'][-1]["pos"]-1 and index-1 > hhDatum['Q'][0]["pos"]: # -1 for len(strand['T']) to keep the "*" character i, pir alignement
                strands['T'][index-1 + shift] = "."
            
            if index-1 < hhDatum['T'][-1]["pos"]-1 and index-1 > hhDatum['Q'][0]["pos"]:
                strands['Q'][index-1 + shift] = "."

            #Reset SHIFT
            shift =  int(hhDatum['Q'][0]["pos"]) - int(hhDatum['T'][0]["pos"])

        strands['Q'] = ''.join(strands['Q'])
        strands['T'] = ''.join(strands['T'])
    
    filePath = filePath if filePath else os.getcwd + 'default.pir'

    i=0
    while not template.hasCoordinates(hhDatum['T'][i], chainID) :
        i += 1

    j = len(hhDatum['T']) - 1
    print hhDatum['T'][j]
    while not template.hasCoordinates(hhDatum['T'][j], chainID) :
        j -= 1


    #print '-->' + str(i) + ' == ' + str(j) + "<--"
    #print '-->' + str(hhDatum['T'][i]['pos']) + ' == ' + str(hhDatum['T'][j]['pos']) + "<--"

    templateStart = hhDatum['T'][i]['pos']
    templateStop  =  hhDatum['T'][j]['pos']

    pirContent = '>P1;' + template.name + '\nstructure:' + template.filePath + ':' + str(templateStart) + ' : ' + chainID
    pirContent += ' : ' + str(templateStop) + ' : ' + chainID + '::::\n'
    pirContent += pyproteins.services.utils.lFormat(strands['T']) + '\n\n'
    # Modif
    pirContent += '>P1;' + query.id + '\nsequence:' + query.id #Condition
    pirContent += ':1: :' + str(len(query.fasta)) +': ::::\n'
    # End
    pirContent += pyproteins.services.utils.lFormat(strands['Q']) + '\n'

    with open (filePath, 'w') as f:
        print os.getcwd()
        f.write(pirContent)


def modellerAll(query, template, data, nModel=1):

     #with drmaa.Session() as s:
        print('modeller drmaa session was started successfully')
            #jobListId = []
            #jt = s.createJobTemplate()
        for i, tObj in enumerate(template):
            workDir = data[i]['workDir']

            sgeScript = workDir + '/runModeller.sh'
            ## Open config.json and fetch the actual call command to lunch Modeller via command line
            with open('/Users/mbachir/dev/pyproteins/src/pyproteins/homology/config.json') as data_file:    
                modellerPath = json.load(data_file)
                modellerPath = str(modellerPath['soft']['modeller']['command_line_call'])

            with open(sgeScript, "w") as f:
                f.write(modellerPath+' '+workDir+'/doModel.py > modeller.log')
            pyproteins.services.utils.chmodX(sgeScript)

            pyString = '# Comparative modeling by the automodel class\n'
            pyString += 'from modeller import *              # Load standard Modeller classes\n'
            pyString += 'from modeller.automodel import *    # Load the automodel class\n'
            pyString += 'log.verbose()    # request verbose output\n'
            pyString += 'env = environ()  # create a new MODELLER environment to build this model in\n'
            pyString += '# directories for input atom files\n'
            pyString += 'env.io.atom_files_directory = [\'' + os.path.dirname(tObj.filePath) + '\']\n'
            pyString += '# Read in HETATM records from template PDBs\n'
            pyString += 'env.io.hetatm = True\n'
            pyString += 'a = automodel(env,\n'
            pyString += '              alnfile  = \'' + data[i]['pirFilePath'] + '\',     # alignment filename\n'
            pyString += '              knowns   = \'' + tObj.name + '\',              # codes of the templates\n'
            pyString += '              sequence = \'' + query.id + '\')              # code of the target\n'
            pyString += 'a.starting_model= 1                 # index of the first model\n'
            pyString += 'a.ending_model  = ' + str(nModel) + '                 # index of the last model\n'
            pyString += '                                    # (determines how many models to calculate)\n'
            pyString += 'a.make()                            # do the actual comparative modeling\n'
            with open(workDir + '/doModel.py', "w") as f:
                f.write(pyString)

            req = subprocess.Popen('cd ' + workDir + ';sh '+workDir+'/runModeller.sh' ,shell=True, stdout=subprocess.PIPE)
            print req.communicate('n\n')[0]
            
            #jt.workingDirectory = workDir
            #jt.joinFiles = True
            #jt.nativeSpecification= "-q short.q";
            #jt.outputPath = ':' + workDir
            #jt.errorPath = ':' + workDir
            #jt.jobEnvironment = {'PATH': os.environ['PATH'], 'KEY_MODELLER6v2' : 'MODELIRANJE' }
            #jt.remoteCommand = sgeScript
            #jobListId.append(s.runJob(jt))

            #for i,curjob in enumerate(jobListId):
                #print('Collecting modeller job ' + curjob)
                #retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                #print('Job: {0} finished with status {1}'.format(retval.jobId,
                                                    #retval.hasExited))


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
