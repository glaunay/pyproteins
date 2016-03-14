#import customCollection

import pyproteins.services.psipredServ as psipredServ
import pyproteins.services.utils as utils

import msa
import json
import os
import random
import re
import collections
import uuid
from Bio.Blast import NCBIXML
from subprocess import call
from Bio import SeqIO


'''
    WARNING Access to a single ss2 property can rely on a psipred prediciton which can be long.
    Prefered way is to iterate over all entries and multi-thread the ss2 property call
'''


'''
Peptide entry are the fundamental entity of the custom NW implementation

'''
ppServ = None

def startPpServ():
    pp = psipredServ.Socket()
    return pp

'''
 indexing use a single catalogue file, to map file name with objects accessors

 listing all element i a single file

'''

# MAin diff w/ previosu collection is that a failed get cannot be resolved from a simple key
# Accessor can be a sequence or even better a Peptide object
# The idea is to store enriched peptide object
# Which can later be accessed using poor (unenriched) peptide Object as accessors


aminoAcidTable = {
    'A' : ['ALA'],
    'C' : ['CYS'],
    'D' : ['ASP'],
    'E' : ['GLU'],
    'F' : ['PHE'],
    'G' : ['GLY'],
    'H' : ['HIS'],
    'I' : ['ILE'],
    'K' : ['LYS'],
    'L' : ['LEU'],
    'M' : ['MET', 'MSE'],
    'N' : ['ASN'],
    'P' : ['PRO'],
    'Q' : ['GLN'],
    'R' : ['ARG'],
    'S' : ['SER'],
    'T' : ['THR'],
    'V' : ['VAL'],
    'W' : ['TRP', 'TRY'],
    'Y' : ['TYR'],
    'X' : ['UNK']
}

def threeToOne(aa):
    for l in aminoAcidTable:
        if aa in  aminoAcidTable[l]:
            return l;
    print "Warning 3 letter code " + aa + " not found";
    return 'X';

def oneToThree(aa):
    if aa in aminoAcidTable:
        return aminoAcidTable[aa][0]
    print "Warning one letter code " + aa + " not found";
    return 'UNK';

class EntrySet(object): #customCollection.EntrySet
    def __init__(self, name=None, dataFile=None):
        self.data = collections.OrderedDict()
        #customCollection.EntrySet.__init__(self, collectionPath="/Users/guillaumelaunay/work/data/peptides", constructor=Entry, typeCheck=isValidID, indexer=strip)
        if not dataFile and not name:
            print "flat file storage missing, Can't creating empty set w/out a name"
            return
        if name:
            self.name = name
            return
        self._parse(dataFile)

    def pluck(self):
        k = random.choice(self.data.keys())
        return self.data[k]

    def mash(self):
        m = hash(self.name + str(len(self)))
        return str(abs(m))

    def _parse(self, dataFile):
        f = open(dataFile, "r")
        dBuffer = json.load(f)
        self.name = dBuffer['id']

        for d in dBuffer['data']:
            e = Entry(datum=d)
            if hash(e) in self.data: # hash string more efficient
                print "Mutliple definition of petide fragment \""+ d['seq'] + "\""
            self.data[hash(e)] = e

    def __len__(self):
        return len(self.data)
    def __repr__(self):
        string = "Peptide set \"" + self.name + "\" ("  + str(len(self)) + " elements)\n"
        for e in self:
            string += e.id + "\n"
        return string

    def enrich(self, _blankShotID=None, *kwargs):
        if not ppServ:
            ppServ = startPpServ()
        stash = [d for d in self if not d.ss2]
        jobid = ppServ.push(peptidesList=stash, _blankShotID=_blankShotID)
        data = ppServ.pull(jobid)

        for i,d in enumerate(stash):
            d.ss2Bind(data[i])

    def __iter__(self):
        for k in self.data:
            yield self.data[k]

    def __getitem__(self, k):
        for x in self.data:
            if self.data[x].id == k:
                return self.data[x]
        return None

    def index(self, string=None, peptideObject=None):
        if string:
            try:
                return [ e.id for e in self.data.values() ].index(string)
            except ValueError:
                return None
        if peptideObject:
            return [ e.id for e in self.data.values() ].index(peptideObject.id)


    def get(self, peptideObject):
        k = hash(peptideObject)
        if k in self.data:
            return self.data[k]
        return None

    def add(self, peptideObject):
        if self.delete(peptideObject):
            print "Following peptide already part of the set"
            print peptideObject
        self.data[hash(peptideObject)] = peptideObject

    def delete(self, peptideObject):
        k = hash(peptideObject)
        if k in self.data:
            del self.data[k]
            return True
        return False

    def serialize(self, targetFile=None, comments="peptide collection"):
        dataOut = {'id' : self.name, 'comments' : comments, 'data' : [ e.dict for e in self ]}
        asJson =  json.JSONEncoder().encode(dataOut)
        if targetFile:
            with open(targetFile, "w") as f:
                f.write(asJson)
        return asJson

    @property # return list of peptides name
    def labels(self):
        return [ e.id for e in self ]




class Entry(object):

    def __getstate__(self): return self.__dict__
    def __setstate__(self, d): self.__dict__.update(d)

    def parse(self, fastaFile):
        handle = open(fastaFile, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        self.id = records[0].id  #first record
        self.seq = str(records[0].seq)  #first record


    def __hash__(self):
        if not self.seed:
            random.seed()
            #self.seed = os.urandom(5)
            self.seed = random.randint(1,1000000)
        return hash(self.seq + str(self.seed))


    def blast(self, xmlPsiBlastOutput): # return an msa Object

        def hit_check(hit):
            for i in range(0, len(hit.hsps) - 1):
                for j in range(i + 1, len(hit.hsps)):
                    if hit.hsps[i].sbjct_start >= hit.hsps[j].sbjct_start:
                        if hit.hsps[i].sbjct_end <= hit.hsps[j].sbjct_end:
                            raise ValueError( "hsp overlap error in hit named " + hit.title)


        def hsp_merge(hsp, array):

            pos = hsp.query_start - 1
            for i in range(0, hsp.align_length):
                if hsp.query[i] == '-':
                    continue
                array[pos] = hsp.sbjct[i]
                pos += 1
            #return ''.join(mergedStrand)


        array=[{ 'id' : self.id,  'seq' : list(self.aaSeq) }]

        if not utils.which("blastpgp"):
            raise initError("Cant find blastpgp executable")

        rootId = uuid.uuid4()
        with open (str(rootId) + '.fasta', 'w') as fOut:
            fOut.write(self.fasta)

        if not xmlPsiBlastOutput:
            call(['blastpgp', '-j', '2', '-i', str(rootId) + '.fasta', '-d', 'nr70', '-m', '7', '-o', str(rootId) + '.xml'])
            f = open(str(rootId) + '.xml')
            records = NCBIXML.parse(f)
        else :
            f = open(xmlPsiBlastOutput)
            records = NCBIXML.parse(f)

        recordList = [x for x in records]
        for psiPass in reversed(recordList):
            if not psiPass.alignments:
                continue
            for hit in psiPass.alignments:
                if len(hit.hsps) > 1:
                    hit_check(hit)
                if hit.hsps[0].sbjct == self.aaSeq:
                    print "Self hit, skipping"
                    continue
                title = str(hit.title) + '[' + ','.join([ str(hsp.sbjct_start) + '-' + str(hsp.sbjct_end) for hsp in hit.hsps ]) + ']'
                datum = { 'id' : title,  'seq' : ['-'] * len(self) }
                for hsp in hit.hsps:
                    hsp_merge(hsp, datum['seq'])
                array.append(datum)

        f.close()
        #return [[ hit['seq'] for hit in array ], [ hit['id'] for hit in array ]]
        msaBean = Msa.MsaBean([ hit['seq'] for hit in array ], [ hit['id'] for hit in array ])
        msa = Msa.Msa(msaBean=msaBean)
        msa.set_id(self.id)
        #return msaBean
        return msa


    def mash(self):
        m = hash(self)
        return str(abs(m))

    def serialize(self):
        return json.JSONEncoder().encode(self.dict) # Check that function name does not create None attribute
        pass # write json formated data structure call be customCollection

    def __init__(self, datum=None, **kwargs):
        self.description = None
        self.ss2Obj = None
        self.ss2Seq = None
        self.seq = None
        self.id = None
        self.seed = None # we need a random generator, bc sequence could be identical

        if datum:
            if 'id' not in datum or 'seq' not in datum:
                raise ValueError("Attribute missing for peptide definition")
            self.id = datum['id']
            self.seq = datum['seq']
            self.description = datum['desc']
            if 'ss2' in datum:
                self.ss2Obj = psipredServ.collection(stream=datum['ss2'])
        else:
            if 'id' in kwargs:
                self.id = kwargs['id']
            if 'seq' in kwargs:
                self.seq = kwargs['seq']
            if 'desc' in kwargs:
                self.description = kwargs['desc']
       #     if 'ss2' in kwargs:
       #         self.ss2 = kwargs['ss2']

        #if not self.seq or not self.id:
        #    raise ValueError("Cant create peptide object")


    def __iter__(self):
        for i in range (0, len(self.seq)):
            yield self[i + 1]

    def __getitem__(self, i):
        if i < 0:
            raise ValueError("minimal amino-acid number is 1")
        v = {'aa' : self.seq[i], 'ss2' : None, 'burial' : None }
        if self.ss2:
            v['ss2'] = self.ss2[i]
        return Msa.Position(v)

    @property
    def hasSse(self):
        if self.ss2:
            return True
        return False

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        s = self.fasta
        if self.ss2Obj:
            s += "\n" + self.ss2Obj.horiz
        return s
    @property
    def dict(self):
        ss2 = None
        if self.ss2Obj:
            ss2 = self.ss2Obj.dict
        return {
            'id' : self.id,
            'seq' : self.seq,
            'desc' : self.description,
            'ss2' : ss2
        }

    def ss2Bind(self, ss2Obj):
        if ss2Obj.aaSeq != self.aaSeq:
            raise ValueError("amino-acid sequence dont match\n" + ss2Obj.fasta + "\n" + self.fasta)
        self.ss2Obj = ss2Obj

    @property
    def aaSeq(self):
        return self.seq.upper()

    @property
    def fasta(self):
        add = ''
        if self.description:
            add = self.description
        return ">" + self.id + add + "\n" + self.aaSeq

    @property
    def ss2(self):
        if self.ss2Obj:
            return self.ss2Obj.horiz
        if self.ss2Seq:
            return self.ss2Seq
        return None
        #else:
        #    tmp = ppServ.push(peptidesList=[self])
        #    return tmp

    @property # try to extract pfam id
    def pfamID(self):
        m = re.search("(PF[\d]+)", self.id)
        if m:
            return m.groups(1)[0]
        return None

