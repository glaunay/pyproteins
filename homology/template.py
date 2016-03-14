from Bio.PDB import *
import Peptide
import copy
import Msa
import os.path
#from Bio.PDB.Polypeptide import PPBuilder

'''
    One or several pdb files
    A pepidic sequence or msa

    Do a msa on the pdb extracted sequences
    Do a msa on the peptidic sequence.
    Do a NW alignement between the two
'''


PDBparser = PDBParser()

def fastaFileToList(file):
    comment = ''
    data = []
    with open(file, 'r') as input:
        for line in input:
            if line.startswith(">"):
                comment += line
            else:
                data += line.split()
    return { 'header' : comment, 'data' : data }



class HomologyModel:
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
            Msa.map(self.msa, template.msa)

class Query: # Possible input, aa sequence string or mfasta file or mfasta bean
    def __init__(self, **kwargs):
        if kwargs is not None:
            if 'sequence' in kwargs:
                print "loading peptide sequence " + kwargs['sequence']
                self.peptide = Peptide.Entry(id='homol-mdl', seq=kwargs['sequence'])
            elif 'fastaFile' in kwargs:
                self.peptide = Peptide.Entry()
                self.peptide.parse(kwargs['fastaFile'])
            if 'msaFile' in kwargs:
                print "loading msa from file"
                self.msa = Msa(fileName=kwargs['msaFile'])



class Template:

    def __repr__(self):
        return 'template string'

    def __init__(self, pdbSource, modelID = None, chain = None, folder=None, id = None):
        self.structure = None
        self.pdbSeq = None
        self.mAli = None
        self.pdbnum = None
        self.see = None
        self.fasta = None
        self.folder = folder
        self.pdbSource = PDBparser.get_structure('mdl', pdbSource)
        model = self.pdbSource[0] if not modelID else self.pdbSource[modelID]
        #By default we use first chain atom coordinates record, user defined alternative chain w/
        # chain argument
        chainIdSorted = [ key for key, value in sorted(model.child_dict.items()) ]
        self.structure = model[chainIdSorted[0]] if not chain else model[chain]
        self.id = id if id else os.path.basename(pdbSource)

        # Extract CA sequence and compute pairwise dist
        self.pdbSeq = [ r['CA'] for r in self.structure if 'CA' in r ]
        self._setFromFolder()
        print 'Done'

    ## initialize object attributes from a make core folder if any provided
    def _setFromFolder(self):
        if not self.folder:
            return False
        # pdbnum
        for file in os.listdir(self.folder):
            #print file
            if file.endswith('.pdbnum'):
                data = fastaFileToList(self.folder + '/' + file)
                self.pdbnum = data['data']
            if file.endswith('.fasta'):
                data = fastaFileToList(self.folder + '/' + file)
                print file
                self.fasta = data['data']
                #print self.fasta

    @property
    def aaSeq(self):
        if self.fasta:
            return self.fasta
        #for atom in self.sequence
        return [ (Peptide.threeToOne(atom.get_parent().resname)) for atom in self.pdbSeq ]

    def msa(self, psiBlastOutputXml=None):
        if not self.mAli:
            peptide = Peptide.Entry(id = "PDB template sequence", seq = ''.join(self.aaSeq))
            self.mAli = peptide.blast(psiBlastOutputXml)
        return self.mAli


