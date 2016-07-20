import paramiko
import re
import uuid
import subprocess
import os
import json
from cStringIO import StringIO
import pyproteins.services.utils
import time

CONF = {
    "migale" : {
        "socket" : "/projet/extern/save/gulaunay/tmp/socket/ppsvr",
        "host" : "migale.jouy.inra.fr",
        "client" : "gulaunay"
    },
    "arwen"  : {
        "socket" : "/home/glaunay/tmp/socket/ppsvr",
        "host" : "arwen.ibcp.fr",
        "client" : "glaunay"
    }
}


INPUT_FORMAT='^[\s]*([\d]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([0-9\.]+)[\s]+([0-9\.]+)[\s]+([0-9\.]+)'




def dumpSlurmBatch(location, njobs, app="psiblast"):

    if app == "psiblast":
        string = "#!/bin/bash\n#SBATCH --job-name=psiPredArray\n#SBATCH -p express-mobi\n#SBATCH --qos express-mobi\n#SBATCH --output=psiPredArray%A_%a.out\n\n#SBATCH --error=psiPredArray%A_%a.err\n"
        string += "#SBATCH --array=1-" + str(njobs) + "\n#SBATCH --time=0:05:00\n#SBATCH --workdir=" + location + "\n#SBATCH --ntasks=1\n#SBATCH -N 1\n"
        string += "~/runpsipred -i " + location + "/peptide_\$SLURM_ARRAY_TASK_ID.fasta -d ~/db/uniprot_sprot_current -r " + location + "/$SLURM_ARRAY_TASK_ID\n"
    elif app == "blastpgp":
        string = "#!/bin/bash\n#SBATCH --job-name=psiBlastArray\n#SBATCH -p express-mobi\n#SBATCH --qos express-mobi\n#SBATCH --output=" + location + "/psiBlastArray%A_%a.out\n\n#SBATCH --error=" + location + "/psiBlastArray%A_%a.err\n"
        string += "#SBATCH --array=1-" + str(njobs) + "\n#SBATCH --time=0:05:00\n#SBATCH --wait\n#SBATCH --workdir=" + location + "\n#SBATCH --ntasks=1\n#SBATCH -N 1\n"
        #string += "echo blastpgp -j 2 -m 7 -i " + location + "/peptide_\$SLURM_ARRAY_TASK_ID.fasta -d ~/db/uniprot_swissprot_current -o " + location + "/peptide_\$SLURM_ARRAY_TASK_ID.blast\n"
        string += "module load /software/mobi/modules/ncbi-blast/2.2.26\n"
        fastaInput = location + "/peptide_\$SLURM_ARRAY_TASK_ID.fasta"
        outputFile = location + "/peptide_\$SLURM_ARRAY_TASK_ID.blast"
        string += "/software/mobi/ncbi-blast/ncbi-blast-2.2.26/bin/blastpgp -b 500 -j 3 -m 7 -i " +  fastaInput + " -d ~/db/uniprot_sprot_current -o " + outputFile + "; touch " + location + "/peptide_\$SLURM_ARRAY_TASK_ID.done\n"
        #string += "/software/mobi/ncbi-blast/ncbi-blast-2.2.26/bin/blastpgp -j 1 -m 7 -i " +  fastaInput + " -d ~/db/nr/nr70 -o " + outputFile + "; touch " + location + "/peptide_\$SLURM_ARRAY_TASK_ID.done\n"

    return string

def isValidRecord(line):
   # print line
    m = re.match(INPUT_FORMAT, line)
    if m:
        return True
    return False






class Socket(object):
    def __init__(self, service="migale", app="psiblast", **kwargs):
        self.client = paramiko.client.SSHClient()
        self.client.load_system_host_keys()
        self.user = CONF[service]["client"]
        self.host = CONF[service]["host"]
        self.socket = CONF[service]["socket"]
        self.serviceName = service
        self.pool = {}
        self.localCache = None
        self.app = app # Intended to specify blastpgp or psiblast
        self.client.connect(hostname=self.host, username=self.user)
        if 'localCache' in kwargs:
            self.localCache = kwargs['localCache']


    def close(self) :
        self.client.close()

    def push(self, fileList=None, peptidesList=None, _blankShotID=None, previous=None):
        if _blankShotID:
            name = _blankShotID
        elif previous:
            name = previous
        else:
            name = uuid.uuid1().get_hex()

        self.pool[name] = {
            'socket' :  self.socket + "/" + name
        }

        if self.localCache :
            self.pool[name]['localCache'] = self.localCache + "/" + name
            if not previous:
                print "creating local cache at " + self.pool[name]['localCache']
            #p = subprocess.Popen(['mkdir ', self.pool[name]['localCache']], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                os.mkdir(self.pool[name]['localCache'])

        if not _blankShotID:
            self.client.exec_command("mkdir " + self.pool[name]['socket'])

        if fileList:
            self.pool[name]['inputs'] = self._pushFiles(name, fileList)
        elif peptidesList:
            self.pool[name]['inputs'] = self._pushPeptides(name, peptidesList)
        else:
            raise ValueError("No input to process")

        if previous:
            print "socket ready to pull previous job " + name
            return name

        #print 'Chunk content ' + str(self.pool[name]['inputs'])
        if self.serviceName == "migale" :
            cmd = 'submitPsipred.sh ' + self.pool[name]['socket']  + ' ' + ' '.join([self.pool[name]['socket'] + '/' + os.path.basename(f) for f in self.pool[name]['inputs'] ])
        elif self.serviceName == "arwen" :
            slurmBatchString = dumpSlurmBatch( self.pool[name]['socket'], len(self.pool[name]['inputs']), app=self.app)
            slurmBatchFile = self.pool[name]['socket'] + '/' + self.app + 'ArraySlurm.sh'
            cmd = "echo -e \"" + slurmBatchString + "\" >" + slurmBatchFile + "\n"
            self.client.exec_command(cmd)
            cmd = "sbatch " + slurmBatchFile

        if _blankShotID:
            print "Not executing remote command:\"" + cmd + "\""
            return name

        print name + " : " + str(len(self.pool[name]['inputs'])) + " " + self.app

        ans = self.client.exec_command(cmd)
        o =  ans[1].read()
        print "---->" + o
        s =  ans[2].read()
        if s:
            raise ValueError(s)

        if self.serviceName == "arwen":
            while True:
                time.sleep( 5 )
                cmd = "ls " + self.pool[name]['socket'] + "/peptide_*.done | wc -l"
                ans = self.client.exec_command(cmd)
                o =  ans[1].read()
                #print o
                #print type(o)
                if len(self.pool[name]['inputs']) == int(o):
                    print "Exiting Remote"
                    break

        return name

    # We implement a generator b/c xml output of several blast cant be all loaded
    # We'll consume then in turn instead
    def pull(self, chunkid):

        print "-->" + self.app + "\n"

        if chunkid not in self.pool:
            raise ValueError("unknwon id \"" + chunkid + "\"")

        if self.app == "psipred":
            data = self._pullPsipred(chunkid)
            for d in data:
                pass
                #yield d

        elif self.app == "blastpgp":
            for i, e in enumerate(self.pool[chunkid]['inputs']):

                cmd = 'cat ' + self.pool[chunkid]['socket'] + "/peptide_" + str(i + 1) + ".blast"
                ans = self.client.exec_command(cmd)
                data = ans[1].read()
                #print ans[1].read()
                yield StringIO(data)



    def _pullPsipred (self, chunkid):
        data = []
        for i, e in enumerate(self.pool[chunkid]['inputs']):
            if self.serviceName == "migale" :
                ss2File = self.pool[chunkid]['socket'] + "/qsub_" + str(i + 1) + "/input.ss2"
            elif self.serviceName == "arwen" :
                ss2File = self.pool[chunkid]['socket'] + "/peptide_" + str(i + 1) + ".ss2"
            else :
                raise ValueError("unknown service name \"" + self.serviceName + "\"\n")
            cmd = "cat " + ss2File
            print cmd
            ans = self.client.exec_command(cmd)
            data.append(collection(stream=ans[1]))
        return data


    def _pushFiles(self, name, fileList): # sending fasta files

        cmd = ['scp'] + fileList + [ self.user + '@' + self.host + ':' + self.pool[name]['socket'] ]
        #print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).wait()
        return [ os.path.basename(f) for f in fileList ]

    def _wrap(self, name, peptidesList):
        for i, e in enumerate(peptidesList):
            e.fastaWrite(path=self.pool[name]['localCache'], fileName="peptide_" + str(i + 1))

        pyproteins.services.utils.make_tarfile(self.pool[name]['localCache'] + '.tar.gz', self.pool[name]['localCache'])
        self._pushFiles(name, [self.pool[name]['localCache'] + '.tar.gz'])

        cmd = 'cd ' + self.pool[name]['socket'] + ';tar -xzf ' + name + '.tar.gz; for ifile in  ' + name + '/peptide_*.fasta;do mv $ifile ./;done;rmdir ' + name
        stdin, stdout, stderr = self.client.exec_command(cmd)
        stdout.channel.recv_exit_status(); # blocking

    def _pushPeptides(self, name, peptidesList): # dumping peptide object content to remote fasta file

        if self.pool[name]['localCache']:
            self._wrap(name, peptidesList)
            print 'wrapped'

        fList = []
        cmd = ''

        for i, e in enumerate(peptidesList):
            fname = self.pool[name]['socket'] + "/peptide_" + str(i + 1) + ".fasta"
            fList.append(fname)

            if self.pool[name]['localCache']:
                continue

            cmd += "(echo \"" + e.fasta.replace("\n", "\" && echo \"") + "\") > " + fname  + "\n"# echo format should be moved to Peptide object
            #print cmd
      #      ans = self.client.exec_command(cmd)

            #print ans[1].read()

            stdin, stdout, stderr = self.client.exec_command(cmd)
            while True:
                if stdout.channel.exit_status_ready():
                    break
            #self.client.exec_command(cmd)

        return fList

class collection():
    def __init__(self, fileName=None, string=None, stream=None):

        if not fileName and not string and not stream:
            raise ValueError("No input provided")
        bufStream = None

        if fileName:
            bufStream = open (fileName, "r")
        elif string:
            bufStream = StringIO('foo')
        else :
            bufStream = stream

        self.data = [ element(line) for line in bufStream if isValidRecord(line) ]

        if pyproteins.services.utils.hasMethod(bufStream, "close"):
            bufStream.close()

    @property
    def horiz(self):
        return ''.join([e.state for e in self.data])
    @property
    def aaSeq(self):
        return ''.join([e.aa for e in self.data]).upper()

    @property
    def dict(self):
        return [ e.pos + " " + e.aa + " " + e.state + " " + e.cProb + " " + e.hProb + " " + e.eProb for e in self.data ]

class element():
    def __init__(self, line):
        lBuffer = line.split()
        self.pos = lBuffer[0]
        self.aa = lBuffer[1]
        self.state = lBuffer[2]
        self.cProb = lBuffer[3]
        self.hProb = lBuffer[4]
        self.eProb = lBuffer[5]
    def __repr__(self):
        return "\t".join([ str(val) for val in [self.pos, self.aa, self.state, self.cProb, self.hProb, self.eProb] ])
