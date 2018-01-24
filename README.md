<h1>GPCR automodel suite : receptor modelling module</h1>
*Last updated GL - 03282016*
<h2>Installation</h2>
Create and move to your project directory
>$ mkdir GPCR
>
>$ cd GPCR

**All commands**, unless otherwise mentioned, **are intended to be executed from this folder**.

1. Set up a dedicated [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/).
We will use the python library already present on the system (the "--system-site-packages" flag), to avoid
the usually painfull installation of numpy.  
To create and activate the virtual environment simply:
>$ virtualenv --system-site-packages venv
>
>$ source venv/bin/activate`

2. Install the pyproteins library
`(venv)$ mkdir venv/modules
(venv)$ pip install pyproteins --target=venv/modules`

3. Install the [pathos library](https://pypi.python.org/pypi/pathos) for multi-threading
This library is not packaged, it is provided in *venv/modules/pyproteins/external.*
First, let's proceed to pathos depedencies installation:  
pip install following packages (eg: `pip install dill`)
    * dill, version >= 0.2.4
    * pox, version >= 0.2.2
    * ppft, version >= 1.6.4.5
    * multiprocess, version >= 0.70.3

	Then install the pathos library itself:
>(venv)$ tar -xjf venv/modules/pyproteins/external/pathos.tar.bz
>
>(venv)$ cd pathos-master
>
>(venv)$ python setup.py build
>
>(venv)$ python setup.py install` 

4. Modify the venv/bin/activate script to indicate python the pyproteins location

	>export PYTHONPATH="\$VIRTUAL_ENV/modules:$PYTHONPATH"
	>
	>export DRMAA_LIBRARY_PATH="/opt/sge/lib/lx24-amd64/libdrmaa.so"
	>
	>export HHLIB="/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh"

	Finally, restart the virtual environment to load all settings,

	>$ deactivate
	>
	>$ source venv/activate

*Note*: your working directories dont have to be in this current directory. But you **must** start the virtual environment to be able to use pyproteins.

<h2>Configuring</h2>
A template configuration file is provided in *venv/modules/pyproteins/conf/confModule1.json*

>{
>    "executable" : {
>       "blastpgp" : "/usr/local/genome/blast/bin/blastpgp"
>  },
> "envVariables" : {
>        "HHLIB" : "/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh",
>        "DRMAA_LIBRARY_PATH" : "/opt/sge/lib/lx24-amd64/libdrmaa.so"
>    },
>    "appendPath" : ["/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/bin"],
>    "pddDir" : ["/projet/extern/save/gulaunay/work/GPCR/data/MAKE_CORE/PDB_DERIVED"],
>    "pdbDir" : ["/projet/extern/save/gulaunay/work/GPCR/data/MAKE_CORE/PDB"],
>    "BLASTDBROOT" : "/projet/extern/save/gulaunay/work/GPCR/data/GPCRDB/",
>    "blastDb" : "GPCRs",
>    "registeredTemplate" : [
>        "4dajA",
>       "4djhA",
>        "4dklA",
>        "4ea3A",
>        "4ej4A",
>        "3emlA",
>        ...
>        "3vw7A"
>    ]
>}

The pddDir must point to the root folder containing your template metadata.
The pdbDir must point to the root folder containing your template pdb files.
The BLASTDBROOT must point to the root folder containing the blast database files (*phr, pin, psq* extensions).
The blastDb is the name of the blast database (eg: the file GPCRs.phr must exist and be in the *$BLASTDBROOT* folder).
<h2>Usage</h2>
The main script is provided in *venv/modules/pyproteins/bin/module1.py*. This file can be moved if needed.

<h4>Listing available templates</h4>
A crude listing of templates can be displayed by typing,   
`python module1.py -c confModule1.json --templateList`
By default, for all applications the entiere set of templates is loaded and use. A limited number of templates can though be specified with the  `--templateSelect` flag.
<h4>Creating a Query</h4>
Just run `$python module1.py --workDir myWorkDir -c confModule1.json -q or1g1.fasta --noTemplate --sge`
The `--sge` is **mandatory** and you have to be logged on to migale. The `--noTemplate` option will avoid the template loading time.  **The query sequence file must have a ".fasta" extension**.

It can take a while and shall produce a *beanQuery* file which stores the path to all files required to reload the same query later-on. The file content is self-explanatory with path to fasta, blast (xml output) files and psipred results folder.
<h4>(Re)building a template MSA</h4>
You can choose to modify the blast database in the configuration file and re-build the multiple sequence alignment of any template. I strongly advise to do it on your first run, anyway.
`python module1.py -c confModule1.json --templateMsaRebuild --sge --workDir myWorkDir`

<h4>Threading a query against a library of templates</h4>
`python module1.py -c confModule1.json -q queryBean.json --sge --workDir myWorkDir --templateSelect=1u19A,2rh1A  --hhThread`
<h4>Build homology-based structures</h3>

`python module1.py -c confModule1.json -q queryBean.json --sge --workDir myWorkDir --templateSelect=1u19A,2rh1A  --hhThread --nModel=10`
Setting the nModel parameter to any value greater than zero will trigger *modeller* downstream and ask it to build the specified number of models **for each template**.



<h3>Required Data</h3>
<h5>PDB structure file</h5>
* naming convention: pdb[*TEMPLATE_ID*].ent
* location: anywhere under  the configuration value of the *pdbDir* key

<h5>PDB metadata file</h5>
* naming convention: PDB_ID + chainID
* location: a folder named as stated above, anywhere under  the configuration value of the *pddDir* key
   
<h5>BLAST sequence databases</h5>
* content : a triplet of files with extensions  *phr, pin, psq* 
* naming convention: none, default is GPCRs
* location: in the BLASTDBROOT folder