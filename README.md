<h1>GPCR automodel suite receptor modelling module</h1>

<h2>Installation</h2>
Create and move to your project directory
`$ mkdir GPCR
$ cd GPCR
`
**All commands**, unless otherwise mentioned, **are intended to be executed from this folder**.

1. Setting up a dedicated [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/)
We will use the python library already present on the system (the "--system-site-packages" flag), to avoid
the usually painfull installation of numpy. To create and activate the virtual environment simply:
`$ virtualenv --system-site-packages
$ source venv/activate`

2. Installing the pyproteins librairy
`(venv)$ mkdir venv/modules
(venv)$ pip install pyproteins --target=venv/modules`

3. Installing the [pathos library](https://pypi.python.org/pypi/pathos) for multi-threading

This library is not packaged, it is provided in venv/modules/pyproteins/external.
First, let's proceed to pathos depedencies installation:

pip install following packages (eg: `pip install dill`)

⋅⋅* dill, version >= 0.2.4
⋅⋅* pox, version >= 0.2.2
⋅⋅* ppft, version >= 1.6.4.5
⋅⋅* multiprocess, version >= 0.70.3

then install the pathos librairy itself
`(venv)$ tar -xjf venv/modules/pyproteins/external/pathos.tar.bz
(venv)$ cd pathos-master
(venv)$ python setup.py build
(venv)$ python setup.py install`

3. Modify the venv/bin/activate script to indicate python the pyproteins location

export PYTHONPATH="$VIRTUAL_ENV/modules:$PYTHONPATH"
export DRMAA_LIBRARY_PATH="/opt/sge/lib/lx24-amd64/libdrmaa.so"
export HHLIB="/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh"

Restart the virtual environment to load all settings,
`$ deactivate
$ source venv/activate`

Your working directories dont have to be in this current directory. But you **must** start the virtual environment to be able to use pyproteins.

<h2>Configuring</h2>
A

<h2>Usage</h2>

<h3>Creating a Query</h3>

<h3>(Re)building a template MSA</h3>

<h3>Creating a Query</h3>

<h3>Threading a query against a library of templates</h3>

<h3>Build homology-based structures</h3>

`python ../module1.py -c ../confModule1.json -q or1g1_queryBean.json --templateMsaRebuild --sge --workDir dev0321`


`python ../module1.py -c ../confModule1.json -q or1g1_queryBean.json --hhThread --sge --workDir dev0321a --templateSelect=1u19A,2rh1A`

