<h1>GPCR automodel suite receptor modelling module</h1>

<h2>Installation</h2>

1. Setting up a dedicated [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/)

In you environment folder create a `venv/modules` folder


create a

edit venv/bin/activate
_OLD_VIRTUAL_PATH="$PATH"
PATH="$VIRTUAL_ENV/bin:$PATH"
export PATH

export PYTHONPATH="$VIRTUAL_ENV/modules:$PYTHONPATH"
export DRMAA_LIBRARY_PATH="/opt/sge/lib/lx24-amd64/libdrmaa.so"
export HHLIB="/usr/local/genome/src/hhsuite-2.0.16-linux-x86_64/lib/hh"





3. Installing the pathos library (multi-threading)

    pip install following packages

⋅⋅* dill, version >= 0.2.4
⋅⋅* pox, version >= 0.2.2
⋅⋅* ppft, version >= 1.6.4.5
⋅⋅* multiprocess, version >= 0.70.3

    then
    `$ cd pathos-master
    $ python setup.py build
    $ python setup.py install`


<h2>Configuring</h2>

<h2>Usage</h2>

<h3>Creating a Query</h3>

<h3>(Re)building a template MSA</h3>

<h3>Creating a Query</h3>

<h3>Threading a query against a library of templates</h3>

<h3>Build homology-based structures</h3>

`python ../module1.py -c ../confModule1.json -q or1g1_queryBean.json --templateMsaRebuild --sge --workDir dev0321`


`python ../module1.py -c ../confModule1.json -q or1g1_queryBean.json --hhThread --sge --workDir dev0321a --templateSelect=1u19A,2rh1A`

