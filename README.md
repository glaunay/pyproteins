

#


## Development purpose

```sh
git clone https://github.com/glaunay/pyproteins.git  LOCAL_PATH_TO_REPO
```

Check for dependencies in `LOCAL_PATH_TO_REPO/setup.py`

```python
sys.path.append(LOCAL_PATH_TO_REPO/src)
```

## Deployment module

```sh
pip install pyproteins
```

## Update pip module
 TO DO

### modules organization

#### container modules

##### Core
Frequently used container and functions.
Downloading and parsing XML in a generic way.
Two-level dictionnaries.

##### customCollection

Managment of a collection of object which can be constructed from flat files (ex: uniprot).

* collectionPath: path to a directory of flat files
* constructor: object constructor
* indexer function to map filename and collection key

#### sequence modules

##### msa

* Load alignemnt files
* Filter alignment
* Compute statistics

##### peptide
Fasta framed informations, like amino-acid types, domain boundaries and secondary structure states.

##### psipred
DEPRECATED

#### sequence alignment Modules

* nw_custom: the Needleman & Wunsch implementaion
* scoringFunctions : AA, AA/SS2, PWA
* scoringMatrix : Blosum62

#### utils modules

test Cecile
