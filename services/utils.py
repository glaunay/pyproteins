import sys
import os
import errno
import StringIO

import os.path
from types import ModuleType

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def toInStream(inputData):
    if hasattr(inputData, 'read'):
        return inputData
    if os.path.isfile(inputData):
        f = open(inputData, 'r')
        return f

    return StringIO.StringIO(inputData)

def hasMethod(obj, askedMethod):
    l = [method for method in dir(obj) if callable(getattr(obj, method))]
    if askedMethod in l:
        return True
    return False

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise





def rreload(module):
    """Recursively reload modules."""
    reload(module)
    for attribute_name in dir(module):
        attribute = getattr(module, attribute_name)
        if type(attribute) is ModuleType:
            rreload(attribute)

'''
    parse a tsv file into a list of dictionary, dictionary keys are extracted from the first line (aka column headers)
'''

def tsvToDictList(fileName):
    buffer = tabularFileToList(fileName, separator = "\t")
    print len(buffer)
    keymap = buffer.pop(0)
    print len(buffer)
    data = []
    for d in buffer:
        data.append({})
        for i,x in enumerate(d):
            data[-1][keymap[i]] = x

    return {'keymap' : keymap, 'data' : data}

def tabularFileToList(fileName, separator = ","):

    with open (fileName, "r") as f:
        data = [line.strip('\n').split(separator) for line in f]

    return data