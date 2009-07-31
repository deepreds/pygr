#!/usr/bin/env python

import sys, os, string

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'dummy'
    sys.exit()

worldbasepath = '/data/server'
downloadablepath = '/data/server/downloadable'
os.environ['WORLDBASEPATH'] = worldbasepath

from pygr import seqdb, metabase
from pygr.downloader import *
from pygr.nlmsa_utils import NLMSABuilder
from step01_register_updated_seqdb import sprotDict, docstringdict

mdb1 = metabase.MetabaseList(worldbasepath)
mdb2 = metabase.MetabaseList(downloadablepath)
genomeList = [ix for ix in mdb1.dir('Bio.MSA.UCSC') if ix[-6:] != '.fasta' and ix[-4:] != '.txt']
txtList = [ix for ix in mdb1.dir('Bio.MSA.UCSC') if ix[-4:] == '.txt']
genomeList.sort()
txtList.sort()

genomeDir = '/data/PYGRDATA'
compressDir = '/data/MAF_TEXT'
scpDir = 'biodb.bioinformatics.ucla.edu:/Volumes/Quadra1/apache/PYGRDATA/'
srcUrl = 'http://biodb.bioinformatics.ucla.edu/PYGRDATA'

registerDict = {}
for pygrname in genomeList:
    if pygrname + '.txt' not in txtList:
        genoname = pygrname.rsplit('.', 1)[-1]
        registerDict[genoname] = pygrname

while 1:
    x = raw_input('If you finish copy: Press 1 Key to Exit')
    if x == '1': break

for genoname, pygrname in registerDict.items():
    newurl = '%s/%s.txt.gz' % (srcUrl, genoname)
    dfile = SourceURL(newurl)
    dfile.__doc__ = genoname + '.txt in textfile dump format'
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname + '.txt', dfile)
    nbuilder = NLMSABuilder(dfile)
    nbuilder.__doc__ = genoname + ' multigenome alignment from UCSC genome browser'
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname, nbuilder)
    mdb2.commit()

