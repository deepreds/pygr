#!/usr/bin/env python

import sys, os, string, glob

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

txtDir = '/data/MAF_TEXT'
genomeList = [os.path.basename(ix).replace('.txt.gz', '') for ix in glob.glob(os.path.join(txtDir, '*.txt.gz'))]

genomeDir = '/data/PYGRDATA'
compressDir = '/data/MAF_TEXT'
scpDir = 'biodb.bioinformatics.ucla.edu:/Volumes/Quadra1/apache/PYGRDATA/'
srcUrl = 'http://biodb.bioinformatics.ucla.edu/PYGRDATA'

registerDict = {}
for genoname in genomeList:
    registerDict[genoname] = 'Bio.MSA.UCSC.%s' % genoname

while 1:
    x = raw_input('If you finish copy: Press 1 Key to Exit')
    if x == '1': break

for genoname, pygrname in registerDict.items():
    newurl = '%s/%s.txt.gz' % (srcUrl, genoname)
    dfile = SourceURL(newurl)
    dfile.__doc__ = genoname + '.txt in textfile dump format'
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname + '.txt', dfile)
    print 'Bio.MSA.UCSC.' + genoname + '.txt'
    nbuilder = NLMSABuilder(dfile)
    nbuilder.__doc__ = genoname + ' multigenome alignment from UCSC genome browser'
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname, nbuilder)
    print 'Bio.MSA.UCSC.' + genoname
    mdb2.commit()

