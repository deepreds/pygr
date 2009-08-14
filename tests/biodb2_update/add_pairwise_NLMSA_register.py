#!/usr/bin/env python

import sys, os, string, glob

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'dummy'
    sys.exit()

worldbasepath = '/data/server'
downloadablepath = '/data/server/downloadable'
os.environ['WORLDBASEPATH'] = worldbasepath

from pygr import seqdb, metabase, cnestedlist
from pygr.downloader import *
from pygr.nlmsa_utils import NLMSABuilder
from step01_register_updated_seqdb import sprotDict, docstringdict

mdb1 = metabase.MetabaseList(worldbasepath)
mdb2 = metabase.MetabaseList(downloadablepath)

txtDir = '/data/axtNet_TEXT'
genomeList = [os.path.basename(ix).replace('.txt.gz', '') for ix in glob.glob(os.path.join(txtDir, '*.txt.gz'))]

genomeDir = '/data/PYGRDATA'
msaDir = '/data/NLMSA'
compressDir = '/data/axtNet_TEXT'
scpDir = 'biodb.bioinformatics.ucla.edu:/Volumes/Quadra1/apache/PYGRDATA/'
srcUrl = 'http://biodb.bioinformatics.ucla.edu/PYGRDATA'

registerDict = {}
for genoname in genomeList:
    registerDict[genoname] = 'Bio.MSA.UCSC.%s' % genoname

while 1:
    x = raw_input('If you finish copy: Press 1 Key to Exit')
    if x == '1': break

for genoname, pygrname in registerDict.items():
    genome1, genome2 = genoname.split('_pairwise')
    genome2 = genome2[0].lower() + genome2[1:]
    newurl = '%s/%s.txt.gz' % (srcUrl, genoname)
    dfile = SourceURL(newurl)
    dfile.__doc__ = genoname + '.txt in textfile dump format'
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname + '.txt', dfile)
    nbuilder = NLMSABuilder(dfile)
    nbuilder.__doc__ = '%s pairwise alignment with %s from UCSC genome browser' % (genome1, genome2)
    mdb2.add_resource('Bio.MSA.UCSC.' + genoname, nbuilder)
    mdb2.commit()

    msa = cnestedlist.NLMSA(os.path.join(msaDir, genoname), 'r')
    msa.__doc__ = '%s pairwise alignment with %s from UCSC genome browser' % (genome1, genome2)
    mdb1.add_resource('Bio.MSA.UCSC.' + genoname, msa)
    mdb1.commit()

