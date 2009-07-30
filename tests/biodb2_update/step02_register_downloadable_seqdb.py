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
from step01_register_updated_seqdb import sprotDict, docstringdict

mdb1 = metabase.MetabaseList(worldbasepath)
mdb2 = metabase.MetabaseList(downloadablepath)
genomeList = [ix for ix in mdb1.dir('Bio.Seq.Genome') if ix[-6:] != '.fasta' and ix[-4:] != '.txt']
fastaList = [ix for ix in mdb2.dir('Bio.Seq.Genome') if ix[-6:] == '.fasta']
genomeList.sort()
fastaList.sort()

genomeDir = '/data/GENOMES'
compressDir = '/data/GENOMES_GZIPPED'
scpDir = 'biodb.bioinformatics.ucla.edu:/Volumes/Quadra1/apache/GENOMES/'
srcUrl = 'http://biodb.bioinformatics.ucla.edu/GENOMES'

registerDict = {}
for pygrname in genomeList:
    if pygrname + '.fasta' not in fastaList:
        genoname = pygrname.rsplit('.', 1)[-1]
        registerDict[genoname] = pygrname
for genoname, pygrname in registerDict.items():
    if os.path.exists(os.path.join(compressDir, genoname)):
        if os.path.exists(os.path.join(compressDir, genoname, genoname + '.gz')):
            continue
    try:
        os.mkdir(os.path.join(compressDir, genoname))
        print 'DIRECTORY CREATED:', os.path.join(compressDir, genoname)
    except:
        pass
    q = 'cat %s/%s|gzip 1> %s/%s/%s.gz' % (genomeDir, genoname, compressDir, genoname, genoname)
    print q
    #os.system(q)

q = 'scp -r %s/* %s' % (compressDir, scpDir)
print q
#os.system(q)

while 1:
    x = raw_input('Copy Finished: Press 1 Key')
    if x == '1': break

for genoname, pygrname in registerDict.items():
    filename = genoname + '.gz'
    newurl = '%s/%s/%s' % (srcUrl, genoname, filename)
    if not docstringdict.has_key(genoname) or not sprotDict.has_key(genoname):
        print geoname, 'NOT REGISTERED: SKIPPED'
        continue

    if '.tar.gz' in filename: mytype = '.tar.gz'
    elif '.gz' in filename: mytype = '.gz'
    elif '.tgz' in filename: mytype = '.tar.gz'
    elif '.zip' in filename: mytype = '.zip'
    else:
        print >> sys.stderr, lines[:-1]
        continue

    if mytype != '.gz':
        src = SourceURL(newurl, filename = genoname + mytype, singleFile = True)
    else:
        src = SourceURL(newurl, filename = genoname + mytype)
    src.__doc__ = docstringdict[genoname] + ' FASTA File'
    mdb2.add_resource('Bio.Seq.Genome.' + sprotDict[genoname] \
        + '.' + genoname + '.fasta', src)
    rsrc = GenericBuilder('BlastDB', src)
    rsrc.__doc__ = docstringdict[genoname]
    mdb2.add_resource('Bio.Seq.Genome.' + sprotDict[genoname] \
        + '.' + genoname, rsrc)
    mdb2.commit()

