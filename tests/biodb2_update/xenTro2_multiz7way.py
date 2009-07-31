#!/usr/bin/env python

# SCRIPT TO BUILD MAF NLMSA AND TEXT BINARIES: WRITTEN BY NAMSHIN KIM N@RNA.KR

import sys, os, string, glob
from pygr import seqdb, cnestedlist, metabase

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'dummy'
    sys.exit()

# SEQUENCE LIST SHOULD BE UPDATED. ALSO, THE NAME OF THE SCRIPT SHOULD BE COINCIDE WITH NLMSA NAME.
# E.G. hg18_multiz44way
seqlist = ['xenTro2', 'galGal2', 'monDom4', 'hg18', 'mm8', 'rn4', 'danRer4']

pygrDir = '/data/server'
os.environ['WORLDBASEPATH'] = pygrDir
dnDb = '/data/server/downloadable/.pygr_data'
seqdir = '/data/GENOMES'
msadir = '/data/NLMSA'
mafdir = os.path.join('/data/MAF', os.path.basename(args[0])[:-3])
txtdir = '/data/MAF_TEXT'

mdb = metabase.MetabaseList(pygrDir)
genomeList = [ix for ix in mdb.dir('Bio.Seq.Genome') if ix[-6:] != '.fasta' and ix[-4:] != '.txt']
genoDict = {}
for pygrstr in genomeList:
    j1, genoname = pygrstr.rsplit('.', 1)
    if genoname not in seqlist: continue
    genoDict[genoname] = pygrstr
genomes = {}
for orgstr in seqlist:
    genomes[orgstr] = mdb(genoDict[orgstr])

genomeUnion = seqdb.PrefixUnionDict(genomes)

maflist = glob.glob(os.path.join(mafdir, '*.maf'))
maflist.sort()

os.system('rm -f ' + os.path.join(msadir, os.path.basename(mafdir)) + '*')

msa = cnestedlist.NLMSA(os.path.join(msadir, os.path.basename(mafdir)), 'w', genomeUnion, maflist)
msa.save_seq_dict()

myorg, mymultiz = string.split(os.path.basename(args[0])[:-3], '_')
msa.__doc__ = myorg + ' referenced ' + mymultiz + ' alignments (' + ', '.join(seqlist) + ')'
mdb.add_resource('Bio.MSA.UCSC.' + os.path.basename(mafdir), msa)
mdb.commit()

cnestedlist.dump_textfile(os.path.join(msadir, os.path.basename(mafdir)), \
    os.path.join(txtdir, os.path.basename(mafdir)) + '.txt')

