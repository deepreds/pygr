#!/usr/bin/env python

import sys, os, string, glob
from pygr import seqdb, cnestedlist

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'dummy'
    sys.exit()

seqlist = ['anoCar1', 'bosTau4', 'calJac1', 'canFam2', 'cavPor3', 'choHof1', 'danRer5', 'dasNov2', 'dipOrd1', 'echTel1', 'equCab2', 'eriEur1', 'felCat3', 'fr2', 'galGal3', 'gasAcu1', 'gorGor1', 'hg18', 'loxAfr2', 'micMur1', 'mm9', 'monDom4', 'myoLuc1', 'ochPri2', 'ornAna1', 'oryCun1', 'oryLat2', 'otoGar1', 'panTro2', 'petMar1', 'ponAbe2', 'proCap1', 'pteVam1', 'rheMac2', 'rn4', 'sorAra1', 'speTri1', 'taeGut1', 'tarSyr1', 'tetNig1', 'tupBel1', 'turTru1', 'vicPac1', 'xenTro2']

os.environ['PYGRDATAPATH'] = 'mysql:PYGRDB_JAN06.pygr_create'
import pygr.Data

seqdir = '/result/pygr_data'
msadir = '/result/pygr_nlmsa'
mafdir = os.path.join('/download/MAF_UCSC', os.path.basename(args[0])[:-3])
txtdir = '/download/NLMSA_Pygr'

genomes = {}
genoDict = {}
for lines in open('2_register_seqdb.txt', 'r').xreadlines():
    genoname, pygrstr = lines.splitlines()[0].split('\t')
    genoDict[genoname] = pygrstr

for orgstr in seqlist:
    genomes[orgstr] = pygr.Data.getResource(genoDict[orgstr]) #seqdb.BlastDB(orgstr)

genomeUnion=seqdb.PrefixUnionDict(genomes)

maflist = glob.glob(os.path.join(mafdir, '*.maf'))
maflist.sort()

os.system('rm -f ' + os.path.join(msadir, os.path.basename(mafdir)) + '*')

msa = cnestedlist.NLMSA(os.path.join(msadir, os.path.basename(mafdir)), 'w', genomeUnion, maflist)
msa.save_seq_dict()

myorg, mymultiz = string.split(os.path.basename(args[0])[:-3], '_')
msa.__doc__ = myorg + ' referenced ' + mymultiz + ' alignments (' + ', '.join(seqlist) + ')'
for genoname, genome in msa.seqDict.prefixDict.items():
    pygr.Data.getResource.addResource(genoDict[genoname], genome)
pygr.Data.getResource.addResource('Bio.MSA.UCSC.' + os.path.basename(mafdir), msa)
pygr.Data.save()

cnestedlist.dump_textfile(os.path.join(msadir, os.path.basename(mafdir)), \
    os.path.join(txtdir, os.path.basename(mafdir)) + '.txt')

