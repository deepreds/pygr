#!/usr/bin/env python

import sys, os, string, glob
from pygr import seqdb, cnestedlist, metabase

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'dummy'
    sys.exit()

seqlist = ['anoCar1', 'bosTau4', 'calJac1', 'canFam2', 'cavPor3', 'choHof1', 'danRer5', 'dasNov2', 'dipOrd1', 'echTel1', 'equCab2', 'eriEur1', 'felCat3', 'fr2', 'galGal3', 'gasAcu1', 'gorGor1', 'hg18', 'loxAfr2', 'micMur1', 'mm9', 'monDom4', 'myoLuc1', 'ochPri2', 'ornAna1', 'oryCun1', 'oryLat2', 'otoGar1', 'panTro2', 'petMar1', 'ponAbe2', 'proCap1', 'pteVam1', 'rheMac2', 'rn4', 'sorAra1', 'speTri1', 'taeGut1', 'tarSyr1', 'tetNig1', 'tupBel1', 'turTru1', 'vicPac1', 'xenTro2']

pygrDir = '/data/server'
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

