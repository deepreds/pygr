#!/usr/bin/env python

# SCRIPT TO BUILD MAF NLMSA AND TEXT BINARIES: WRITTEN BY NAMSHIN KIM N@RNA.KR

import sys, os, string, glob
from pygr import seqdb, cnestedlist, metabase

args = sys.argv
if len(args) != 2:
    print 'Usage:', args[0], 'inputDir'
    sys.exit()

# SEQUENCE LIST SHOULD BE UPDATED. ALSO, THE NAME OF THE SCRIPT SHOULD BE COINCIDE WITH NLMSA NAME.
inputDir = args[1]
if inputDir[-1] == os.sep:
    inputDir = inputDir[:-len(os.sep)]
pathList = inputDir.split(os.sep)
seqdb1 = pathList[-2]
seqdb2 = pathList[-1][2].lower() + pathList[-1][3:]
if seqdb2 == 'self':
    seqlist = [seqdb1]
else:
    seqlist = [seqdb1, seqdb2]

pygrDir = '/data/server'
os.environ['WORLDBASEPATH'] = pygrDir
dnDb = '/data/server/downloadable/.pygr_data'
dnDir = '/data/server/downloadable'
seqdir = '/data/GENOMES'
msadir = '/data/NLMSA'
txtdir = '/data/axtNet_TEXT'

mdb = metabase.MetabaseList(pygrDir)
genomeList = [ix for ix in mdb.dir('Bio.Seq.Genome') if ix[-6:] != '.fasta' and ix[-4:] != '.txt']
genoDict = {}
for pygrstr in genomeList:
    j1, genoname = pygrstr.rsplit('.', 1)
    if genoname not in seqlist: continue
    genoDict[genoname] = pygrstr

genomes = {}
for orgstr in seqlist:
    try:
        genomes[orgstr] = mdb(genoDict[orgstr])
    except:
        print 'GENOME ASSEMBLY NOT REGISTERED: %s' % orgstr
        pass
if len(genomes) not in (1, 2): sys.exit()

genomeUnion = seqdb.PrefixUnionDict(genomes)

os.system('gzip -d %s' % os.path.join(inputDir, 'axtNet/*.net.axt.gz'))
axtList = glob.glob(os.path.join(inputDir, 'axtNet/*.net.axt'))
if len(axtList) == 0:
    os.system('gzip -d %s' % os.path.join(inputDir, '*.net.axt.gz'))
    axtList = glob.glob(os.path.join(inputDir, '*.net.axt'))

if len(axtList) == 0: sys.exit('AXTNET FILE NOT FOUND: %s' % args[1])

if len(seqlist) == 1:
    msaname = '%s_pairwise%s' % (seqdb1, 'Self')
else:
    msaname = '%s_pairwise%s' % (seqdb1, seqdb2[0].upper() + seqdb2[1:])

os.system('rm -f ' + os.path.join(msadir, msaname) + '*')

msa = cnestedlist.NLMSA(os.path.join(msadir, msaname), 'w', genomeUnion, axtFiles = axtList)
msa.save_seq_dict()

#myorg, mymultiz = msaname.split(msaname, '_')
#msa.__doc__ = myorg + ' referenced ' + mymultiz + ' alignments (' + ', '.join(seqlist) + ')'
#mdb.add_resource('Bio.MSA.UCSC.' + msaname, msa)
#mdb.commit()

cnestedlist.dump_textfile(os.path.join(msadir, msaname), \
    os.path.join(txtdir, msaname) + '.txt')


"""
[deepreds@mbi136-219 axtNet]$ pwd
/data/axtNet/hg18/vsSelf/axtNet

[deepreds@mbi136-219 biodb2_update]$ l /data/axtNet/hg18
total 16
drwxr-xr-x 3 deepreds users  111 Aug  2 05:29 vsAnoCar1
drwxr-xr-x 3 deepreds users  111 Aug  2 05:29 vsBosTau2
drwxr-xr-x 3 deepreds users 4096 Aug  2 05:34 vsBosTau3
drwxr-xr-x 3 deepreds users  111 Aug  2 05:38 vsBosTau4
drwxr-xr-x 4 deepreds users  132 Aug  2 05:40 vsBraFlo1
"""

