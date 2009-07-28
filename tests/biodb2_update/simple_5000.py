#!/usr/bin/env python

# THIS SCRIPT RUNS BIODB2 PYGR XMLRPC SERVER: WRITTEN BY NAMSHIN KIM N@RNA.KR

import os, sys
from pygr import cnestedlist, seqdb, metabase

pygrDir = '/data/server'

mdb = metabase.MetabaseList(pygrDir)

for ix in mdb.dir('Bio'):
    if '.txt' in ix: continue
    mdb(ix)

server = metabase.ResourceServer(mdb, 'biodb2_5000', withIndex = True, \
    port = 5000, host = 'biodb2.bioinformatics.ucla.edu', \
    downloadDB = '/data/server/downloadable/.pygr_data')
server.serve_forever()

