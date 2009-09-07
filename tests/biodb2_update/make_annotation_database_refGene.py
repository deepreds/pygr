#!/usr/bin/env python

import ConfigParser, sys, os, string, gzip
from pygr.mapping import Collection
from pygr import seqdb, cnestedlist, metabase

#config = ConfigParser.ConfigParser({'testOutputBaseDir' : '.', 'smallSampleKey': ''})
#config.read([ os.path.join(os.path.expanduser('~'), '.pygrrc'), os.path.join(os.path.expanduser('~'), 'pygr.cfg'), '.pygrrc', 'pygr.cfg' ])
#msaDir = config.get('megatests_hg18', 'msaDir')
#seqDir = config.get('megatests_hg18', 'seqDir')
#smallSampleKey = config.get('megatests_hg18', 'smallSampleKey')
#testInputDB = config.get('megatests', 'testInputDB')
#testInputDir = config.get('megatests', 'testInputDir')
#testOutputBaseDir = config.get('megatests', 'testOutputBaseDir')
#
#if smallSampleKey:
#    smallSamplePostfix = '_' + smallSampleKey
#else:
#    smallSamplePostfix = ''

# DIRECTIONARY FOR DOC STRING OF SEQDB
genoName = 'hg18'
docString = 'Human Genome (May 2006)'

txtDir = os.path.abspath('ucsc')
msaDir = os.path.abspath('ucsc')

# exon -> intron -> exon schema
# annotation of exon and introns

os.environ['WORLDBASEPATH'] = os.path.abspath('ucsc')
from pygr import worldbase
hg18 = seqdb.SequenceFileDB('/data/GENOMES/hg18')
hg18.__doc__ = docString
worldbase.Bio.Seq.Genome.HUMAN.hg18 = hg18
hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()

exon_slices = Collection(filename=os.path.join(msaDir, 'refGene_exons.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
intron_slices = Collection(filename=os.path.join(msaDir, 'refGene_introns.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
start_slices = Collection(filename=os.path.join(msaDir, 'refGene_starts.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
stop_slices = Collection(filename=os.path.join(msaDir, 'refGene_stops.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
tu_slices = Collection(filename=os.path.join(msaDir, 'refGene_tus.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
cds_slices = Collection(filename=os.path.join(msaDir, 'refGene_cdss.cdb'), \
    intKeys=True, protocol=2, mode='cr', writeback=False)
exon_db = seqdb.AnnotationDB(exon_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, exon_id=6, exon_type=7, bin=8))
intron_db = seqdb.AnnotationDB(intron_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, intron_id=6, splice5=7,splice3=8,bin=9))
start_db = seqdb.AnnotationDB(start_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, is_splice=6, start_size=7, bin=8))
stop_db = seqdb.AnnotationDB(stop_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, is_splice=6, stop_size=7, bin=8))
tu_db = seqdb.AnnotationDB(tu_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, bin=6))
cds_db = seqdb.AnnotationDB(cds_slices, hg18, \
    sliceAttrDict=dict(id=0, start=1, stop=2, orientation=3, int_id=4, gbacc=5, bin=6))
exon_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_exons'), 'w', pairwiseMode=True, bidirectional=False)
intron_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_introns'), 'w', pairwiseMode=True, bidirectional=False)
start_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_starts'), 'w', pairwiseMode=True, bidirectional=False)
stop_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_stops'), 'w', pairwiseMode=True, bidirectional=False)
tu_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_tus'), 'w', pairwiseMode=True, bidirectional=False)
cds_msa = cnestedlist.NLMSA(os.path.join(msaDir, 'refGene_cdss'), 'w', pairwiseMode=True, bidirectional=False)

intExonId, intIntronId, intStartId, intStopId, intTuId, intCdsId = 0, 0, 0, 0, 0, 0

graphfile = open(os.path.join('refGene_splicegraph.txt'), 'w')

infile = gzip.GzipFile(os.path.join(txtDir, 'refGene.txt.gz'), 'r')
while 1:
    lines = infile.readline()
    if lines == '': break
    #print lines
    bin, qName, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds \
        = lines.split('\t')[:11]
    if strand == '+': orientation = 1
    elif strand == '-': orientation = -1
    else: continue
    bin, txStart, txEnd, cdsStart, cdsEnd = int(bin), int(txStart), int(txEnd), int(cdsStart), int(cdsEnd)
    exonEnds = exonEnds.strip()
    exStartList = [int(ix) for ix in exonStarts[:-1].split(',')]
    exEndList = [int(ix) for ix in exonEnds[:-1].split(',')]
    exons = [(ix, iy) for ix, iy in zip(exStartList, exEndList)]
    introns = [(ix, iy) for ix, iy in zip(exEndList[:-1], exStartList[1:])]
    exonPaths = []
    for iEx in range(len(exons)):
        iExon = iEx + 1
        exon_id = '%s.E%d' % (qName, iExon)
        iStart, iEnd = exons[iEx]
        if iEnd <= cdsStart: exon_type = 'UTR'
        elif iStart >= cdsEnd: exon_type = 'UTR'
        elif cdsStart <= iStart < iEnd <= cdsEnd: exon_type = 'CDS'
        else: exon_type = 'CDS/UTR'
        oneRow = (chrom, iStart, iEnd, orientation, intExonId, qName, exon_id, exon_type, bin)
        #print oneRow
        exon_slices[intExonId] = oneRow
        exon = exon_db[intExonId]
        exon_msa.addAnnotation(exon)
        exonPaths.append((intExonId, exon_id, qName, iStart, iEnd))
        intExonId += 1
    intronPaths = []
    for iInt in range(len(introns)):
        iIntron = iInt + 1
        intron_id = '%s.I%d' % (qName, iIntron)
        jStart, jEnd = introns[iInt]
        jSplice5, jSplice3 = str(hg18[chrom][jStart:jStart+2]).upper(), str(hg18[chrom][jEnd-2:jEnd]).upper()
        oneRow = (chrom, jStart, jEnd, orientation, intIntronId, qName, intron_id, jSplice5, jSplice3, bin)
        #print oneRow
        intron_slices[intIntronId] = oneRow
        intron = intron_db[intIntronId]
        intron_msa.addAnnotation(intron)
        intronPaths.append((intIntronId, intron_id, qName, jStart, jEnd))
        intIntronId += 1
    if len(intronPaths):
        for iInt in range(len(intronPaths)):
            exon1, exon2 = exonPaths[iInt], exonPaths[iInt+1]
            graphfile.write('\t'.join(map(str, exon1 + intronPaths[iInt] + exon2)) + '\n')
    if cdsStart < cdsEnd:
        startExons = [ix for ix in exons if ix[0] <= cdsStart <= ix[1]]
        if len(startExons) != 1:
            print >> sys.stderr, 'wrong exon', lines, startExons
        startExon = startExons[0]
        iStartExon = exons.index(startExon)
        if startExon[1] - cdsStart < 3: # spliced start codon
            startExon2 = exons[iStartExon+1][0], exons[iStartExon+1][0] + 3 - (startExon[1] - cdsStart)
        else:
            startExon2 = []
        stopExons = [ix for ix in exons if ix[0] <= cdsEnd <= ix[1]]
        if len(stopExons) != 1:
            print >> sys.stderr, 'wrong exon', lines, stopExons
        stopExon = stopExons[0]
        iStopExon = exons.index(stopExon)
        if cdsEnd - stopExon[0] < 3: # spliced stop codon
            stopExon2 = exons[iStopExon-1][1] - 3 + cdsEnd - stopExon[0], exons[iStopExon-1][1]
        else:
            stopExon2 = []
        if len(startExon2):
            start1, start2, start3, start4 = cdsStart, startExon[1], startExon2[0], startExon2[1]
            oneRow1 = (chrom, start1, start2, orientation, intStartId, qName, 1, start2-start1, bin)
            start_slices[intStartId] = oneRow1
            start = start_db[intStartId]
            start_msa.addAnnotation(start)
            intStartId += 1
            oneRow2 = (chrom, start3, start4, orientation, intStartId, qName, 1, start4-start3, bin)
            start_slices[intStartId] = oneRow2
            start = start_db[intStartId]
            start_msa.addAnnotation(start)
            intStartId += 1
        else:
            start1, start2, start3, start4 = cdsStart, cdsStart+3, 0, 0
            oneRow1 = (chrom, start1, start2, orientation, intStartId, qName, 0, start2-start1, bin)
            start_slices[intStartId] = oneRow1
            start = start_db[intStartId]
            start_msa.addAnnotation(start)
            intStartId += 1
        if len(stopExon2):
            stop1, stop2, stop3, stop4 = stopExon2[0], stopExon2[1], stopExon[0], cdsEnd
            oneRow1 = (chrom, stop3, stop4, orientation, intStopId, qName, 1, stop4-stop3, bin)
            stop_slices[intStopId] = oneRow1
            stop = stop_db[intStopId]
            stop_msa.addAnnotation(stop)
            intStopId += 1
            oneRow2 = (chrom, stop1, stop2, orientation, intStopId, qName, 1, stop2-stop1, bin)
            stop_slices[intStopId] = oneRow2
            stop = stop_db[intStopId]
            stop_msa.addAnnotation(stop)
            intStopId += 1
        else:
            stop1, stop2, stop3, stop4 = 0, 0, cdsEnd-3, cdsEnd
            oneRow1 = (chrom, stop3, stop4, orientation, intStopId, qName, 0, stop4-stop3, bin)
            stop_slices[intStopId] = oneRow1
            stop = stop_db[intStopId]
            stop_msa.addAnnotation(stop)
            intStopId += 1
        #print startExon, startExon2, start1, start2, start3, start4
        #print stopExon, stopExon2, stop1, stop2, stop3, stop4
        #if len(startExon2) or len(stopExon2): raw_input('Press Start/Stop Key')

    oneRow = (chrom, txStart, txEnd, orientation, intTuId, qName, bin)
    #print oneRow
    tu_slices[intTuId] = oneRow
    tu = tu_db[intTuId]
    tu_msa.addAnnotation(tu)
    intTuId += 1
    if cdsStart < cdsEnd:
        oneRow = (chrom, cdsStart, cdsEnd, orientation, intCdsId, qName, bin)
        #print oneRow
        cds_slices[intCdsId] = oneRow
        cds = cds_db[intCdsId]
        cds_msa.addAnnotation(cds)
        intCdsId += 1

    #raw_input('Press Enter Key')

infile.close()
graphfile.close()

exon_db.clear_cache()
intron_db.clear_cache()
start_db.clear_cache()
stop_db.clear_cache()
tu_db.clear_cache()
cds_db.clear_cache()

exon_slices.close()
intron_slices.close()
start_slices.close()
stop_slices.close()
tu_slices.close()
cds_slices.close()

exon_msa.build()
intron_msa.build()
start_msa.build()
stop_msa.build()
tu_msa.build()
cds_msa.build()

exon_msa.save_seq_dict()
intron_msa.save_seq_dict()
start_msa.save_seq_dict()
stop_msa.save_seq_dict()
tu_msa.save_seq_dict()
cds_msa.save_seq_dict()

# schema binding: Bio.Annotation.UCSC.refGene.HUMAN.hg18.exons()

mdb = worldbase._mdb

exon_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC exons from refGene'
intron_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC introns from refGene'
start_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC start codons from refGene'
stop_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC stop codons from refGene'
tu_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC transcription units from refGene'
cds_db.__doc__ = 'AnnotationDB for hg18 referenced UCSC coding units from refGene'
            
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.exons', exon_db)
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.introns', intron_db)
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.starts', start_db)
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.stops', stop_db)
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.tus', tu_db)
mdb.add_resource('Bio.Annotation.UCSC.refGene.HUMAN.hg18.cdss', cds_db)

exon_msa.__doc__ = 'NLMSA for hg18 referenced UCSC exons from refGene'
intron_msa.__doc__ = 'NLMSA for hg18 referenced UCSC introns from refGene'
start_msa.__doc__ = 'NLMSA for hg18 referenced UCSC start codons from refGene'
stop_msa.__doc__ = 'NLMSA for hg18 referenced UCSC stop codons from refGene'
tu_msa.__doc__ = 'NLMSA for hg18 referenced UCSC transcription units from refGene'
cds_msa.__doc__ = 'NLMSA for hg18 referenced UCSC coding units from refGene'

mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.exons', exon_msa)
mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.introns', intron_msa)
mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.starts', start_msa)
mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.stops', stop_msa)
mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.tus', tu_msa)
mdb.add_resource('Bio.MSA.UCSC.refGene.HUMAN.hg18.cdss', cds_msa)

exon_schema = metabase.ManyToManyRelation(hg18, exon_db, bindAttrs = ('exon', ))
intron_schema = metabase.ManyToManyRelation(hg18, intron_db, bindAttrs = ('intron', ))
start_schema = metabase.ManyToManyRelation(hg18, start_db, bindAttrs = ('start_codon', ))
stop_schema = metabase.ManyToManyRelation(hg18, stop_db, bindAttrs = ('stop_codon', ))
tu_schema = metabase.ManyToManyRelation(hg18, tu_db, bindAttrs = ('tu', ))
cds_schema = metabase.ManyToManyRelation(hg18, cds_db, bindAttrs = ('cds', ))

exon_schema.__doc__ = 'Schema for hg18 referenced UCSC exons from refGene'
intron_schema.__doc__ = 'Schema for hg18 referenced UCSC introns from refGene'
start_schema.__doc__ = 'Schema for hg18 referenced UCSC start codons from refGene'
stop_schema.__doc__ = 'Schema for hg18 referenced UCSC stop codons from refGene'
tu_schema.__doc__ = 'Schema for hg18 referenced UCSC transcription units from refGene'
cds_schema.__doc__ = 'Schema for hg18 referenced UCSC coding units from refGene'

mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.exons', exon_schema)
mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.introns', intron_schema)
mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.starts', start_schema)
mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.stops', stop_schema)
mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.tus', tu_schema)
mdb.add_schema('Bio.Annotation.UCSC.refGene.HUMAN.hg18.cdss', cds_schema)

mdb.commit()



"""

"""

