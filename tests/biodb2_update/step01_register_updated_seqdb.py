#!/usr/bin/env python

import sys, os, string, glob
from pygr import seqdb, cnestedlist, metabase

# A SCRIPT TO UPDATE WORLDBASE SEQDB RESOURCES. ALL UNIPROT MNEMONIC AND DOCSTRING ARE UPDATED MANUALLY
# WRITTEN BY NAMSHIN KIM: N@RNA.KR

# UNIPROT MNEMONIC KEYWORDS FOR EACH SPECIES: MNEMONIC FROM UNIPROT TAXONOMY BROWSER
sprotDict = {
    'AaegL1': 'AEDAE', 'aedes1': 'AEDAE', 'amel2': 'APIME', 'anoCar1': 'ANOCA', 'anoGam1': 'ANOGA', 'apiMel2': 'APIME',
    'apiMel3': 'APIME', 'bomMor0': 'BOMMO', 'bosTau2': 'BOVIN', 'bosTau3': 'BOVIN', 'canFam2': 'CANFA',
    'cavPor2': 'CAVPO', 'ce2': 'CAEEL', 'ci2': 'CIOIN', 'CSAV2': 'CIOSA', 'danRer1': 'DANRE', 'danRer2': 'DANRE',
    'danRer3': 'DANRE', 'danRer4': 'DANRE', 'dasNov1': 'DASNO', 'dm2': 'DROME', 'dp3': 'DROPS', 'dp4': 'DROPS',
    'droAna1': 'DROAN', 'droAna2': 'DROAN', 'droAna3': 'DROAN', 'droEre1': 'DROER', 'droEre2': 'DROER',
    'droGri1': 'DROGR', 'droGri2': 'DROGR', 'droMoj1': 'DROMO', 'droMoj2': 'DROMO', 'droMoj3': 'DROMO',
    'droPer1': 'DROPE', 'droSec1': 'DROSE', 'droSim1': 'DROSI', 'droVir1': 'DROVI', 'droVir2': 'DROVI',
    'droVir3': 'DROVI', 'droWil1': 'DROWI', 'droYak1': 'DROYA', 'droYak2': 'DROYA', 'echTel1': 'ECHTE',
    'equCab1': 'HORSE', 'eriEur1': 'ERIEU', 'felCat3': 'FELCA', 'fr1': 'FUGRU', 'fr2': 'FUGRU', 'galGal2': 'CHICK',
    'galGal3': 'CHICK', 'gasAcu1': 'GASAC', 'hg17': 'HUMAN', 'hg18': 'HUMAN', 'loxAfr1': 'LOXAF', 'mm5': 'MOUSE',
    'mm6': 'MOUSE', 'mm7': 'MOUSE', 'mm8': 'MOUSE', 'monDom1': 'MONDO', 'monDom2': 'MONDO', 'monDom4': 'MONDO',
    'ornAna1': 'ORNAN', 'oryCun1': 'RABIT', 'oryLat1': 'ORYLA', 'otoGar1': 'OTOGA', 'panTro1': 'PANTR',
    'panTro2': 'PANTR', 'rheMac1': 'MACMU', 'rheMac2': 'MACMU', 'rn3': 'RAT', 'rn4': 'RAT', 'sacCer1': 'YEAST',
    'sorAra1': 'SORAR', 'Sscrofa2': 'PIG', 'strPur1': 'STRPU', 'tetNig1': 'TETNG',
    'triCas2': 'TRICA', 'tupBel1': 'TUPGB', 'xenTro1': 'XENTR', 'xenTro2': 'XENTR', 'priPac1': 'PRIPA',
    'mm9': 'MOUSE', 'dm3': 'DROME', 'ce4': 'CAEEL', 'cb3': 'CAEBR', 'caeRem2': 'CAERE',
    'danRer5': 'DANRE', 'bosTau4': 'BOVIN', 'calJac1': 'CALJA', 'dasNov2': 'DASNO', 'dipOrd1': 'DIPOR',
    'ponAbe2': 'PONAB', 'vicPac1': 'LAMPA', 'choHof1': 'CHOHO', 'micMur1': 'MICMU', 'speTri1': 'SPETR',
    'loxAfr2': 'LOXAF', 'turTru1': 'TURTR', 'ochPri2': 'OCHPR', 'gorGor1': 'GORGO', 'proCap1': 'PROCA',
    'tarSyr1': 'TARSY', 'pteVam1': 'PTEVA', 'taeGut1': 'TAEGU', 'myoLuc1': 'MYOLU', 'equCab2': 'HORSE',
    'oryLat2': 'ORYLA', 'cavPor3': 'CAVPO', 'petMar1': 'PETMA', 'braFlo1': 'BRAFL', 'caeJap1': 'CAEJA',
    'caePb1': 'CAEPB', 'caePb2': 'CAEPB', 'caeRem3': 'CAERE', 'ce6': 'CAEEL'
    }

# NAME OF THE GENOME AND ASSEMBLY DATE FROM UCSC GENOME BROWSER, MULTIGENOME ALIGNMENT INFORMATION
docstringdict = {
    'anoCar1':'Lizard Genome (January 2007)',
    'anoGam1':'A. gambiae Genome (February 2003)',
    'apiMel2':'A. mellifera Genome (January 2005)',
    'apiMel3':'A. mellifera Genome (May 2005)',
    'bosTau2':'Cow Genome (March 2005)',
    'bosTau3':'Cow Genome (August 2006)',
    'bosTau4':'Cow Genome (October 2007)',
    'braFlo1':'Lancelet genome (March 2006)',
    'caeJap1':'C. japonica genome (March 2008)',
    'caePb1':'C. brenneri Genome (July 2007)',
    'caePb2':'C. brenneri Genome (February 2008)',
    'caeRem2':'C. remanei Genome (March 2006)',
    'caeRem3':'C. remanei Genome (May 2007)',
    'calJac1':'Marmoset Genome (June 2007)',
    'canFam2':'Dog Genome (May 2005)',
    'cavPor2':'Guinea Pig (October 2005)',
    'cavPor3':'Guinea Pig (February 2008)',
    'cb3':'C. briggsae Genome (January 2007)',
    'ce2':'C. elegans Genome (March 2004)',
    'ce4':'C. elegans Genome (January 2007)',
    'ce6':'C. elegans Genome (May 2008)',
    'choHof1':'Sloth Genome (July 2008)',
    'ci2': 'C. intestinalis Genome (March 2005)',
    'danRer1':'Zebrafish Genome (November 2003)',
    'danRer2':'Zebrafish Genome (Junuary 2004)',
    'danRer3':'Zebrafish Genome (May 2005)',
    'danRer4':'Zebrafish Genome (March 2006)',
    'danRer5':'Zebrafish Genome (July 2007)',
    'dasNov1':'Armadillo Genome (May 2005)',
    'dasNov2':'Armadillo Genome (July 2008)',
    'dipOrd1':'Kangaroo rat Genome (July 2008)',
    'dm2':'D. melanogaster Genome (April 2004)',
    'dm3':'D. melanogaster Genome (April 2006)',
    'dp3':'D. pseudoobscura Genome (November 2004)',
    'dp4':'D. pseudoobscura Genome (February 2006)',
    'droAna1':'D. ananassae Genome (July 2004)',
    'droAna2':'D. ananassae Genome (August 2005)',
    'droAna3':'D. ananassae Genome (February 2006)',
    'droEre1':'D. erecta Genome (August 2005)',
    'droEre2':'D. erecta Genome (February 2006)',
    'droGri1':'D. grimshawi Genome (August 2005)',
    'droGri2':'D. grimshawi Genome (February 2006)',
    'droMoj1':'D. mojavensis Genome (August 2004)',
    'droMoj2':'D. mojavensis Genome (August 2005)',
    'droMoj3':'D. mojavensis Genome (February 2006)',
    'droPer1':'D. persimilis Genome (October 2005)',
    'droSec1':'D. sechellia Genome (October 2005)',
    'droSim1':'D. simulans Genome (April 2005)',
    'droVir1':'D. virilis Genome (July 2004)',
    'droVir2':'D. virilis Genome (August 2005)',
    'droVir3':'D. virilis Genome (February 2006)',
    'droWil1':'D. willistoni Genome (February 2006)',
    'droYak1':'D. yakuba Genome (April 2004)',
    'droYak2':'D. yakuba Genome (November 2005)',
    'echTel1':'Tenrec Genome (July 2005)',
    'eriEur1':'European Hedgehog (Junuary 2006)',
    'equCab1':'Horse Genome (January 2007)',
    'equCab2':'Horse Genome (September 2007)',
    'felCat3':'Cat Genome (March 2006)',
    'fr1':'Fugu Genome (August 2002)',
    'fr2':'Fugu Genome (October 2004)',
    'galGal2':'Chicken Genome (February 2004)',
    'galGal3':'Chicken Genome (May 2006)',
    'gasAcu1':'Stickleback Genome (February 2006)',
    'gorGor1':'Gorilla Genome (October 2008)',
    'hg17':'Human Genome (May 2004)',
    'hg18':'Human Genome (May 2006)',
    'loxAfr1':'Elephant Genome (May 2005)',
    'loxAfr2':'Elephant Genome (July 2008)',
    'micMur1':'Mouse lemur Genome (June 2003)',
    'mm5':'Mouse Genome (May 2004)',
    'mm6':'Mouse Genome (March 2005)',
    'mm7':'Mouse Genome (August 2005)',
    'mm8':'Mouse Genome (March 2006)',
    'mm9':'Mouse Genome (July 2007)',
    'monDom1':'Opossum Genome (October 2004)',
    'monDom2':'Opossum Genome (June 2005)',
    'monDom4':'Opossum Genome (January 2006)',
    'myoLuc1':'Microbat Genome (March 2006)',
    'ochPri2':'Pika Genome (July 2008)',
    'ornAna1':'Platypus Genome (March 2007)',
    'oryCun1':'Rabbit Genome (May 2005)',
    'oryLat1':'Medaka Genome (April 2006)',
    'oryLat2':'Medaka Genome (October 2005)',
    'otoGar1':'Bushbaby Genome (December 2006)',
    'panTro1':'Chimpanzee Genome (November 2003)',
    'panTro2':'Chimpanzee Genome (March 2006)',
    'petMar1':'Lamprey Genome (March 2007)',
    'ponAbe2':'Orangutan Genome (July 2007)',
    'priPac1':'P. pacificus Genome (February 2007)',
    'proCap1':'Rock hyrax Genome (July 2008)',
    'pteVam1':'Megabat Genome (July 2008)',
    'rheMac1':'Rhesus Genome (January 2005)',
    'rheMac2':'Rhesus Genome (January 2006)',
    'rn3':'Rat Genome (June 2003)',
    'rn4':'Rat Genome (November 2004)',
    'sacCer1':'Yeast (S. cerevisiae) Genome (October 2003)',
    'sorAra1':'Shrew (Junuary 2006)',
    'speTri1':'Squirrel Genome (Febuary 2008)',
    'strPur1':'S. purpuratus Genome (April 2005)',
    'taeGut1':'Zebra finch Genome (July 2008)',
    'tarSyr1':'Tarsier Genome (August 2008)',
    'tetNig1':'Tetraodon Genome (February 2004)',
    'tupBel1':'Tree Shrew (December 2006)',
    'triCas2':'T. castaneum Genome (September 2005)',
    'turTru1':'Dolphin Genome (Febuary 2008)',
    'vicPac1':'Alpaca Genome (July 2008)',
    'xenTro1':'X. tropicalis Genome (October 2004)',
    'xenTro2':'X. tropicalis Genome (August 2005)'
    }

if __name__ == '__main__':

    pygrDir = '/data/server'
    dnDb = '/data/server/downloadable/.pygr_data'
    seqdir = '/data/GENOMES'
    msadir = '/data/NLMSA'

    mdb = metabase.MetabaseList(pygrDir)
    genomeList = [ix for ix in mdb.dir('Bio.Seq.Genome') if ix[-6:] != '.fasta' and ix[-4:] != '.txt']

    seqlist = glob.glob(os.path.join(seqdir, '*.seqlen'))
    for seqname in seqlist:
        genoname = os.path.basename(seqname)[:-7]
        # UPDATE EACH TIME IT EXITS
        if not docstringdict.has_key(genoname) or not sprotDict.has_key(genoname):
            sys.exit('NO DOCSTRING: %s' % genoname)
        resourceName = 'Bio.Seq.Genome.' + sprotDict[genoname] + '.' + genoname
        if resourceName in genomeList: continue
        print 'NOT REGISTERED: %s' % genoname
        genome = seqdb.SequenceFileDB(os.path.join(seqdir, genoname))
        genome.__doc__ = docstringdict[genoname]
        mdb.add_resource(resourceName, genome)
        print genoname + '\t' + 'Bio.Seq.Genome.' + sprotDict[genoname] + '.' + genoname + '\tREGISTERED'
        mdb.commit()

    print 'SEQDB REGISTRATION DONE'

