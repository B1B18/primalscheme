global_args = {
    'PRIMER_OPT_SIZE': 22,
    'PRIMER_MIN_SIZE': 22,
    'PRIMER_MAX_SIZE': 30,
    'PRIMER_OPT_TM': 61.5,
    'PRIMER_MIN_TM': 60.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 30.0,
    'PRIMER_MAX_GC': 55.0,
    'PRIMER_MAX_POLY_X': 5,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY_TH': 47.0,
    'PRIMER_MAX_SELF_END_TH': 47.0,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 47.0,
    'PRIMER_PAIR_MAX_COMPL_END_TH': 47.0,
    'PRIMER_MAX_HAIRPIN_TH': 47.0,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
}

seq_args = {
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [-1, -1, -1, -1],
    }

MATCHES = [
    set(['A', 'T']),
    set(['C', 'G']),
    set(['G', 'T']),
    set(['C', 'T']),
    set(['T', 'T'])
],

MISMATCHES = [
    set(['A', 'A']),
    set(['A', 'C']),
    set(['C', 'C']),
    set(['G', 'A']),
    set(['G', 'G']),
]

NATIVE_DICT = {
    'NB01': 'AAGAAAGTTGTCGGTGTCTTTGTG',
    'NB02': 'TCGATTCCGTTTGTAGTCGTCTGT',
    'NB03': 'GAGTCTTGTGTCCCAGTTACCAGG',
    'NB04': 'TTCGGATTCTATCGTGTTTCCCTA',
    'NB05': 'CTTGTCCAGGGTTTGTGTAACCTT',
    'NB06': 'TTCTCGCAAAGGCAGAAAGTAGTC',
    'NB07': 'GTGTTACCGTGGGAATGAATCCTT',
    'NB08': 'TTCAGGGAACAAACCAAGTTACGT',
    'NB09': 'AACTAGGCACAGCGAGTCTTGGTT',
    'NB10': 'AAGCGTTGAAACCTTTGTCCTCTC',
    'NB11': 'GTTTCATCTATCGGAGGGAATGGA',
    'NB12': 'CAGGTAGAAAGAAGCAGAATCGGA',
}

SISPA_PRIMER = 'TTAACCGGCAGTGACACTGC'

RLBseq = 'TTTTTCGTGCGCCGCTTCAAC'

# the tail sequences get appended to the target primer sequences
TAIL_P5 = 'CCCTACACGACGCTCTTCCGATCT'  # T_m 67C at 1mM Mg2+, 0.2 mM dNTP
TAIL_N7 = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
ENRICHMENT_P5 = 'AATGATACGGCGACCACCGAGATCTACACAAAAAAAAACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ENRICHMENT_N7 = 'CAAGCAGAAGACGGCATACGAGATAAAAAAAAGTCTCGTGGGCTCGG'

ENRICHMENT_N7_STUB = 'GTCTCGTGGGCTCGGAGATGTG*T*A*T'
ENRICHMENT_F7_STUB = 'CAAGCAGAAGACGGCATACGAGAT'
ENRICHMENT_P5_STUB = 'ACACTCTTTCCCTACACGACGCTCTTCCGA*T*C*T'
ENRICHMENT_F5_STUB = 'AATGATACGGCGACCACCGAGATCTACAC'

INLINE_BARCODES = [
'TCGTTA','CATGAT','GTCATT','TGCGCA','GTCTCG','ACTTAG','GGTTAA','CTAAGC',
'GATATG','ACGAGT','CGATGA','AGAGGT','CTCCGT','GACGTC','GCATTG','GCTGCT',
'AAGGTT','TCAATC','TAACCA','TGCTAT','CCTCAA','AAGTAC','TCTCGG','ATGATG',
'CGGTTC','GGAACG','GCAGGC','GCCGAA','AGGTCG','AGTCTT','AGCTTA','CTCGTA',
'TAGCTC','GGACTC','AACCGA','TATAGA','TCCTCC','CCAGCG','CCGCTT','ATTCGC',
'ATTGCG','GCTAAC','CGTAGG','ATAGAA','AGCAAC','TTGCCG','TGCCTG','CTATAG',
'TACGAG','CTGGCT','GTTGGA','ATTAAT','AATACC','CAGTCA','GCGCCA','TGAGAC',
'TTCCAA','CAGAAG','GAACGG','CATCCG','GAATCC','GAGAGC','CGTTCT','CCGGAC',
'TTCGGC','AGGCAA','GCGTAT','TGGACT','TTGAAC','TTAGTT','ATACCT','GCCAGG',
'AACTCT','GACCAT','ACTATA','GTTCAG','GAGGCG','ATATTC','CCTTGC','CGAATT']


ENRICHMENT_BARCODES_N7 = {
'ATCAATCG':'CGATTGAT','TACCTCAA':'TTGAGGTA','CGCATAGA':'TCTATGCG',
'AAGCTATA':'TATAGCTT','ATAACGTA':'TACGTTAT','TTCCTTAC':'GTAAGGAA',
'TTACCTCT':'AGAGGTAA','CCTAGCTC':'GAGCTAGG','CGGACTTC':'GAAGTCCG',
'TTGATTAT':'ATAATCAA','TGCTCATA':'TATGAGCA','AACTATAA':'TTATAGTT',
'CTATCAGG':'CCTGATAG','GAAGGAGG':'CCTCCTTC','CCATTATA':'TATAATGG',
'ATATCTAA':'TTAGATAT','AGTTCTAT':'ATAGAACT','TATTATAC':'GTATAATA',
'TCAATGAG':'CTCATTGA','CAACGGAC':'GTCCGTTG','TCATTCGG':'CCGAATGA',
'ATAAGATT':'AATCTTAT','AGAGTCGT':'ACGACTCT','CTATGGTA':'TACCATAG',
'GGAATGAA':'TTCATTCC','CGATGATG':'CATCATCG','TAGATCTG':'CAGATCTA',
'AACTCAAC':'GTTGAGTT','CGACTCGA':'TCGAGTCG','TCAAGCTA':'TAGCTTGA',
'GTTAGGAT':'ATCCTAAC','CTCATACT':'AGTATGAG','GAATACTG':'CAGTATTC',
'CTCTCCAA':'TTGGAGAG','GAATATGA':'TCATATTC','CTACTATC':'GATAGTAG',
'TAGCTTGG':'CCAAGCTA','GCAATCCG':'CGGATTGC','ACGTACGA':'TCGTACGT',
'TCATATCA':'TGATATGA','ACGACTCG':'CGAGTCGT','CATCTACG':'CGTAGATG',
'AGCTAGTC':'GACTAGCT','AATACCGT':'ACGGTATT','GAGGTAAG':'CTTACCTC',
'CTATAACT':'AGTTATAG','GAGTATAG':'CTATACTC','CCTCATCA':'TGATGAGG',
'TCAGTAGT':'ACTACTGA','TAACTAAG':'CTTAGTTA','ACCTATGC':'GCATAGGT',
'TCATCAAG':'CTTGATGA','CTAGAAGA':'TCTTCTAG','GGATCAGA':'TCTGATCC',
'TGGAAGAC':'GTCTTCCA','CAATGAGT':'ACTCATTG','ATACTACT':'AGTAGTAT',
'CTTCCTTC':'GAAGGAAG','TCTATATC':'GATATAGA','GAGAGTAT':'ATACTCTC',
'CGTTGAAC':'GTTCAACG','GGAAGAGC':'GCTCTTCC','GAAGTCTT':'AAGACTTC',
'GGCTATTG':'CAATAGCC','ATGAACTA':'TAGTTCAT','TAAGCCTG':'CAGGCTTA',
'GGTCTAGC':'GCTAGACC','CATCTTAA':'TTAAGATG','AATCCTTA':'TAAGGATT',
'TATTCAAT':'ATTGAATA','TAGACCAT':'ATGGTCTA','AATGAGGT':'ACCTCATT',
'CAGGCTAA':'TTAGCCTG','ACTTGATT':'AATCAAGT','AATCTCCA':'TGGAGATT',
'GTCCGAAT':'ATTCGGAC','AATTCGTT':'AACGAATT','ATCTCTCT':'AGAGAGAT',
'CATCAGAG':'CTCTGATG','GTCAACGG':'CCGTTGAC','ATCCGATC':'GATCGGAT',
'CATTAATC':'GATTAATG','TAGTTATC':'GATAACTA','AAGAATCT':'AGATTCTT',
'CTCTCATC':'GATGAGAG','ATATAATC':'GATTATAT','TTATATAT':'ATATATAA',
'GATCGAAG':'CTTCGATC','CGTTAGGT':'ACCTAACG','CCATTCCT':'AGGAATGG',
'GCTTGGAG':'CTCCAAGC','TAATTACT':'AGTAATTA','CGGAATAG':'CTATTCCG',
'TATAGGTT':'AACCTATA','GTTCTCAA':'TTGAGAAC','GCTAATAT':'ATATTAGC'
}

ENRICHMENT_BARCODES_P5 = {
'AGTAGGTC':'AGTAGGTC','AGATGGAG':'AGATGGAG','CGATAACC':'CGATAACC',
'GAGTACGT':'GAGTACGT','GTATATGA':'GTATATGA','AAGTTATG':'AAGTTATG',
'GAGATCCT':'GAGATCCT','CAGATCAA':'CAGATCAA','AGGTCTTC':'AGGTCTTC',
'AGTATGAA':'AGTATGAA','CTGAGATT':'CTGAGATT','CTAGATCC':'CTAGATCC',
'CTCTTATA':'CTCTTATA','CGATTCAT':'CGATTCAT','TTCGGAGC':'TTCGGAGC',
'GGCCTAGG':'GGCCTAGG','ATACTGAA':'ATACTGAA','CTTGAAGC':'CTTGAAGC',
'GCGGCCTA':'GCGGCCTA','CTAGTATT':'CTAGTATT','GATCTTGC':'GATCTTGC',
'TACTGCTT':'TACTGCTT','TGCAGTCA':'TGCAGTCA','AGCCTCCG':'AGCCTCCG',
'CAATTGCT':'CAATTGCT','AGAAGAAT':'AGAAGAAT','GCATATTG':'GCATATTG',
'GCTGAGAA':'GCTGAGAA','TCCTCGCG':'TCCTCGCG','TTGCCGTT':'TTGCCGTT',
'CCAACCTT':'CCAACCTT','CATCTGCC':'CATCTGCC','GTAGGAGG':'GTAGGAGG',
'ACTAATGA':'ACTAATGA','GTCTTCGG':'GTCTTCGG','ATAATTCT':'ATAATTCT',
'ATCTGGCC':'ATCTGGCC','CCGAATCA':'CCGAATCA','CCGAGGAA':'CCGAGGAA',
'GTATGAAT':'GTATGAAT','CAGCATAA':'CAGCATAA','ATCAAGCG':'ATCAAGCG',
'AGGTTAGC':'AGGTTAGC','CGTAGGCA':'CGTAGGCA','AGATCATG':'AGATCATG',
'GCATTCAA':'GCATTCAA','TACCTGAC':'TACCTGAC','GTAGCTCA':'GTAGCTCA',
'TACTTAGA':'TACTTAGA','CCGATATA':'CCGATATA','CAGATTCC':'CAGATTCC',
'TTGATATC':'TTGATATC','GCGAATAC':'GCGAATAC','GCATCATA':'GCATCATA',
'AGCTTCGA':'AGCTTCGA','TAATTATT':'TAATTATT','AAGAATCG':'AAGAATCG',
'AAGCAATT':'AAGCAATT','CGACCAGG':'CGACCAGG','ACTTCGCC':'ACTTCGCC',
'TGCCTTAA':'TGCCTTAA','TTATATTC':'TTATATTC','CCATCTCA':'CCATCTCA',
'GCCTGGTT':'GCCTGGTT','AGCTAGTA':'AGCTAGTA','TAAGATTG':'TAAGATTG',
'GCCGTACC':'GCCGTACC','CAGTCAGG':'CAGTCAGG','AGAGGATC':'AGAGGATC',
'GTACTTAT':'GTACTTAT','ACTTCTAA':'ACTTCTAA','CTCCTCAA':'CTCCTCAA',
'GTAGTAAC':'GTAGTAAC','GGTTGATA':'GGTTGATA','TTATCGTA':'TTATCGTA',
'TAGAACTA':'TAGAACTA','CTAGCGAC':'CTAGCGAC','CGGAACTC':'CGGAACTC',
'GTAGTCTA':'GTAGTCTA','GAGTAATC':'GAGTAATC','TGAATCGG':'TGAATCGG',
'TCTATGAT':'TCTATGAT','GCAATAAT':'GCAATAAT','AGAACTAC':'AGAACTAC',
'GAAGACGA':'GAAGACGA','CTGCCTAC':'CTGCCTAC','ACTTGAGG':'ACTTGAGG',
'AGGAGTAA':'AGGAGTAA','TATCCAAG':'TATCCAAG','TTCCTAAG':'TTCCTAAG',
'GTTAGGTT':'GTTAGGTT','TTATCCAG':'TTATCCAG','GACCAGGC':'GACCAGGC',
'TCATTGAC':'TCATTGAC','CGATAGGT':'CGATAGGT','CTCCAACC':'CTCCAACC'
}


MAX_HAIRPIN_TM = 50
MAX_HETERODIMER_TM = 45
MAX_HOMODIMER_TM = 45

# if more than this fraction of reference sequences differ from the oligo in
# a position then attempt to introduce degeneracy in the oligo to accommodate
ACCEPTABLE_DIVERGENCE = 0.02

# Length 22, 30, 24.43
# Tm 59.96, 62.72, 61.33
# GC 33.33, 54.55, 45.71
# Longest homopolymer 5
# Ambiguous 0
