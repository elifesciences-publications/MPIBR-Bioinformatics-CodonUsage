#!/usr/bin/python3

from Bio import SeqIO
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

obj = CodonUsage.CodonAdaptationIndex();
obj.generate_index("mouse.rna.cds_clean.fa");

obj.print_index();
