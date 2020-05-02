import os
import subprocess

from pprint import pprint,pformat

class GATKUtils:

    def __init__(self):
       self.callbackURL = os.environ['SDK_CALLBACK_URL']
       pass 

    def deinterleave(self, fastq_file):
        path = "/kb/module/work/tmp/"
        fastq_1 = open(path + "r.fastq",'w')
        fastq_2 = open(path + "f.fastq",'w')
        [fastq_1.write(line) if (i % 8 < 4) else fastq_2.write(line) for i, line in enumerate(open(fastq_file))]
        fastq_1.close()
        fastq_2.close()
 
    def build_genome(assembly_file,  ):
        return 

    def mapping_genome():
        return 

    def duplicate_marking():
        return

    def indel_realignment():
        return

    def base_quality_score_recalibration():
        return

    def variant_calling():
        return

    def genotype_refinement():
        return

    def annotation():
        return

