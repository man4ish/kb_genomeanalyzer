import os
import subprocess
import shutil

from pprint import pprint,pformat

class GATKUtils:

    def __init__(self):
       self.path = "/kb/module/data/"
       self.scratch = "/kb/module/work/tmp/"
       #self.path = "/home/manish/Desktop/apps/kb_genomeanalyzer/data/"
       #self.scratch = "/home/manish/Desktop/apps/kb_genomeanalyzer/test_local/workdir/tmp/"
       pass 

    def run_cmd(self, cmd):
        try:
           process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
           stdout, stderr = process.communicate()
           if stdout:
               print ("ret> ", process.returncode)
               print ("OK> output ", stdout)
           if stderr:
               print ("ret> ", process.returncode)
               print ("Error> error ", stderr.strip())

        except OSError as e:
           print ("OSError > ", e.errno)
           print ("OSError > ", e.strerror)
           print ("OSError > ", e.filename)
    
    def build_genome(self, assembly_file ):
        cmd = "bwa index -a bwtsw " + assembly_file
        self.run_cmd(cmd)
       
    def index_assembly(self, assembly_file):
        cmd = "samtools faidx " + assembly_file
        self.run_cmd(cmd)
      
    def generate_sequence_dictionary(self, assembly_file):
        cmd = "java -jar "+ self.path + "picard.jar CreateSequenceDictionary REFERENCE=" + assembly_file +" OUTPUT=" + assembly_file.replace("fa","dict")
        self.run_cmd(cmd)
               
    def mapping_genome(self, ref_genome, rev_fastq, fwd_fastq ):
        cmd = "bwa mem -t 32 -M -R " + "\"@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1\" " + ref_genome +" "+ rev_fastq + " " + fwd_fastq + " > " + self.scratch + "aligned_reads.sam"
        self.run_cmd(cmd)
    
    def duplicate_marking(self):
        cmd = "java -jar "+ self.path + "picard.jar SortSam  INPUT= " + self.scratch + "aligned_reads.sam   OUTPUT=" + self.scratch + "aligned_reads.bam  SORT_ORDER=coordinate"
        self.run_cmd(cmd)
       
    def sort_bam_index(self):
        cmd = "samtools index " + self.scratch + "aligned_reads.bam"
        self.run_cmd(cmd)
        
    def collect_alignment_and_insert_size_metrics(self, assembly_file):
        cmd1 = "java -jar "+ self.path + "picard.jar CollectAlignmentSummaryMetrics R="+ assembly_file  +" I=" + self.scratch + "aligned_reads.bam O=" + self.scratch + "alignment_metrics.txt"
        self.run_cmd(cmd1)
        cmd2 = "java -jar "+ self.path + "picard.jar CollectInsertSizeMetrics INPUT=" + self.scratch + "aligned_reads.bam OUTPUT=" + self.scratch + "insert_metrics.txt HISTOGRAM_FILE=" + self.scratch + "insert_size_histogram.pdf"
        self.run_cmd(cmd2)
        cmd3 = "samtools depth -a " + self.scratch + "aligned_reads.bam > " + self.scratch + "depth_out.txt"
        self.run_cmd(cmd3)

    def variant_calling(self, assembly_file):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar HaplotypeCaller -R "+ assembly_file  + " -I " + self.scratch + "aligned_reads.bam -O " + self.scratch + "raw_variants.vcf"
        self.run_cmd(cmd)

    def extract_variants(self, assembly_file ):
        cmd1 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file  + " -V " + self.scratch + "raw_variants.vcf --select-type SNP -O " + self.scratch + "raw_snps.vcf"
        self.run_cmd(cmd1)
        cmd2 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file  + " -V " + self.scratch + "raw_variants.vcf --select-type INDEL -O " + self.scratch + "raw_indels.vcf"
        self.run_cmd(cmd2)

    def filter_SNPs(self, assembly_file, output_file):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file  + " -V " + self.scratch + "raw_snps.vcf -O "  + self.scratch +  output_file +" -filter-name 'QD_filter' -filter 'QD < 2.0' -filter-name 'FS_filter' -filter 'FS > 60.0' -filter-name 'MQ_filter' -filter 'MQ < 40.0' -filter-name 'SOR_filter' -filter 'SOR > 4.0' -filter-name 'MQRankSum_filter' -filter 'MQRankSum < -12.5' -filter-name 'ReadPosRankSum_filter' -filter 'ReadPosRankSum < -8.0'"
        self.run_cmd(cmd)

    def filter_Indels(self, assembly_file, output_file):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file  + " -V " + self.scratch + "raw_indels.vcf -O " + self.scratch + output_file +" -filter-name 'QD_filter' -filter 'QD < 2.0' -filter-name 'FS_filter' -filter 'FS > 200.0' -filter-name 'SOR_filter' -filter 'SOR > 10.0'"
        self.run_cmd(cmd)

    def exclude_filtered_variants(self):
        cmd1 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V "  + self.scratch +"filtered_snps.vcf -O "  + self.scratch +"bqsr_snps.vcf"
        self.run_cmd(cmd1)

        cmd2 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V "  + self.scratch + "filtered_indels.vcf -O "  + self.scratch +"bqsr_indels.vcf"
        self.run_cmd(cmd2)

    def base_quality_score_recalibration(self, assembly_file, data_table):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar BaseRecalibrator -R " + assembly_file  + " -I "  + self.scratch + "aligned_reads.bam --known-sites "  + self.scratch + "bqsr_snps.vcf --known-sites "  + self.scratch + "bqsr_indels.vcf -O "  + self.scratch + data_table
        self.run_cmd(cmd)

    def apply_BQSR(self, assembly_file, data_table):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ApplyBQSR -R " + assembly_file  + "  -I "  + self.scratch +"aligned_reads.bam -bqsr "+ self.scratch + data_table + " -O "  + self.scratch +"recal_reads.bam"
        self.run_cmd(cmd)

    def analyze_covariates(self):
        #Error in library(gplots) : there is no package called ‘gplots’
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar AnalyzeCovariates -before "  + self.scratch +"recal_data.table -after "  + self.scratch +"post_recal_data.table -plots "  + self.scratch +"recalibration_plots.pdf"
        self.run_cmd(cmd)

    def annotate_SNPs_and_predict_effects(self):
        cmd = "java -jar "+ self.path + "snpEff/snpEff.jar -v Arabidopsis_thaliana "  + self.scratch +"filtered_snps_final.vcf > "  + self.scratch +"filtered_snps_final.ann.vcf"
        self.run_cmd(cmd)
    
    def compile_statistics(self):
        cmd = "parse_metrics.sh sample_id > sample_id_report.csv"
        self.run_cmd(cmd)
