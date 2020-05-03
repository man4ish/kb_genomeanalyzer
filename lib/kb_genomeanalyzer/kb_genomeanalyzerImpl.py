# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import shutil
from kb_genomeanalyzer.Utils.GATKUtils import GATKUtils
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class kb_genomeanalyzer:
    '''
    Module Name:
    kb_genomeanalyzer

    Module Description:
    A KBase module: kb_genomeanalyzer
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/man4ish/kb_genomeanalyzer.git"
    GIT_COMMIT_HASH = "44ecb6fc4ade9f0269bd9d162aebc29e669cfd8d"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.gu = GATKUtils()
        #END_CONSTRUCTOR
        pass


    def run_kb_genomeanalyzer(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_genomeanalyzer

        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})

        self.scratch = "/kb/module/work/tmp"
        src_assembly_file = "/kb/module/data/Athaliana_TAIR10.assembly.fa"
        shutil.copy(src_assembly_file, self.scratch)
        assembly_file = os.path.join(self.scratch,"Athaliana_TAIR10.assembly.fa")

        self.gu.build_genome(assembly_file)

        self.gu.index_assembly(assembly_file)

        self.gu.generate_sequence_dictionary(assembly_file)

        fwd_fastq = "/kb/module/data/arabidposis.1.fastq"
        rev_fastq = "/kb/module/data/arabidposis.2.fastq"

        self.gu.mapping_genome(assembly_file, fwd_fastq, rev_fastq )
  
        self.gu.duplicate_marking()

        self.gu.sort_bam_index()

        self.gu.collect_alignment_and_insert_size_metrics(assembly_file)
   
        self.gu.analyze_covariates()
   
        self.gu.variant_calling(assembly_file)
   
        self.gu.extract_variants(assembly_file)
   
        self.gu.filter_SNPs(assembly_file, "filtered_snps.vcf")
   
        self.gu.filter_Indels(assembly_file, "filtered_indels.vcf")
   
        self.gu.exclude_filtered_variants()
   
        self.gu.base_quality_score_recalibration(assembly_file, "recal_data.table")
   
        self.gu.apply_BQSR(assembly_file, "recal_data.table")
   
        self.gu.base_quality_score_recalibration(assembly_file, "post_recal_data.table")
   
        self.gu.apply_BQSR(assembly_file,  "post_recal_data.table")
   
        self.gu.filter_SNPs(assembly_file, "filtered_snps_final.vcf")
   
        self.gu.filter_Indels(assembly_file, "filtered_indels_final.vcf")

        self.gu.annotate_SNPs_and_predict_effects()

        self.gu.compile_statistics()

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kb_genomeanalyzer

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_genomeanalyzer return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
