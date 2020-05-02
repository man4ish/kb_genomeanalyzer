import os
import subprocess

from pprint import pprint,pformat

from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil


class DownloadFastqUtils:

    def __init__(self):
       self.callbackURL = os.environ['SDK_CALLBACK_URL']
       self.au = AssemblyUtil(self.callbackURL)
       pass 


    def _stage_input_file(self, ref, reads_type):

        ru = ReadsUtils(self.callbackURL)
        if reads_type == 'KBaseFile.PairedEndLibrary' or 'KBaseAssembly.PairedEndLibrary':
            input_file_info = ru.download_reads({
                    'read_libraries': [ref],
                    'interleaved': 'true'
                    })['files'][ref]
        elif reads_type == 'KBaseFile.SingleEndLibrary' or 'KBaseAssembly.SingleEndLibrary':
            input_file_info = ru.download_reads({
                    'read_libraries': [ref]
                    })['files'][ref]
        else:
            raise ValueError ("Can't download_reads() for object type: '"+str(reads_type)+"'")
        input_file_info['input_ref'] = ref
        file_location = input_file_info['files']['fwd']


        interleaved = False
        if input_file_info['files']['type'] == 'interleaved':
            interleaved = True
       
        return input_file_info

    def download_genome(self, genomeref):
        file = self.au.get_assembly_as_fasta({
          'ref': genomeref
        })
        return file

