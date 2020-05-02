/*
A KBase module: kb_genomeanalyzer
*/

module kb_genomeanalyzer {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_kb_genomeanalyzer(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
