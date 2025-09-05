import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { bam2fastq_subworkflow } from '/home/nmahfel/nf-training/nf-core-nanoraredx/tests/../subworkflows/local/bam2fastq.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()

workflow {

    // run dependencies
    

    // workflow mapping
    def input = []
    
                input[0] = [
                    [ id:'test_methyl_single' ],
                    [
                        file("/home/nmahfel/nf-training/nf-core-nanoraredx/assets/test_data/bam_pass/test_pass_fc1e9677_854ad362_0.bam", checkIfExists: true)
                    ]
                ]
                // ch_fasta - empty channel
                input[1] = [[:], []]
                // ch_fai - empty channel  
                input[2] = [[:], []]
                
    //----

    //run workflow
    bam2fastq_subworkflow(*input)
    
    if (bam2fastq_subworkflow.output){

        // consumes all named output channels and stores items in a json file
        for (def name in bam2fastq_subworkflow.out.getNames()) {
            serializeChannel(name, bam2fastq_subworkflow.out.getProperty(name), jsonOutput)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = bam2fastq_subworkflow.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
}


def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
