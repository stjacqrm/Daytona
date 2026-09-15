#!/usr/bin/env nextflow

/*
Note:
Before running the script, please set the parameters in the config file params.yaml
*/

//Step1:input data files
nextflow.enable.dsl=2
def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

myDir.eachFileMatch ~/.*.fastq/, {L001R1Lst << it.name}
L001R1Lst.sort()
L001R1Lst.each{
	//println it
   def x =[]
   x = it.split('\\.')
   //println x
   sampleNames.add(x[0])
}
//println L001R1Lst
//println sampleNames

//create a new directory
myDir2 = file("${params.output}/assemblies_pass")
myDir2.mkdirs()
myDir3 = file("${params.output}/vadr_error_reports_clearlabs")
myDir3.mkdirs()


//Step2: process the inputed data
A = Channel.fromList(sampleNames)
//A.view()


include { pystats2 } from './modules/pystats2.nf'

workflow {
    pystats2(A) | view
}
