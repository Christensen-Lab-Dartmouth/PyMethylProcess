class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: meffil_encode
baseCommand: []
inputs:
  - id: input_sample_sheet
    type: File
  - id: idat_dir_input
    type: Directory
  - id: encode
    type: boolean?
outputs:
  - id: output_sample_sheet
    type: File
    outputBinding:
      glob: geo_idats/$(inputs.input_sample_sheet.basename)
  - id: idat_dir
    type: Directory
    outputBinding:
      glob: geo_idats
label: meffil_encode
arguments:
  - position: 0
    prefix: mkdir geo_idats && mv
    shellQuote: false
    valueFrom: $(inputs.idat_dir_input.path+'/*') geo_idats
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${
          if (inputs.encode) {

             return " && pymethyl-preprocess meffil_encode -os geo_idats/"+inputs.input_sample_sheet.basename + " -is " + "geo_idats/"+inputs.input_sample_sheet.basename
              
          }
          else {
              return ""
              
          }
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 2000
  - class: DockerRequirement
    dockerPull: 'pymethylprocess:latest'
  - class: InlineJavascriptRequirement
