class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: preprocess_pipeline
baseCommand: []
inputs:
  - id: idat_csv
    type: File
    inputBinding:
      position: 1
      prefix: '-i'
      shellQuote: false
      valueFrom: samplesheet_dir/
  - id: idat_dir
    type: Directory
  - id: n_cores
    type: int?
    inputBinding:
      position: 1
      prefix: '-n'
      shellQuote: false
  - id: meffil
    type: boolean?
    inputBinding:
      position: 1
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.meffil){
                return '-m'
            }
            else {
                return ''
            }
        }
outputs:
  - id: output_pkl
    type: File
    outputBinding:
      glob: '*/*array.pkl'
label: preprocess_pipeline
arguments:
  - position: 0
    prefix: mkdir geo_idats samplesheet_dir && cp
    shellQuote: false
    valueFrom: '$(inputs.idat_dir.path)/*.idat geo_idats/ && '
  - position: 1
    prefix: '&& pymethyl-preprocess preprocess_pipeline -o'
    shellQuote: false
    valueFrom: ./preprocess_outputs/methyl_array.pkl
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: cp $(inputs.idat_csv.path) samplesheet_dir/
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: |-
      ${
          if (inputs.n_cores){
          return inputs.n_cores+1
          }
          else {
              return 1
          }
      }
  - class: DockerRequirement
    dockerPull: 'pymethylprocess:latest'
  - class: InlineJavascriptRequirement
