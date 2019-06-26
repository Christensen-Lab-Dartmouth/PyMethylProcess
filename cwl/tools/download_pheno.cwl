class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: download_pheno
baseCommand: []
inputs:
  - id: query
    type: string
    inputBinding:
      position: 0
      prefix: '-g'
      shellQuote: false
      valueFrom: $(inputs.query) -of $(inputs.query)_pheno.csv
outputs:
  - id: output_pheno_sheet
    type: File
    outputBinding:
      glob: '*pheno.csv'
label: download_pheno
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: 'pymethyl-preprocess download_geo_clinical_info '
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pymethylprocess:latest'
  - class: InlineJavascriptRequirement
