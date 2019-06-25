class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: mad_filter
baseCommand: []
inputs:
  - id: input_pkl
    type: File
    inputBinding:
      position: 0
      prefix: '-i'
      shellQuote: false
  - id: n_top_cpgs
    type: int?
    inputBinding:
      position: 0
      prefix: '-n'
      shellQuote: false
outputs:
  - id: output_pkl
    type: File
    outputBinding:
      glob: '*/*array.pkl'
label: mad_filter
arguments:
  - position: 0
    prefix: pymethyl-preprocess feature_select -o
    shellQuote: false
    valueFrom: '*/*array.pkl'
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 2000
  - class: DockerRequirement
    dockerPull: 'pymethylprocess:latest'
