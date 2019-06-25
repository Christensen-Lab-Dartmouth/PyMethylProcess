class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: train_test_split
baseCommand:
  - pymethyl-utils train_test_val_split
inputs:
  - id: methyl_array
    type: File
    inputBinding:
      position: 0
      prefix: '-i'
      shellQuote: false
  - id: train_percent
    type: float
    inputBinding:
      position: 0
      prefix: '-tp'
      shellQuote: false
  - id: val_percent
    type: float
    inputBinding:
      position: 0
      prefix: '-vp'
      shellQuote: false
  - id: categorical
    type: boolean?
  - id: column_split_on
    type: string?
    inputBinding:
      position: 0
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.categorical) {
                return "-cat -k "+inputs.column_split_on
            }
            else {
               return "" 
            }
            return ""
            
        }
outputs:
  - id: train_methyl_arr
    type: File
    outputBinding:
      glob: '*/train*.pkl'
  - id: val_methyl_arr
    type: File
    outputBinding:
      glob: '*/val*.pkl'
  - id: test_methyl_arr
    type: File
    outputBinding:
      glob: '*/test*.pkl'
label: train_test_split
arguments:
  - position: 0
    prefix: '-o'
    shellQuote: false
    valueFrom: train_val_test_sets/
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pymethylprocess:latest'
  - class: InlineJavascriptRequirement
