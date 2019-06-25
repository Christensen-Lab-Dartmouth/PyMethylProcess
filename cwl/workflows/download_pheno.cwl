class: Workflow
cwlVersion: v1.0
id: download_pheno
label: download_pheno
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: query
    type: string
    'sbg:x': -276
    'sbg:y': -74
outputs:
  - id: output_pheno_sheet
    outputSource:
      - download_pheno/output_pheno_sheet
    type: File
    'sbg:x': 34
    'sbg:y': -71
steps:
  - id: download_pheno
    in:
      - id: query
        source: query
    out:
      - id: output_pheno_sheet
    run: ../tools/download_pheno.cwl
    label: download_pheno
    'sbg:x': -133
    'sbg:y': -73
requirements: []
