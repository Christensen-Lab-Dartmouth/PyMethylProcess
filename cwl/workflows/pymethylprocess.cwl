class: Workflow
cwlVersion: v1.0
id: pymethylprocess
label: pymethylprocess
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: disease_class_column
    type: string?
    'sbg:x': -116.5560073852539
    'sbg:y': 339.029541015625
  - id: base_name_col
    type: string?
    'sbg:x': -133.2068634033203
    'sbg:y': 464.0591125488281
  - id: encode
    type: boolean?
    'sbg:x': 232.37869262695312
    'sbg:y': 419.82196044921875
  - id: n_neighbor
    type: int?
    'sbg:x': 2017.2542724609375
    'sbg:y': -40.810298919677734
  - id: include_columns_file
    type: File?
    'sbg:x': -92.7690658569336
    'sbg:y': 77.45567321777344
  - id: solver
    type: string?
    'sbg:x': 1232.29150390625
    'sbg:y': -228.28501892089844
  - id: imputation_method
    type: string?
    'sbg:x': 1244.4619140625
    'sbg:y': 287.550537109375
  - id: n_neighbors
    type: int?
    'sbg:x': 1242.911376953125
    'sbg:y': 27.75119400024414
  - id: query
    type: string
    'sbg:x': -246.3263397216797
    'sbg:y': 163.10784912109375
  - id: cpg_threshold
    type: float
    'sbg:x': 1238.1484375
    'sbg:y': 417.4735412597656
  - id: sample_threshold
    type: float
    'sbg:x': 1239.0537109375
    'sbg:y': -113.20138549804688
  - id: n_top_cpgs
    type: int?
    'sbg:x': 1672.260498046875
    'sbg:y': 38.22473907470703
  - id: val_percent
    type: float
    'sbg:x': 1929.2720947265625
    'sbg:y': 386.50311279296875
  - id: train_percent
    type: float
    'sbg:x': 1929.562255859375
    'sbg:y': 501.3965759277344
  - id: categorical
    type: boolean?
    'sbg:x': 1929.2900390625
    'sbg:y': 780.8756103515625
  - id: column_split_on
    type: string?
    'sbg:x': 1931.828125
    'sbg:y': 629.609619140625
  - id: n_cores
    type: int?
    'sbg:x': 716.2692260742188
    'sbg:y': -77.57911682128906
outputs:
  - id: output_visual
    outputSource:
      - transform_plot_1/output_visual
    type: File
    'sbg:x': 2355.1826171875
    'sbg:y': -32.67634963989258
  - id: test_methyl_arr
    outputSource:
      - train_test_split/test_methyl_arr
    type: File
    'sbg:x': 2435.85791015625
    'sbg:y': 383.20208740234375
  - id: train_methyl_arr
    outputSource:
      - train_test_split/train_methyl_arr
    type: File
    'sbg:x': 2438.12451171875
    'sbg:y': 255.71066284179688
  - id: val_methyl_arr
    outputSource:
      - train_test_split/val_methyl_arr
    type: File
    'sbg:x': 2423.485107421875
    'sbg:y': 130.1484832763672
steps:
  - id: create_sample_sheet
    in:
      - id: input_sample_sheet
        source: download_geo/initial_sample_sheet
      - id: source_type
        default: geo
      - id: idat_dir
        source: download_geo/idat_dir
      - id: header_line
        default: 0
      - id: disease_class_column
        source: disease_class_column
      - id: base_name_col
        source: base_name_col
      - id: include_columns_file
        source: include_columns_file
    out:
      - id: idat_dir_out
      - id: final_csv
    run: ../tools/create_sample_sheet.cwl
    label: Generate Sample Sheet
    'sbg:x': 230.09375
    'sbg:y': 266.6796875
  - id: meffil_encode
    in:
      - id: input_sample_sheet
        source: create_sample_sheet/final_csv
      - id: idat_dir_input
        source: create_sample_sheet/idat_dir_out
      - id: encode
        source: encode
    out:
      - id: output_sample_sheet
      - id: idat_dir
    run: ../tools/meffil_encode.cwl
    label: Format Sample Sheet
    'sbg:x': 462.36700439453125
    'sbg:y': 229.9230194091797
  - id: split_by_subtype
    in:
      - id: idat_csv
        source: meffil_encode/output_sample_sheet
    out:
      - id: output_csv
    run: ../tools/split_by_subtype.cwl
    label: 'Split Up Sample Sheet by Subtypes, Scatter Sample Sheets'
    'sbg:x': 718.5984497070312
    'sbg:y': 53.45710372924805
  - id: preprocess_pipeline_1
    in:
      - id: idat_csv
        source: split_by_subtype/output_csv
      - id: idat_dir
        source: meffil_encode/idat_dir
      - id: n_cores
        default: 8
        source: n_cores
      - id: meffil
        default: true
    out:
      - id: output_pkl
    run: ../tools/preprocess_pipeline.cwl
    label: Preprocessing Pipeline
    scatter:
      - idat_csv
    'sbg:x': 945.2371215820312
    'sbg:y': 189.73951721191406
  - id: combine_methyl_arrays
    in:
      - id: input_pkls
        source:
          - preprocess_pipeline_1/output_pkl
    out:
      - id: output_methyl_array
    run: ../tools/combine_methyl_arrays.cwl
    label: Combine Scattered MethylationArrays
    'sbg:x': 1235.7044677734375
    'sbg:y': 162.34913635253906
  - id: transform_plot_1
    in:
      - id: input_pkl
        source: mad_filter/output_pkl
      - id: column_of_interest
        default: disease
      - id: n_neighbor
        default: 8
        source: n_neighbor
      - id: axes_off
        default: false
    out:
      - id: output_visual
    run: ../tools/transform_plot.cwl
    label: Transform via UMAP and Plot via Plotly
    'sbg:x': 2150.485107421875
    'sbg:y': 13.538844108581543
  - id: imputation
    in:
      - id: input_pkl
        source: combine_methyl_arrays/output_methyl_array
      - id: imputation_method
        default: Mean
        source: imputation_method
      - id: solver
        default: sklearn
        source: solver
      - id: n_neighbors
        default: 5
        source: n_neighbors
      - id: sample_threshold
        default: 1
        source: sample_threshold
      - id: cpg_threshold
        default: 1
        source: cpg_threshold
    out:
      - id: imputed_methylarray
    run: ../tools/imputation.cwl
    label: Impute MethylationArray
    'sbg:x': 1684.01611328125
    'sbg:y': 213.3125
  - id: mad_filter
    in:
      - id: input_pkl
        source: imputation/imputed_methylarray
      - id: n_top_cpgs
        default: 300000
        source: n_top_cpgs
    out:
      - id: output_pkl
    run: ../tools/mad_filter.cwl
    label: Feature Selection
    'sbg:x': 1934
    'sbg:y': 208.91134643554688
  - id: download_geo
    in:
      - id: query
        source: query
    out:
      - id: idat_dir
      - id: initial_sample_sheet
    run: ../tools/download_geo.cwl
    label: Download GEO Methylation IDATs
    'sbg:x': -101.9169692993164
    'sbg:y': 213.77859497070312
  - id: train_test_split
    in:
      - id: methyl_array
        source: mad_filter/output_pkl
      - id: train_percent
        source: train_percent
      - id: val_percent
        source: val_percent
      - id: categorical
        source: categorical
      - id: column_split_on
        source: column_split_on
    out:
      - id: train_methyl_arr
      - id: val_methyl_arr
      - id: test_methyl_arr
    run: ../tools/train_test_split.cwl
    label: 'Split MethylationArray into Train, Val, Test'
    'sbg:x': 2230
    'sbg:y': 399.59796142578125
requirements:
  - class: ScatterFeatureRequirement
