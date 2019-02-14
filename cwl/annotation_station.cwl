class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: annotation-station
baseCommand:
  - python
  - /annotation-station/annotation-station/annotation_station.py
inputs:
  - id: input_type
    type: string?
    inputBinding:
      position: 0
      prefix: '--input-type'
  - id: input_header
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--input-header'
  - id: annotate_repeats
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--annotate-repeats'
  - id: repeats_table
    type: File?
    inputBinding:
      position: 0
      prefix: '--repeats-table'
  - id: annotate_transvar
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--annotate-transvar'
  - id: primary_transcripts
    type: File?
    inputBinding:
      position: 0
      prefix: '--primary-transcripts'
  - id: input_file
    type: File
    inputBinding:
      position: 99
outputs:
  - id: output_file
    type: File?
    outputBinding:
      glob: output.tsv
label: annotation-station
arguments:
  - position: 0
    prefix: '--output'
    valueFrom: output.tsv
requirements:
  - class: DockerRequirement
    dockerPull: 'estorrs/annotation-station:0.0.2'
