class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: mytool
baseCommand:
  - python
  - /skeleton-cwl/mytool/mytool.py
inputs:
  - id: input_int
    type: int?
    inputBinding:
      position: 0
      prefix: '--input-int'
  - id: input_bool
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--input-bool'
  - id: input_string
    type: string?
    inputBinding:
      position: 0
      prefix: '--input-string'
  - id: input_file_a
    type: File
    inputBinding:
      position: 99
outputs:
  - id: output_a
    type: File?
    outputBinding:
      glob: *.output.ext
  - id: output_b
    type: File?
    outputBinding:
      glob: output_b.ext
label: mytool
arguments:
  - position: 0
    prefix: '--output-fp'
    valueFrom: output.ext
requirements:
  - class: DockerRequirement
    dockerPull: 'repo/mytool:0.0.1'
