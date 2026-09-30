#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

inputs:
  csv:
    type: File
    inputBinding:
      position: 1
      prefix: --csv
    doc: "Path to CSV file with paired read files (columns: sample, read)"

baseCommand: [/usr/local/bin/_entrypoint.sh, python, -c, "from pangenometools.cli.csv_cli import main; main()"]

outputs:
  sample:
    type:
      type: array
      items: string
    outputBinding:
      glob:
      - sample.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents))
  read:
    type:
      type: array
      items: string
    outputBinding:
      glob:
      - read.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents))

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 1024
