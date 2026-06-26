#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

inputs:
  pangenome_index:
    type: File
    inputBinding:
      position: 1
      prefix: --pangenome-index
    doc: Path to pangenome index file

  genotypes:
    type: string[]?
    inputBinding:
      prefix: --genotypes
    doc: Genotype(s) to run the analysis for (all by default)

baseCommand: [/usr/local/bin/_entrypoint.sh, python, -c, "from pangenometools.cli.genotypes_cli import genotypes_cli; genotypes_cli()"]

outputs:
  genotypes:
    type:
      type: array
      items: string
    outputBinding:
      glob:
      - genotypes.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents))

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl:latest"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 1024