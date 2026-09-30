cwlVersion: v1.2
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement

# Metadata
label: "Collect an array of Files as a Directory object"
doc: |
  Turns an array of File objects into a Directory object.
  When used as the output from a Workflow, this results in the files being returned in a folder.

inputs:
  files: File[]
  dir_name:
    type: string
    doc: "Name of the returned folder"
outputs:
  out: Directory
expression: |
  ${
    return {"out": {
      "class": "Directory", 
      "basename": inputs.dir_name,
      "listing": inputs.files
    } };
  }