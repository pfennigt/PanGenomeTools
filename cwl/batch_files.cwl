cwlVersion: v1.2
class: ExpressionTool

label: Batch files

inputs:
  files:
    type: File[]

  files2:
    type: File[]?

  batch_size:
    type: int

outputs:
  batches:
    type:
      type: array
      items:
        type: array
        items: File

expression: |
  ${
    var files = inputs.files.concat(inputs.files2 || []);
    var batches = [];

    for (var i = 0; i < files.length; i += inputs.batch_size) {
      batches.push(files.slice(i, i + inputs.batch_size));
    }

    return {"batches": batches};
  }

requirements:
  - class: InlineJavascriptRequirement