# PanGenomeTools CWL Tools

This directory contains CWL (Common Workflow Language) tool definitions for PanGenomeTools. These tools can be used in CWL workflows to create reproducible, portable genomic analysis pipelines.

## Available Tools

### 1. `extract_genes_from_fasta.cwl`
Extracts sequences from FASTA files using GFF coordinates.

**Key Features:**
- Uses PanGenomeTools FastaHandler for efficient sequence extraction
- Supports upstream/downstream regions, inner regions, and padding
- Multiple merge strategies for handling multiple features
- Can include coordinates in output headers

### 2. `extract_genes_from_gff.cwl`
Extracts GFF features for target genes.

**Key Features:**
- Uses PanGenomeTools GFFHandler for efficient feature extraction
- Multiple search modes (strict, children, pattern)
- Multiple merge strategies for handling multiple features

### 3. `extract_genes_from_bigwig.cwl`
Extracts signal values from BigWig files using GFF coordinates.

**Key Features:**
- Uses PanGenomeTools BigWigHandler for efficient signal extraction
- Supports upstream/downstream regions and inner regions
- Can extract entire feature sequences
- Outputs signals in JSON format

### 4. `make_pseudo_annotation.cwl`
Creates pseudo GFF annotations from FASTA files.

**Key Features:**
- Useful when no annotations are available
- Creates minimal GFF entries for each sequence
- Configurable upstream/downstream regions

### 5. `extract_all_genes_from_fasta.cwl`
Extracts all genes from FASTA files using AGAT for multiple genomes.

**Key Features:**
- Runs AGAT transcriptome extraction for each genome
- Processes multiple genomes from a pangenome index
- Handles logging and error reporting

## Usage Examples

### Basic Usage

To run a tool directly:

```bash
cwl-runner extract_genes_from_fasta.cwl input.json
```

### Example Input JSON

Here's an example input JSON for `extract_genes_from_fasta.cwl`:

```json
{
  "pangenome_folder": {
    "class": "Directory",
    "path": "/path/to/pangenome"
  },
  "pangenome_index": {
    "class": "File",
    "path": "/path/to/pangenome_index.csv"
  },
  "target_genes": {
    "class": "File",
    "path": "/path/to/target_genes.csv"
  },
  "output": "extracted_sequences.fasta",
  "feature_type": "gene",
  "upstream": 1000,
  "downstream": 500,
  "merge_strategy": "merge"
}
```

### Using in Workflows

These tools can be combined in CWL workflows. Here's a simple example:

```cwl
cwlVersion: v1.0
class: Workflow

inputs:
  pangenome_folder: Directory
  pangenome_index: File
  target_genes: File

steps:
  extract_gff:
    run: extract_genes_from_gff.cwl
    in:
      pangenome_folder: pangenome_folder
      pangenome_index: pangenome_index
      target_genes: target_genes
      output: "extracted_features.gff"
    out: [extracted_features]

  extract_fasta:
    run: extract_genes_from_fasta.cwl
    in:
      pangenome_folder: pangenome_folder
      pangenome_index: pangenome_index
      target_genes: target_genes
      output: "extracted_sequences.fasta"
    out: [extracted_sequences]

outputs:
  gff_features:
    type: File
    outputSource: extract_gff/extracted_features
  fasta_sequences:
    type: File
    outputSource: extract_fasta/extracted_sequences
```

## Implementation Details

### Direct Function Calls

Most tools use **direct Python function calls** rather than command-line invocation. This approach:

- **Improves performance** by avoiding CLI overhead
- **Enhances reliability** by using the same code path as the Python API
- **Maintains compatibility** with existing CLI interfaces
- **Simplifies debugging** with direct Python error messages

### Tool Structure

Each CWL tool follows this structure:

1. **Metadata**: Clear labels and documentation
2. **Inputs**: Well-documented parameters with defaults
3. **Outputs**: Clearly defined output files
4. **Implementation**: Direct function calls where possible

### Error Handling

Tools inherit the error handling from the underlying Python functions:
- Missing files are reported clearly
- Invalid parameters are validated
- Progress can be suppressed with `--silent` flag

## Best Practices

### Input Validation

Always validate your inputs:
- Ensure paths are correct and accessible
- Check that required files exist
- Validate parameter ranges

### Output Management

- Use descriptive output filenames
- Consider using separate directories for different output types
- Be aware of file size limitations

### Performance Considerations

- For large datasets, consider breaking into smaller batches
- Use appropriate merge strategies for your use case
- Be mindful of memory usage with large upstream/downstream regions

## Support

For issues or questions about these CWL tools:

1. Check the main PanGenomeTools documentation
2. Review the underlying Python code for detailed parameter descriptions
3. Test with small datasets before running large analyses

## License

These CWL tools are licensed under the same license as PanGenomeTools.