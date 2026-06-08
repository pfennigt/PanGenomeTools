# PanGenomeTools CWL Implementation Summary

## Overview

This document summarizes the CWL (Common Workflow Language) implementation for PanGenomeTools. The implementation provides CWL tool definitions that can be used in reproducible genomic workflows.

## Implementation Approach

### Design Principles

1. **Direct Function Calls**: Most tools use direct Python function calls rather than CLI invocation for better performance and reliability.

2. **Consistent Naming**: CWL tool names match the original Python scripts for easy identification.

3. **Comprehensive Documentation**: Each tool includes detailed metadata, parameter descriptions, and usage examples.

4. **Backward Compatibility**: All tools maintain compatibility with existing CLI interfaces.

### Tool Implementation Strategy

For tools that use the refactored CLI (fasta, gff, bigwig):
- Use `baseCommand: python`
- Call the underlying function directly using `stdin` with inline Python
- Create proper `Namespace` objects to pass parameters
- Maintain all original functionality and error handling

For tools that call external programs (extract_all_genes_from_fasta):
- Use traditional command-line invocation
- Properly handle file paths and arguments
- Include appropriate requirements (ShellCommandRequirement)

## File Structure

```
workflows/PanGenomeTools/cwl/
├── extract_genes_from_fasta.cwl          # FASTA sequence extraction
├── extract_genes_from_gff.cwl           # GFF feature extraction
├── extract_genes_from_bigwig.cwl        # BigWig signal extraction
├── make_pseudo_annotation.cwl           # Pseudo annotation creation
├── extract_all_genes_from_fasta.cwl     # Batch AGAT extraction
├── example_workflow.cwl                 # Example workflow
├── README.md                            # Usage documentation
├── CWL_IMPLEMENTATION_SUMMARY.md       # This file
└── examples/
    ├── extract_genes_from_fasta_input.json
    ├── extract_genes_from_gff_input.json
    ├── extract_genes_from_bigwig_input.json
    ├── make_pseudo_annotation_input.json
    ├── extract_all_genes_from_fasta_input.json
    └── workflow_input.json
```

## Tool Details

### 1. extract_genes_from_fasta.cwl

**Purpose**: Extract sequences from FASTA files using GFF coordinates

**Key Features**:
- Direct call to `extract_fasta_sequences()` function
- Supports all original parameters (upstream, downstream, merge strategies, etc.)
- Outputs FASTA file with extracted sequences

**Implementation**:
```cwl
baseCommand: python
stdin: |
  from pangenometools.cli.fasta_cli import extract_fasta_sequences
  from argparse import Namespace
  # Create Namespace and call function directly
```

### 2. extract_genes_from_gff.cwl

**Purpose**: Extract GFF features for target genes

**Key Features**:
- Direct call to `extract_gff_features()` function
- Supports search modes and merge strategies
- Outputs GFF file with extracted features

**Implementation**:
```cwl
baseCommand: python
stdin: |
  from pangenometools.cli.gff_cli import extract_gff_features
  from argparse import Namespace
  # Create Namespace and call function directly
```

### 3. extract_genes_from_bigwig.cwl

**Purpose**: Extract signal values from BigWig files

**Key Features**:
- Direct call to `extract_bigwig_signals()` function
- Supports coordinate adjustments and padding
- Outputs JSON file with signal data

**Implementation**:
```cwl
baseCommand: python
stdin: |
  from pangenometools.cli.bigwig_cli import extract_bigwig_signals
  from argparse import Namespace
  # Create Namespace and call function directly
```

### 4. make_pseudo_annotation.cwl

**Purpose**: Create pseudo GFF annotations from FASTA files

**Key Features**:
- Calls the standalone script directly
- Useful when no annotations are available
- Configurable feature regions

**Implementation**:
```cwl
baseCommand: python
arguments:
  - valueFrom: workflows/PanGenomeTools/scripts/make_pseudo_annotation.py
  # Pass arguments via command line
```

### 5. extract_all_genes_from_fasta.cwl

**Purpose**: Batch extraction using AGAT for multiple genomes

**Key Features**:
- Calls external AGAT tool via subprocess
- Processes multiple genomes from index
- Handles logging and error reporting

**Implementation**:
```cwl
baseCommand: python
arguments:
  - valueFrom: workflows/PanGenomeTools/scripts/extract_all_genes_from_fasta.py
  # Pass arguments via command line
```

## Usage Patterns

### Single Tool Execution

```bash
# Run a single tool
cwl-runner extract_genes_from_fasta.cwl examples/extract_genes_from_fasta_input.json
```

### Workflow Execution

```bash
# Run the example workflow
cwl-runner example_workflow.cwl examples/workflow_input.json
```

### Integration with Other Workflows

```cwl
# Import and use in other CWL workflows
steps:
  my_step:
    run: workflows/PanGenomeTools/cwl/extract_genes_from_fasta.cwl
    in:
      # ... inputs
    out: [output]
```

## Performance Considerations

### Advantages of Direct Function Calls

1. **Reduced Overhead**: No CLI parsing or subprocess creation
2. **Better Error Messages**: Direct Python exceptions with full stack traces
3. **Type Safety**: Parameters are passed as native Python types
4. **Consistent Behavior**: Same code path as Python API usage

### Memory Management

- Large datasets may require significant memory
- Consider batching for very large analyses
- Monitor memory usage with upstream/downstream parameters

## Error Handling

Tools inherit comprehensive error handling from the underlying Python code:

- **File Validation**: Checks for missing or invalid files
- **Parameter Validation**: Validates ranges and types
- **Progress Reporting**: Can be suppressed with `--silent` flag
- **Detailed Error Messages**: Clear reporting of issues

## Testing Recommendations

1. **Start Small**: Test with small datasets first
2. **Validate Inputs**: Ensure all paths are correct
3. **Check Outputs**: Verify output files are generated correctly
4. **Monitor Performance**: Watch memory and CPU usage
5. **Test Edge Cases**: Try with missing data, empty files, etc.

## Integration with Other Tools

These CWL tools are designed to be used with other genomic tools:

- **Input**: Can accept outputs from alignment, assembly, or annotation tools
- **Output**: Produces standard file formats (FASTA, GFF, JSON) for downstream analysis
- **Compatibility**: Works with any CWL-compliant workflow engine

## Future Enhancements

Potential improvements for future versions:

1. **Docker Integration**: Add Docker requirements for containerized execution
2. **Resource Requirements**: Add CPU/memory specifications
3. **Scatter/Gather**: Support for parallel processing of multiple genes
4. **Additional Tools**: CWL definitions for any new PanGenomeTools features
5. **Workflow Templates**: More complex example workflows

## Best Practices

### File Management

- Use absolute paths for reliability
- Create separate output directories for different analyses
- Clean up temporary files when done

### Parameter Tuning

- Adjust upstream/downstream based on your analysis needs
- Choose appropriate merge strategies for your data
- Consider padding for sequence analysis

### Reproducibility

- Version control your CWL files and input JSONs
- Document parameter choices and rationale
- Use consistent naming conventions

## Support and Troubleshooting

### Common Issues

1. **Missing Files**: Ensure all input files exist and paths are correct
2. **Permission Issues**: Check file permissions and directory access
3. **Memory Errors**: Reduce dataset size or increase available memory
4. **Parameter Errors**: Validate all parameters against documentation

### Debugging Tips

- Use `--debug` flag with CWL runner for detailed output
- Check intermediate files for partial results
- Test individual tools before combining in workflows
- Review error messages carefully for clues

## Conclusion

This CWL implementation provides a robust, reproducible way to use PanGenomeTools in workflow environments. The direct function call approach ensures optimal performance while maintaining full compatibility with existing functionality.

For questions or issues, refer to the main PanGenomeTools documentation or examine the underlying Python code for detailed implementation specifics.