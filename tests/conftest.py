"""
Pytest fixtures for PanGenomeTools tests.
"""

import pytest
from pathlib import Path
import tempfile
import shutil
import pyBigWig
import numpy as np

@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for testing."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)

@pytest.fixture
def setup_test_files(temp_test_dir):
    """Set up test files in a temporary directory."""
    # Create GFF file
    gff_content = """##gff-version 3
chr1	.	gene	100	500	.	+	.	ID=GENE001
chr1	.	mRNA	100	500	.	+	.	ID=GENE001.1;Parent=GENE001
chr1	.	CDS	150	200	.	+	.	ID=CDS001;Parent=GENE001.1
chr1	.	CDS	300	350	.	+	.	ID=CDS002;Parent=GENE001.1
chr1	.	CDS	400	450	.	+	.	ID=CDS003;Parent=GENE001.1
chr1	.	gene	600	1000	.	-	.	ID=GENE002
chr1	.	mRNA	600	1000	.	-	.	ID=GENE002.1;Parent=GENE002
chr1	.	CDS	650	700	.	-	.	ID=CDS004;Parent=GENE002.1
chr1	.	CDS	800	850	.	-	.	ID=CDS005;Parent=GENE002.1
"""

    # Create index file
    index_content = """genotype,annotation,assembly
test1,test.gff3,test.fa
"""

    # Create a pseudo FASTA file 
    # fasta_content = ">chr1\n" + "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" * 20
    _fasta_seq = "".join([c*i for i in range(1,25) for c in "ACGT"])
    fasta_content = ">chr1\n" + "\n".join([_fasta_seq[i*80:(i+1)*80] for i in range(int(len(_fasta_seq)/80))])


    # Create target genes file
    target_content = """gene_name,gene_ID_test1
GENE001,GENE001
GENE002,GENE002"""

    # Create BigWig index file
    bigwig_index_content = """genotype,condition,bigwig
test1,test_condition,test.bw
"""

    # Write files
    (temp_test_dir / "test.gff3").write_text(gff_content)
    (temp_test_dir / "index.csv").write_text(index_content)
    (temp_test_dir / "test.fa").write_text(fasta_content)
    (temp_test_dir / "target_genes.csv").write_text(target_content)
    (temp_test_dir / "bigwig_index.csv").write_text(bigwig_index_content)

    # Create BigWig folder
    bigwig_folder = temp_test_dir / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)

    # Create actual BigWig file using pyBigWig
    bigwig_file = bigwig_folder / "test.bw"
    
    # Create chromosome sizes
    chrom_sizes = [("chr1", 2000)]

    # Create BigWig file with some test data
    bw = pyBigWig.open(str(bigwig_file), "w")
    bw.addHeader(chrom_sizes)

    # Add some test values
    # Create signal values for the gene regions
    start = 100
    end = 500
    values = np.linspace(0,5, end - start)  # Test values between 0-5
    bw.addEntries("chr1", start, values=values, span=1, step=1)

    # Add values for the second gene
    start = 600
    end = 1000
    values = np.linspace(5,10, end - start)  # Test values between 5-10
    bw.addEntries("chr1", start, values=values, span=1, step=1)

    bw.close()

    return temp_test_dir
