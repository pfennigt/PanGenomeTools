"""
Pytest fixtures for PanGenomeTools tests.
"""

import pytest
from pathlib import Path
import tempfile
import shutil

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
chr1	.	CDS	150	200	.	+	.	ID=CDS001;Parent=GENE001
chr1	.	CDS	300	350	.	+	.	ID=CDS002;Parent=GENE001
chr1	.	CDS	400	450	.	+	.	ID=CDS003;Parent=GENE001
chr1	.	gene	600	1000	.	-	.	ID=GENE002
chr1	.	CDS	650	700	.	-	.	ID=CDS004;Parent=GENE002
chr1	.	CDS	800	850	.	-	.	ID=CDS005;Parent=GENE002
"""

    # Create index file
    index_content = """genotype,annotation,assembly
test1,test.gff,test.fa
"""

    # Create FASTA file
    fasta_content = ">chr1\n" + "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" * 20

    # Write files
    (temp_test_dir / "test.gff").write_text(gff_content)
    (temp_test_dir / "index.csv").write_text(index_content)
    (temp_test_dir / "test.fa").write_text(fasta_content)

    return temp_test_dir
