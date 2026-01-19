"""
Unit tests for BioSeqAnalyzer.

Run with: pytest tests/test_analyzer.py -v
"""

import pytest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analyzer import BioSeqAnalyzer, FileNotFoundError, InvalidFileFormatError


# Path to test data
TEST_FASTA = Path(__file__).parent.parent / "data" / "test.fasta"


@pytest.fixture
def analyzer():
    """Create a BioSeqAnalyzer instance with test data loaded."""
    ana = BioSeqAnalyzer(str(TEST_FASTA))
    ana.load_sequences()
    return ana


class TestGCContent:
    """Tests for GC content calculation."""
    
    def test_gc_content_calculation(self, analyzer):
        """Test GC content on sequence with known 50% GC content.
        
        seq1 = ATGCATGCATGCATGCATGC (20 bp)
        A=5, T=5, G=5, C=5 → GC = 10/20 = 50%
        """
        gc_results = analyzer.calculate_gc_content()
        
        assert "seq1" in gc_results
        assert gc_results["seq1"] == pytest.approx(50.0, abs=0.1)
    
    def test_gc_returns_dict_for_all_sequences(self, analyzer):
        """Test that GC content returns results for all loaded sequences."""
        gc_results = analyzer.calculate_gc_content()
        
        assert len(gc_results) == len(analyzer.sequences)
        assert all(isinstance(v, float) for v in gc_results.values())
        assert all(0 <= v <= 100 for v in gc_results.values())
    
    def test_gc_raises_error_when_no_sequences(self):
        """Test that GC calculation raises ValueError if no sequences loaded."""
        ana = BioSeqAnalyzer(str(TEST_FASTA))
        # Don't call load_sequences()
        
        with pytest.raises(ValueError, match="No sequences loaded"):
            ana.calculate_gc_content()


class TestTranslation:
    """Tests for DNA to protein translation."""
    
    def test_translation_known_sequence(self, analyzer):
        """Test translation of a known DNA sequence.
        
        seq2 = ATGAAATTTGGGTAA (15 bp)
        ATG = M (Met), AAA = K (Lys), TTT = F (Phe), GGG = G (Gly), TAA = * (Stop)
        Expected protein: MKFG*
        """
        translations = analyzer.get_translations(table=1)
        
        assert "seq2" in translations
        assert translations["seq2"] == "MKFG*"
    
    def test_translation_with_different_table(self, analyzer):
        """Test that translation works with different codon tables."""
        translations_std = analyzer.get_translations(table=1)
        translations_bac = analyzer.get_translations(table=11)
        
        # Both should return results
        assert len(translations_std) == len(analyzer.sequences)
        assert len(translations_bac) == len(analyzer.sequences)
    
    def test_translation_returns_dict(self, analyzer):
        """Test that translation returns a dictionary with string values."""
        translations = analyzer.get_translations()
        
        assert isinstance(translations, dict)
        assert all(isinstance(v, str) for v in translations.values())


class TestReverseComplement:
    """Tests for reverse complement generation."""
    
    def test_reverse_complement_known_sequence(self, analyzer):
        """Test reverse complement on a known sequence.
        
        seq3 = AAACCCGGGTTTATGC (16 bp)
        Complement: TTTGGGCCCAAATACG
        Reverse complement: GCATAAACCCGGGTTT
        """
        rc_records = analyzer.reverse_complement_all()
        
        # Find seq3_rc
        seq3_rc = None
        for record in rc_records:
            if record.id == "seq3_rc":
                seq3_rc = record
                break
        
        assert seq3_rc is not None
        assert str(seq3_rc.seq) == "GCATAAACCCGGGTTT"
    
    def test_reverse_complement_returns_seqrecords(self, analyzer):
        """Test that reverse complement returns list of SeqRecord objects."""
        from Bio.SeqRecord import SeqRecord
        
        rc_records = analyzer.reverse_complement_all()
        
        assert isinstance(rc_records, list)
        assert len(rc_records) == len(analyzer.sequences)
        assert all(isinstance(r, SeqRecord) for r in rc_records)
    
    def test_reverse_complement_ids_have_suffix(self, analyzer):
        """Test that reverse complement IDs have '_rc' suffix."""
        rc_records = analyzer.reverse_complement_all()
        
        for record in rc_records:
            assert record.id.endswith("_rc")
    
    def test_reverse_complement_preserves_length(self, analyzer):
        """Test that reverse complement sequences have the same length."""
        rc_records = analyzer.reverse_complement_all()
        
        for orig, rc in zip(analyzer.sequences, rc_records):
            assert len(rc.seq) == len(orig.seq)


class TestFileHandling:
    """Tests for file handling and error cases."""
    
    def test_file_not_found_raises_error(self):
        """Test that missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            BioSeqAnalyzer("nonexistent_file.fasta")
    
    def test_load_sequences_returns_list(self, analyzer):
        """Test that load_sequences returns the loaded sequences."""
        # analyzer fixture already has sequences loaded
        assert len(analyzer.sequences) == 4
    
    def test_analyzer_repr(self, analyzer):
        """Test the string representation of the analyzer."""
        repr_str = repr(analyzer)
        assert "BioSeqAnalyzer" in repr_str
        assert "test.fasta" in repr_str


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
