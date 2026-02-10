"""Tests for raccoon.utils.io module."""
import os
import tempfile
import pytest
from raccoon.utils import io


class TestEnsureOutputDirectory:
    """Tests for ensure_output_directory function."""

    def test_creates_nonexistent_directory(self, tmp_path):
        """Test that function creates a directory if it doesn't exist."""
        new_dir = tmp_path / "new_output"
        assert not new_dir.exists()
        
        result = io.ensure_output_directory(str(new_dir))
        
        assert result is True
        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_accepts_existing_directory(self, tmp_path):
        """Test that function accepts an existing writable directory."""
        result = io.ensure_output_directory(str(tmp_path))
        
        assert result is True

    def test_rejects_readonly_directory(self, tmp_path):
        """Test that function rejects a read-only directory."""
        readonly_dir = tmp_path / "readonly"
        readonly_dir.mkdir()
        os.chmod(str(readonly_dir), 0o444)
        
        try:
            result = io.ensure_output_directory(str(readonly_dir))
            assert result is False
        finally:
            # restore permissions for cleanup
            os.chmod(str(readonly_dir), 0o755)


class TestValidateInputFile:
    """Tests for validate_input_file function."""

    def test_accepts_existing_readable_file(self, tmp_path):
        """Test that function accepts an existing readable file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("content")
        
        result = io.validate_input_file(str(test_file))
        
        assert result is True

    def test_rejects_nonexistent_file(self, tmp_path):
        """Test that function rejects a nonexistent file."""
        result = io.validate_input_file(str(tmp_path / "nonexistent.txt"))
        
        assert result is False

    def test_accepts_empty_filepath(self):
        """Test that function accepts empty/None filepath for optional files."""
        result = io.validate_input_file("")
        
        assert result is True

    def test_rejects_readonly_file(self, tmp_path):
        """Test that function rejects a read-only file."""
        test_file = tmp_path / "readonly.txt"
        test_file.write_text("content")
        os.chmod(str(test_file), 0o000)
        
        try:
            result = io.validate_input_file(str(test_file))
            assert result is False
        finally:
            # restore permissions for cleanup
            os.chmod(str(test_file), 0o644)


class TestValidateAlignmentFile:
    """Tests for validate_alignment_file function."""

    def test_accepts_existing_alignment_file(self, tmp_path):
        """Test that function accepts an existing alignment file."""
        aln_file = tmp_path / "alignment.fasta"
        aln_file.write_text(">seq1\nACGT")
        
        result = io.validate_alignment_file(str(aln_file))
        
        assert result is True

    def test_rejects_nonexistent_alignment_file(self, tmp_path):
        """Test that function rejects a nonexistent alignment file."""
        result = io.validate_alignment_file(str(tmp_path / "nonexistent.fasta"))
        
        assert result is False


class TestValidateGenbankFile:
    """Tests for validate_genbank_file function."""

    def test_accepts_existing_genbank_file(self, tmp_path):
        """Test that function accepts an existing GenBank file."""
        gb_file = tmp_path / "reference.gb"
        gb_file.write_text("LOCUS test")
        
        result = io.validate_genbank_file(str(gb_file))
        
        assert result is True

    def test_accepts_empty_filepath_optional(self):
        """Test that function accepts empty filepath (optional)."""
        result = io.validate_genbank_file("")
        
        assert result is True

    def test_rejects_nonexistent_genbank_file(self, tmp_path):
        """Test that function rejects a nonexistent GenBank file."""
        result = io.validate_genbank_file(str(tmp_path / "nonexistent.gb"))
        
        assert result is False


class TestValidateReferenceFile:
    """Tests for validate_reference_file function."""

    def test_accepts_existing_reference_file(self, tmp_path):
        """Test that function accepts an existing reference file."""
        ref_file = tmp_path / "reference.fasta"
        ref_file.write_text(">ref\nACGT")
        
        result = io.validate_reference_file(str(ref_file))
        
        assert result is True

    def test_accepts_empty_filepath_optional(self):
        """Test that function accepts empty filepath (optional)."""
        result = io.validate_reference_file("")
        
        assert result is True

    def test_rejects_nonexistent_reference_file(self, tmp_path):
        """Test that function rejects a nonexistent reference file."""
        result = io.validate_reference_file(str(tmp_path / "nonexistent.fasta"))
        
        assert result is False
