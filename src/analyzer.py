"""
BioSeqAnalyzer: A class for analyzing biological sequences from FASTA files.
"""

import os
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional
from urllib.error import HTTPError, URLError
import socket
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML


class FileNotFoundError(Exception):
    """Raised when the specified FASTA file is not found."""
    pass


class InvalidFileFormatError(Exception):
    """Raised when the file format is invalid or unreadable."""
    pass


class BioSeqAnalyzer:
    """
    A bioinformatics analyzer for FASTA sequence files.
    
    Attributes:
        file_path (Path): Path to the FASTA file.
        sequences (List[SeqRecord]): List of parsed sequence records.
    
    Example:
        >>> analyzer = BioSeqAnalyzer("sequences.fasta")
        >>> analyzer.load_sequences()
        >>> print(f"Loaded {len(analyzer.sequences)} sequences")
    """
    
    def __init__(self, file_path: str) -> None:
        """
        Initialize the BioSeqAnalyzer with a path to a FASTA file.
        
        Args:
            file_path: Path to the FASTA file to analyze.
        
        Raises:
            FileNotFoundError: If the specified file does not exist.
        """
        self.file_path = Path(file_path)
        self.sequences: List[SeqRecord] = []
        
        # Validate file exists at initialization
        if not self.file_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found: {self.file_path}"
            )
        
        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"Path is not a file: {self.file_path}"
            )
    
    def load_sequences(self) -> List[SeqRecord]:
        """
        Parse the FASTA file and load sequences into memory.
        
        Uses Biopython's SeqIO to parse the FASTA file and stores
        all sequence records in the `sequences` attribute.
        
        Returns:
            List of SeqRecord objects parsed from the file.
        
        Raises:
            InvalidFileFormatError: If the file cannot be parsed as FASTA.
        """
        try:
            self.sequences = list(SeqIO.parse(str(self.file_path), "fasta"))
            
            if not self.sequences:
                raise InvalidFileFormatError(
                    f"No valid FASTA sequences found in: {self.file_path}"
                )
            
            return self.sequences
            
        except Exception as e:
            if isinstance(e, InvalidFileFormatError):
                raise
            raise InvalidFileFormatError(
                f"Error parsing FASTA file '{self.file_path}': {e}"
            ) from e
    
    def calculate_gc_content(self) -> Dict[str, float]:
        """
        Calculate GC content percentage for each loaded sequence.
        
        Uses Biopython's gc_fraction which correctly handles ambiguous
        nucleotides (e.g., S = G or C counts as 1, N = any counts as 0.25).
        
        Returns:
            Dictionary mapping sequence ID to GC percentage (0-100).
        
        Raises:
            ValueError: If no sequences have been loaded.
        
        Example:
            >>> analyzer.load_sequences()
            >>> gc_dict = analyzer.calculate_gc_content()
            >>> for seq_id, gc in gc_dict.items():
            ...     print(f"{seq_id}: {gc:.2f}%")
        """
        if not self.sequences:
            raise ValueError(
                "No sequences loaded. Call load_sequences() first."
            )
        
        gc_content = {}
        for record in self.sequences:
            # gc_fraction handles ambiguous nucleotides correctly
            # Returns fraction (0-1), multiply by 100 for percentage
            gc_content[record.id] = gc_fraction(record.seq) * 100
        
        return gc_content
    
    def get_translations(self, table: int = 1) -> Dict[str, str]:
        """
        Translate nucleotide sequences to protein sequences.
        
        Uses Biopython's translate method with the specified genetic code table.
        Handles sequences that are not multiples of 3 by trimming extra nucleotides.
        Ambiguous codons that cannot be translated will result in 'X' in the output.
        
        Args:
            table: NCBI genetic code table number (default=1, Standard Code).
                   Common tables:
                   - 1: Standard Code
                   - 2: Vertebrate Mitochondrial
                   - 11: Bacterial/Archaeal/Plant Plastid
        
        Returns:
            Dictionary mapping sequence ID to translated protein string.
        
        Raises:
            ValueError: If no sequences have been loaded.
        
        Example:
            >>> translations = analyzer.get_translations(table=1)
            >>> for seq_id, protein in translations.items():
            ...     print(f"{seq_id}: {protein[:50]}...")
        """
        if not self.sequences:
            raise ValueError(
                "No sequences loaded. Call load_sequences() first."
            )
        
        translations = {}
        for record in self.sequences:
            seq = record.seq
            
            # Trim sequence to multiple of 3 for translation
            remainder = len(seq) % 3
            if remainder:
                seq = seq[:-remainder]
            
            # Translate with specified table
            # to_stop=False to get full translation including stop codons (*)
            protein = str(seq.translate(table=table))
            translations[record.id] = protein
        
        return translations
    
    def reverse_complement_all(self) -> List[SeqRecord]:
        """
        Generate reverse complement for all loaded sequences.
        
        Creates new SeqRecord objects containing the reverse complement
        of each original sequence. Correctly handles ambiguous nucleotides
        (e.g., R -> Y, M -> K, etc.) using IUPAC conventions.
        
        Returns:
            List of new SeqRecord objects with reverse complemented sequences.
            IDs are appended with '_rc' suffix.
        
        Raises:
            ValueError: If no sequences have been loaded.
        
        Example:
            >>> rc_seqs = analyzer.reverse_complement_all()
            >>> for rc in rc_seqs:
            ...     print(f"{rc.id}: {rc.seq[:50]}...")
        """
        if not self.sequences:
            raise ValueError(
                "No sequences loaded. Call load_sequences() first."
            )
        
        rc_records = []
        for record in self.sequences:
            # Create reverse complement sequence
            rc_seq = record.seq.reverse_complement()
            
            # Create new SeqRecord with updated metadata
            rc_record = SeqRecord(
                seq=rc_seq,
                id=f"{record.id}_rc",
                name=f"{record.name}_rc" if record.name else "",
                description=f"Reverse complement of {record.description}"
            )
            rc_records.append(rc_record)
        
        return rc_records
    
    def run_blast_search(
        self,
        sequence_id: str,
        database: str = "nr",
        program: str = "blastn",
        hitlist_size: int = 5,
        max_retries: int = 3,
        timeout: int = 120
    ) -> List[Dict[str, str]]:
        """
        Perform a remote BLAST search using NCBI's QBLAST service.
        
        Extracts the specified sequence by ID and submits it to NCBI BLAST.
        Parses and returns the top hits with their accession IDs and descriptions.
        
        Args:
            sequence_id: ID of the sequence to search (must be loaded).
            database: NCBI database to search (default='nr').
                      Common options: 'nr', 'nt', 'refseq_rna', 'swissprot'
            program: BLAST program to use (default='blastn').
                     Options: 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'
            hitlist_size: Maximum number of hits to return (default=5).
            max_retries: Number of retry attempts on network failure (default=3).
            timeout: Timeout in seconds for the BLAST request (default=120).
        
        Returns:
            List of dictionaries, each containing:
            - 'accession': Hit accession ID
            - 'description': Hit description/title
            - 'e_value': E-value of the alignment
            - 'score': Bit score of the alignment
        
        Raises:
            ValueError: If sequence_id is not found or no sequences loaded.
            ConnectionError: If BLAST request fails after all retries.
        
        Example:
            >>> results = analyzer.run_blast_search('seq1', database='nt')
            >>> for hit in results:
            ...     print(f"{hit['accession']}: {hit['description']}")
        """
        if not self.sequences:
            raise ValueError(
                "No sequences loaded. Call load_sequences() first."
            )
        
        # Find the sequence by ID
        target_seq = None
        for record in self.sequences:
            if record.id == sequence_id:
                target_seq = record
                break
        
        if target_seq is None:
            available_ids = [r.id for r in self.sequences]
            raise ValueError(
                f"Sequence ID '{sequence_id}' not found. "
                f"Available IDs: {available_ids}"
            )
        
        # Perform BLAST with retry logic for network errors
        result_handle = None
        last_error = None
        
        for attempt in range(max_retries):
            try:
                print(f"Submitting BLAST search (attempt {attempt + 1}/{max_retries})...")
                print(f"  Sequence: {sequence_id} ({len(target_seq.seq)} bp)")
                print(f"  Database: {database}, Program: {program}")
                
                result_handle = NCBIWWW.qblast(
                    program=program,
                    database=database,
                    sequence=str(target_seq.seq),
                    hitlist_size=hitlist_size,
                    expect=10.0,  # E-value threshold
                )
                print("BLAST search completed successfully.")
                break
                
            except (HTTPError, URLError, socket.timeout, socket.error) as e:
                last_error = e
                print(f"Network error on attempt {attempt + 1}: {e}")
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 5  # Exponential backoff
                    print(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
            except Exception as e:
                # Re-raise unexpected errors
                raise ConnectionError(f"Unexpected BLAST error: {e}") from e
        
        if result_handle is None:
            raise ConnectionError(
                f"BLAST search failed after {max_retries} attempts. "
                f"Last error: {last_error}"
            )
        
        # Parse BLAST results
        try:
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
        except Exception as e:
            raise ConnectionError(f"Failed to parse BLAST results: {e}") from e
        finally:
            result_handle.close()
        
        # Extract top hits
        hits = []
        for alignment in blast_record.alignments[:hitlist_size]:
            # Get the best HSP (High-scoring Segment Pair)
            best_hsp = alignment.hsps[0] if alignment.hsps else None
            
            hit_info = {
                'accession': alignment.accession,
                'description': alignment.hit_def,
                'e_value': str(best_hsp.expect) if best_hsp else 'N/A',
                'score': str(best_hsp.bits) if best_hsp else 'N/A',
            }
            hits.append(hit_info)
            
            # Print hit information
            print(f"\n  Hit: {hit_info['accession']}")
            print(f"  Description: {hit_info['description'][:80]}...")
            print(f"  E-value: {hit_info['e_value']}, Score: {hit_info['score']}")
        
        if not hits:
            print("No significant hits found.")
        else:
            print(f"\nFound {len(hits)} hit(s).")
        
        return hits
    
    def __len__(self) -> int:
        """Return the number of loaded sequences."""
        return len(self.sequences)
    
    def __repr__(self) -> str:
        """Return a string representation of the analyzer."""
        return (
            f"BioSeqAnalyzer(file='{self.file_path}', "
            f"sequences_loaded={len(self.sequences)})"
        )
