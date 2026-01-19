#!/usr/bin/env python3
"""
Biopython-SeqAnalyzer CLI

A command-line interface for analyzing biological sequences from FASTA files.

Usage:
    python run_analyzer.py sequences.fasta --gc
    python run_analyzer.py sequences.fasta --translate
    python run_analyzer.py sequences.fasta --reverse-complement
    python run_analyzer.py sequences.fasta --blast SEQ_ID --database nt
"""

import argparse
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from src.analyzer import BioSeqAnalyzer, FileNotFoundError, InvalidFileFormatError


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog="run_analyzer.py",
        description="Analyze biological sequences from FASTA files using Biopython.",
        epilog="""
Examples:
  %(prog)s sequences.fasta --gc
      Calculate GC content for all sequences
  
  %(prog)s sequences.fasta --translate --table 11
      Translate sequences using bacterial codon table
  
  %(prog)s sequences.fasta --reverse-complement --output rc_seqs.fasta
      Generate reverse complements and save to file
  
  %(prog)s sequences.fasta --blast SEQ001 --database nt
      Run BLAST search for sequence SEQ001 against nucleotide database
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required positional argument
    parser.add_argument(
        "fasta_file",
        type=str,
        help="Path to the input FASTA file"
    )
    
    # Analysis action flags (mutually exclusive group)
    actions = parser.add_argument_group("Analysis Actions (choose one)")
    action_group = actions.add_mutually_exclusive_group(required=True)
    
    action_group.add_argument(
        "--gc",
        action="store_true",
        help="Calculate GC content percentage for each sequence"
    )
    
    action_group.add_argument(
        "--translate",
        action="store_true",
        help="Translate nucleotide sequences to protein sequences"
    )
    
    action_group.add_argument(
        "--reverse-complement", "--rc",
        action="store_true",
        dest="reverse_complement",
        help="Generate reverse complement of all sequences"
    )
    
    action_group.add_argument(
        "--blast",
        type=str,
        metavar="SEQ_ID",
        help="Run BLAST search for the specified sequence ID"
    )
    
    action_group.add_argument(
        "--info",
        action="store_true",
        help="Display basic information about all sequences"
    )
    
    # Optional parameters
    options = parser.add_argument_group("Optional Parameters")
    
    options.add_argument(
        "--table", "-t",
        type=int,
        default=1,
        help="Genetic code table for translation (default: 1 = Standard). "
             "Common: 1=Standard, 2=Vertebrate Mito, 11=Bacterial"
    )
    
    options.add_argument(
        "--database", "-d",
        type=str,
        default="nr",
        help="NCBI database for BLAST search (default: nr). "
             "Options: nr, nt, refseq_rna, swissprot"
    )
    
    options.add_argument(
        "--program", "-p",
        type=str,
        default="blastn",
        choices=["blastn", "blastp", "blastx", "tblastn", "tblastx"],
        help="BLAST program to use (default: blastn)"
    )
    
    options.add_argument(
        "--hits", "-n",
        type=int,
        default=5,
        help="Number of BLAST hits to return (default: 5)"
    )
    
    options.add_argument(
        "--output", "-o",
        type=str,
        help="Output file path (for reverse-complement action)"
    )
    
    return parser


def print_header(title: str) -> None:
    """Print a formatted section header."""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}\n")


def run_gc_analysis(analyzer: BioSeqAnalyzer) -> None:
    """Run GC content analysis and display results."""
    print_header("GC Content Analysis")
    
    gc_results = analyzer.calculate_gc_content()
    
    print(f"{'Sequence ID':<30} {'GC %':>10}")
    print("-" * 42)
    
    for seq_id, gc_pct in gc_results.items():
        print(f"{seq_id:<30} {gc_pct:>10.2f}%")
    
    # Summary statistics
    if gc_results:
        avg_gc = sum(gc_results.values()) / len(gc_results)
        print("-" * 42)
        print(f"{'Average':<30} {avg_gc:>10.2f}%")
        print(f"{'Min':<30} {min(gc_results.values()):>10.2f}%")
        print(f"{'Max':<30} {max(gc_results.values()):>10.2f}%")


def run_translation(analyzer: BioSeqAnalyzer, table: int) -> None:
    """Run translation and display results."""
    print_header(f"Translation (Genetic Code Table {table})")
    
    translations = analyzer.get_translations(table=table)
    
    for seq_id, protein in translations.items():
        print(f"\n> {seq_id}")
        # Wrap protein sequence at 60 characters
        for i in range(0, len(protein), 60):
            print(f"  {protein[i:i+60]}")
        print(f"  Length: {len(protein)} aa")


def run_reverse_complement(analyzer: BioSeqAnalyzer, output_file: str = None) -> None:
    """Generate reverse complements and optionally save to file."""
    print_header("Reverse Complement")
    
    rc_records = analyzer.reverse_complement_all()
    
    for record in rc_records:
        print(f"\n> {record.id}")
        seq_str = str(record.seq)
        # Wrap sequence at 60 characters
        for i in range(0, len(seq_str), 60):
            print(f"  {seq_str[i:i+60]}")
        print(f"  Length: {len(record.seq)} bp")
    
    # Save to file if requested
    if output_file:
        from Bio import SeqIO
        SeqIO.write(rc_records, output_file, "fasta")
        print(f"\n✓ Saved {len(rc_records)} sequences to: {output_file}")


def run_blast_search(analyzer: BioSeqAnalyzer, seq_id: str, 
                     database: str, program: str, hits: int) -> None:
    """Run BLAST search for a specific sequence."""
    print_header(f"BLAST Search: {seq_id}")
    
    print(f"Database: {database}")
    print(f"Program: {program}")
    print(f"Max hits: {hits}")
    print()
    
    try:
        results = analyzer.run_blast_search(
            sequence_id=seq_id,
            database=database,
            program=program,
            hitlist_size=hits
        )
        
        if results:
            print("\n" + "-" * 60)
            print("Summary of Results:")
            print("-" * 60)
            for i, hit in enumerate(results, 1):
                print(f"\n{i}. {hit['accession']}")
                print(f"   {hit['description'][:70]}...")
                print(f"   E-value: {hit['e_value']}, Score: {hit['score']}")
                
    except ConnectionError as e:
        print(f"\n✗ BLAST search failed: {e}")
        sys.exit(1)


def run_info(analyzer: BioSeqAnalyzer) -> None:
    """Display basic sequence information."""
    print_header("Sequence Information")
    
    print(f"File: {analyzer.file_path}")
    print(f"Total sequences: {len(analyzer.sequences)}")
    print()
    
    print(f"{'Sequence ID':<30} {'Length (bp)':>12} {'Description'}")
    print("-" * 80)
    
    total_bp = 0
    for record in analyzer.sequences:
        length = len(record.seq)
        total_bp += length
        desc = record.description[:30] + "..." if len(record.description) > 30 else record.description
        print(f"{record.id:<30} {length:>12,} {desc}")
    
    print("-" * 80)
    print(f"{'Total':<30} {total_bp:>12,} bp")


def main() -> int:
    """Main entry point for the CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Initialize analyzer
    try:
        analyzer = BioSeqAnalyzer(args.fasta_file)
        analyzer.load_sequences()
        print(f"✓ Loaded {len(analyzer.sequences)} sequence(s) from {args.fasta_file}")
    except FileNotFoundError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        return 1
    except InvalidFileFormatError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        return 1
    
    # Execute the requested action
    try:
        if args.gc:
            run_gc_analysis(analyzer)
        elif args.translate:
            run_translation(analyzer, args.table)
        elif args.reverse_complement:
            run_reverse_complement(analyzer, args.output)
        elif args.blast:
            run_blast_search(analyzer, args.blast, args.database, args.program, args.hits)
        elif args.info:
            run_info(analyzer)
            
    except ValueError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
        return 130
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
