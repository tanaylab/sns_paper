#!/usr/bin/env python
"""
Create a perturbed genome FASTA by applying splice/sequence perturbations.

This script reads a reference genome FASTA and perturbation specifications,
applies the perturbations, and writes a new FASTA file that can be used
for training with perturbed sequences.

Usage:
    python scripts/create_perturbed_genome.py \
        --genome_fasta data/mm10.fa \
        --splice_pert data/splice_pert.tsv \
        --output data/mm10_perturbed.fa

Then use in training config:
    data:
      genome_fasta: "data/mm10_perturbed.fa"
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Union

import pysam

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.perturbation_simple import SplicePert, SeqPert, load_splice_pert, load_seq_pert


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a perturbed genome FASTA by applying perturbations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Apply splice perturbations only
    python scripts/create_perturbed_genome.py \\
        --genome_fasta data/mm10.fa \\
        --splice_pert data/splice_pert.tsv \\
        --output data/mm10_perturbed.fa

    # Apply both splice and sequence perturbations
    python scripts/create_perturbed_genome.py \\
        --genome_fasta data/mm10.fa \\
        --splice_pert data/splice_pert.tsv \\
        --seq_pert data/seq_pert.tsv \\
        --output data/mm10_perturbed.fa

    # Process specific chromosomes only
    python scripts/create_perturbed_genome.py \\
        --genome_fasta data/mm10.fa \\
        --splice_pert data/splice_pert.tsv \\
        --output data/mm10_perturbed.fa \\
        --chromosomes chr1,chr2,chr3
        """
    )

    parser.add_argument(
        '--genome_fasta', type=str, required=True,
        help='Input reference genome FASTA file'
    )
    parser.add_argument(
        '--output', type=str, required=True,
        help='Output perturbed genome FASTA file'
    )
    parser.add_argument(
        '--splice_pert', type=str, default=None,
        help='Path to splice perturbation TSV file'
    )
    parser.add_argument(
        '--seq_pert', type=str, default=None,
        help='Path to sequence perturbation TSV file'
    )
    parser.add_argument(
        '--chromosomes', type=str, default=None,
        help='Comma-separated list of chromosomes to include (default: all)'
    )
    parser.add_argument(
        '--line_width', type=int, default=60,
        help='Line width for FASTA output (default: 60)'
    )

    return parser.parse_args()


def index_perturbations_by_chrom(
    splice_perts: List[SplicePert],
    seq_perts: List[SeqPert]
) -> tuple:
    """Index perturbations by chromosome for efficient lookup."""
    splice_by_chrom: Dict[str, List[SplicePert]] = defaultdict(list)
    seq_by_chrom: Dict[str, List[SeqPert]] = defaultdict(list)

    for pert in splice_perts:
        splice_by_chrom[pert.chrom].append(pert)

    for pert in seq_perts:
        seq_by_chrom[pert.chrom].append(pert)

    return splice_by_chrom, seq_by_chrom


def apply_perturbations_to_chrom(
    seq: str,
    splice_perts: List[SplicePert],
    seq_perts: List[SeqPert],
    genome: pysam.FastaFile
) -> tuple:
    """
    Apply perturbations to a chromosome sequence.

    Perturbations are applied from end to start (by position descending)
    to avoid coordinate invalidation issues.

    Args:
        seq: Original chromosome sequence
        splice_perts: Splice perturbations for this chromosome
        seq_perts: Sequence perturbations for this chromosome
        genome: Genome FASTA handle for fetching donor sequences

    Returns:
        Tuple of (perturbed_sequence, num_perturbations_applied, bases_modified)
    """
    # Convert to bytearray for efficient in-place modification
    seq_arr = bytearray(seq.upper().encode('ascii'))

    # Combine all perturbations
    all_perts: List[Union[SplicePert, SeqPert]] = list(splice_perts) + list(seq_perts)

    # Sort by start position descending (apply from end to preserve coordinates)
    all_perts.sort(key=lambda p: p.start, reverse=True)

    num_applied = 0
    bases_modified = 0

    for pert in all_perts:
        # Get donor sequence
        if pert.seq is not None:
            donor = pert.seq.upper()
        else:
            try:
                donor = genome.fetch(pert.chrom_src, pert.start_src, pert.end_src).upper()
            except Exception as e:
                print(f"  Warning: Could not fetch donor for {pert.id}: {e}")
                continue

        # Validate coordinates
        if pert.start < 0 or pert.end > len(seq_arr):
            print(f"  Warning: Perturbation {pert.id} out of bounds "
                  f"(start={pert.start}, end={pert.end}, chrom_len={len(seq_arr)})")
            continue

        # Apply perturbation
        target_len = pert.end - pert.start
        donor_bytes = donor.encode('ascii')

        if len(donor_bytes) != target_len:
            print(f"  Warning: Donor length ({len(donor_bytes)}) != target length ({target_len}) "
                  f"for {pert.id}. Perturbation will change chromosome length.")

        seq_arr[pert.start:pert.end] = donor_bytes
        num_applied += 1
        bases_modified += target_len

    return seq_arr.decode('ascii'), num_applied, bases_modified


def write_fasta_sequence(f, chrom: str, seq: str, line_width: int = 60):
    """Write a single chromosome to FASTA format."""
    f.write(f">{chrom}\n")
    for i in range(0, len(seq), line_width):
        f.write(seq[i:i+line_width] + "\n")


def main():
    args = parse_args()

    # Validate inputs
    if not os.path.exists(args.genome_fasta):
        print(f"Error: Genome FASTA not found: {args.genome_fasta}")
        sys.exit(1)

    if args.splice_pert is None and args.seq_pert is None:
        print("Error: Must specify at least one of --splice_pert or --seq_pert")
        sys.exit(1)

    # Load perturbations
    splice_perts: List[SplicePert] = []
    seq_perts: List[SeqPert] = []

    if args.splice_pert:
        if not os.path.exists(args.splice_pert):
            print(f"Error: Splice perturbation file not found: {args.splice_pert}")
            sys.exit(1)
        print(f"Loading splice perturbations from {args.splice_pert}...")
        splice_perts = load_splice_pert(args.splice_pert)
        print(f"  Loaded {len(splice_perts)} splice perturbations")

    if args.seq_pert:
        if not os.path.exists(args.seq_pert):
            print(f"Error: Sequence perturbation file not found: {args.seq_pert}")
            sys.exit(1)
        print(f"Loading sequence perturbations from {args.seq_pert}...")
        seq_perts = load_seq_pert(args.seq_pert)
        print(f"  Loaded {len(seq_perts)} sequence perturbations")

    # Index perturbations by chromosome
    splice_by_chrom, seq_by_chrom = index_perturbations_by_chrom(splice_perts, seq_perts)

    # Get chromosomes with perturbations
    perturbed_chroms = set(splice_by_chrom.keys()) | set(seq_by_chrom.keys())
    print(f"\nPerturbations affect {len(perturbed_chroms)} chromosomes: "
          f"{', '.join(sorted(perturbed_chroms)[:5])}{'...' if len(perturbed_chroms) > 5 else ''}")

    # Open genome FASTA
    print(f"\nOpening genome FASTA: {args.genome_fasta}")
    genome = pysam.FastaFile(args.genome_fasta)

    # Determine chromosomes to process
    all_chroms = list(genome.references)
    if args.chromosomes:
        requested_chroms = [c.strip() for c in args.chromosomes.split(',')]
        chroms_to_process = [c for c in requested_chroms if c in all_chroms]
        missing = set(requested_chroms) - set(chroms_to_process)
        if missing:
            print(f"Warning: Requested chromosomes not found in genome: {missing}")
    else:
        chroms_to_process = all_chroms

    print(f"Processing {len(chroms_to_process)} chromosomes")

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Process chromosomes and write output
    total_perturbations = 0
    total_bases_modified = 0

    print(f"\nWriting perturbed genome to {args.output}...")
    with open(args.output, 'w') as f:
        for chrom in chroms_to_process:
            # Fetch original sequence
            seq = genome.fetch(chrom)

            # Get perturbations for this chromosome
            chrom_splice = splice_by_chrom.get(chrom, [])
            chrom_seq = seq_by_chrom.get(chrom, [])

            num_perts = len(chrom_splice) + len(chrom_seq)

            if num_perts > 0:
                # Apply perturbations
                perturbed_seq, num_applied, bases_modified = apply_perturbations_to_chrom(
                    seq, chrom_splice, chrom_seq, genome
                )
                print(f"  {chrom}: {num_applied} perturbations applied, "
                      f"{bases_modified:,} bases modified")
                total_perturbations += num_applied
                total_bases_modified += bases_modified
            else:
                # No perturbations - use original sequence
                perturbed_seq = seq.upper()

            # Write to output
            write_fasta_sequence(f, chrom, perturbed_seq, args.line_width)

    genome.close()

    # Create FASTA index
    print(f"\nCreating FASTA index...")
    pysam.faidx(args.output)

    # Print summary
    output_size = output_path.stat().st_size
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Output file: {args.output}")
    print(f"  Output size: {output_size / 1e9:.2f} GB")
    print(f"  Chromosomes processed: {len(chroms_to_process)}")
    print(f"  Total perturbations applied: {total_perturbations}")
    print(f"  Total bases modified: {total_bases_modified:,}")
    print(f"  Index file: {args.output}.fai")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
