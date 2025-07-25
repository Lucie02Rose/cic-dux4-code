import sys
import pysam
from collections import defaultdict

if len(sys.argv) != 4:
    print(f"Usage: python {sys.argv[0]} <contig_bam> <input_regions> <output_mapped_coords>")
    sys.exit(1)

contig_bam = sys.argv[1]
input_regions = sys.argv[2]
output_coords = sys.argv[3]
MAPQ_THRESHOLD = 20

bam = pysam.AlignmentFile(contig_bam, "rb")
contig_alignments = defaultdict(list)

# Step 1: Read all alignments into memory per contig
for aln in bam.fetch(until_eof=True):
    if aln.is_unmapped or aln.mapping_quality < MAPQ_THRESHOLD:
        continue
    if not aln.is_supplementary and not aln.is_secondary:
        aln_type = "primary"
    elif aln.is_supplementary and aln.mapping_quality >=50:
        aln_type = "supplementary"
    else:
        continue

    contig = aln.query_name
    ref = bam.get_reference_name(aln.reference_id)
    strand = '-' if aln.is_reverse else '+'
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    query_start = aln.query_alignment_start
    query_end = aln.query_alignment_end

    contig_alignments[contig].append({
        "ref": ref,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "strand": strand,
        "query_start": query_start,
        "query_end": query_end,
        "mapq": aln.mapping_quality,
        "aln_len": query_end - query_start,
        "type": aln_type
    })

bam.close()

# Step 2: Map input regions using the best alignment
# Step 2: Map input regions using *all alignments* that sufficiently overlap the D4Z4 region
with open(input_regions) as f_in, open(output_coords, "w") as f_out:
    for line in f_in:
        line = line.strip()
        if not line:
            continue

        contig, coords = line.split(":")
        start, end = map(int, coords.split("-"))
        d4z4_len = end - start

        if contig not in contig_alignments:
            print(f"Warning: {contig} not found in BAM alignments.")
            continue

        found = False
        for aln in contig_alignments[contig]:
            # Compute overlap between this alignment and the D4Z4 region
            overlap_start = max(aln["query_start"], start)
            overlap_end = min(aln["query_end"], end)
            overlap_len = max(0, overlap_end - overlap_start)

            if overlap_len / d4z4_len >= 0.5:  # ≥50% of D4Z4 region overlaps
                found = True

                ref = aln["ref"]
                strand = aln["strand"]
                ref_start = aln["ref_start"]
                ref_end = aln["ref_end"]
                q_start = aln["query_start"]
                mapq = aln["mapq"]

                offset = start - q_start

                if strand == '+':
                    ref_d4z4_start = ref_start + offset
                else:
                    ref_d4z4_start = ref_end - offset - d4z4_len

                ref_d4z4_end = ref_d4z4_start + d4z4_len

                # Output each valid overlapping alignment
                f_out.write(f"{contig}:{start}-{end}\t{ref}:{ref_d4z4_start+1}-{ref_d4z4_end}\t{strand}\tMAPQ={mapq}\n")

        if not found:
            print(f"Warning: No alignment overlaps D4Z4 region in {contig}:{start}-{end}")
