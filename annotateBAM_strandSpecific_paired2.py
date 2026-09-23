import pysam
import sys
from intervaltree import IntervalTree

# Parameters
# 1: input_bam, 2: te_bed, 3: output_bam, 4: min_overlap_pct (0.0-1.0), 
# 5: mapq_thr, 6: strand_mode (Unstranded, Forward, Reverse)
inputfile = sys.argv[1]
te_bed_file = sys.argv[2]
outfilename = sys.argv[3]
overlap_threshold = float(sys.argv[4])
mapq_thr = int(sys.argv[5])
strand_mode = sys.argv[6]

# 1. Build IntervalTrees for fast TE lookup
# We store TEs in a dict: {chrom: IntervalTree(start, end, [name, strand, locus])}
te_trees = {}
with open(te_bed_file, "r") as f:
    for line in f:
        row = line.strip().split("\t")
        chrom, start, end, locus, _, strand = row[0], int(row[1]), int(row[2]), row[3], row[4], row[5]
        if chrom not in te_trees:
            te_trees[chrom] = IntervalTree()
        # Data stores both the specific locus and the general name
        te_name = locus.split("|")[3] if "|" in locus else locus
        te_trees[chrom].addi(start, end, {"locus": locus, "name": te_name, "strand": strand})

samfile = pysam.AlignmentFile(inputfile, "rb")
outfile = pysam.AlignmentFile(outfilename, "wb", template=samfile)

# Dictionary to cache pair information: {query_name: [read1_annotation, read2_annotation]}
pair_buffer = {}

def get_te_info(sam_record):
    """Returns (locus, name) if overlap meets threshold, else (None, None)."""
    if sam_record.reference_name not in te_trees:
        return None, None
    
    # Strand Logic
    if strand_mode != "Unstranded":
        read_is_rev = sam_record.is_reverse
        if sam_record.is_paired:
            frag_is_rev = read_is_rev if sam_record.is_read1 else not read_is_rev
        else:
            frag_is_rev = read_is_rev
        
        # Determine the strand of the original RNA molecule
        if strand_mode == "Reverse":
            # 10x/dUTP: If fragment is reverse-mapped, it's a (+) strand RNA
            actual_strand = "+" if frag_is_rev else "-"
        else:
            # Forward: If fragment is reverse-mapped, it's a (-) strand RNA
            actual_strand = "-" if frag_is_rev else "+"
    else:
        actual_strand = None

    best_locus, best_name, max_overlap = None, None, 0
    potential_tes = te_trees[sam_record.reference_name].overlap(sam_record.reference_start, sam_record.reference_end)
    read_blocks = sam_record.get_blocks()
    
    for te in potential_tes:
        if actual_strand and te.data['strand'] != actual_strand:
            continue
            
        overlap_bp = 0
        for b_start, b_end in read_blocks:
            inter_start = max(te.begin, b_start)
            inter_end = min(te.end, b_end)
            if inter_start < inter_end:
                overlap_bp += (inter_end - inter_start)
        
        if overlap_bp > max_overlap:
            max_overlap = overlap_bp
            best_locus = te.data['locus']
            best_name = te.data['name']

    if max_overlap >= (sam_record.query_alignment_length * overlap_threshold):
        return best_locus, best_name
    return None, None


# Main Loop
for read in samfile.fetch(until_eof=True):
    if not (read.has_tag("CB") and read.has_tag("UB")):
        outfile.write(read)
        continue

    # 1. Individual Read Annotation
    locus, name = get_te_info(read)
    
    # 2. MAPQ Determination for this specific read
    is_unique = read.mapping_quality >= mapq_thr
    read_assignment = locus if is_unique else name

    # 3. Paired-End Reconciliation
    if read.is_paired:
        qname = read.query_name
        if qname not in pair_buffer:
            pair_buffer[qname] = {"read": read, "assign": read_assignment, "chrom": read.reference_name}
            continue
        else:
            mate_data = pair_buffer.pop(qname)
            mate = mate_data["read"]
            
            # --- CHIMERIC CHECK ---
            # If reads map to different chromosomes, we mark as Ambiguous
            # if read.reference_name != mate_data["chrom"]:
            #     final_assignment = "Ambiguous"
            # --- CONSISTENCY CHECK ---
            if read_assignment == mate_data["assign"]:
                final_assignment = read_assignment
            elif read_assignment and mate_data["assign"]:
                # Both hit TEs, but different ones
                final_assignment = "AmbiguousTE"
            else:
                # One hit a TE, the other was None
                final_assignment = read_assignment or mate_data["assign"]

            # --- TAGGING BOTH READS ---
            for r in [read, mate]:
                existing_gn = r.get_tag("GN") if r.has_tag("GN") else None
                if final_assignment:
                    te_tag = "SoloTE|" + final_assignment
                    new_tag = f"{existing_gn};{te_tag}" if (existing_gn and existing_gn != "-") else te_tag
                    r.set_tag("GN", new_tag)
                outfile.write(r)
    else:
        # Single-End Logic
        if read_assignment:
            existing_gn = read.get_tag("GN") if read.has_tag("GN") else None
            te_tag = "SoloTE|" + read_assignment
            new_tag = f"{existing_gn};{te_tag}" if (existing_gn and existing_gn != "-") else te_tag
            read.set_tag("GN", new_tag)
        outfile.write(read)

samfile.close()
outfile.close()