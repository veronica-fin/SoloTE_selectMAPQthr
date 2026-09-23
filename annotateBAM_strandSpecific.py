import pysam
import sys

# Argument Mapping
# 1: input_bam, 2: te_bed, 3: output_bam, 4: min_overlap_pct (0.0 to 1.0), 
# 5: mapq_thr, 6: strand_mode (Unstranded, Forward, Reverse)

inputfile = sys.argv[1]
samfile = pysam.AlignmentFile(inputfile, "rb")

te_bed_iterator = pysam.tabix_iterator(open(sys.argv[2], "r"), pysam.asBed())

outfile = pysam.AlignmentFile(sys.argv[3], "wb", template=samfile)

overlap_threshold = float(sys.argv[4]) # e.g., 0.5 for 50% overlap
thr = int(sys.argv[5])

strand_mode = sys.argv[6]

for te in te_bed_iterator:
    te_chrom, te_start, te_end = te[0], int(te[1]), int(te[2])
    te_locusname, te_strand = te[3], te[5]
    te_name = te_locusname.split("|")[3] if "|" in te_locusname else te_locusname

    for sam_record in samfile.fetch(te_chrom, te_start, te_end):
        
        # --- 1. STRANDEDNESS CHECK ---
        if strand_mode != "Unstranded":
            read_is_rev = sam_record.is_reverse
            # Logic: 10x is 'Reverse' (Read 2 is antisense to RNA)
            if strand_mode == "Forward":
                read_strand = "-" if read_is_rev else "+"
            else: # Reverse
                read_strand = "+" if read_is_rev else "-"
            
            if read_strand != te_strand:
                continue

        # --- 2. PROPORTIONAL INTERSECTION ---
        # Get exact mapped coordinates (handles CIGAR N/D/M automatically)
        read_blocks = sam_record.get_blocks() 
        total_overlap_bp = 0
        
        for b_start, b_end in read_blocks:
            # Calculate intersection of this CIGAR block and the TE
            inter_start = max(te_start, b_start)
            inter_end = min(te_end, b_end)
            if inter_start < inter_end:
                total_overlap_bp += (inter_end - inter_start)

        # Stellarscope-style threshold: (Overlap / Read Length)
        # Using sam_record.query_alignment_length (length of mapped portion)
        if total_overlap_bp < (sam_record.query_alignment_length * overlap_threshold):
            continue

        # --- 3. TAGGING ---
        sam_record.set_tag("GX", te_locusname)
        if sam_record.mapping_quality >= thr:
            sam_record.set_tag("GN", "SoloTE|" + te_locusname)
        else:
            sam_record.set_tag("GN", "SoloTE|" + te_name)
            
        outfile.write(sam_record)

samfile.close()
outfile.close()


# # --- IMPROVED STRAND CHECK FOR PAIRED-END ---
# if strand_mode != "Unstranded":
#     # 1. Get the raw orientation
#     read_is_rev = sam_record.is_reverse
    
#     # 2. If paired, we must reconcile R1 and R2 to find the Fragment Strand
#     if sam_record.is_paired:
#         if sam_record.is_read1:
#             # For most kits (Reverse/10x), R1 is antisense to the original RNA
#             frag_is_rev = read_is_rev
#         else:
#             # R2 is the opposite of R1
#             frag_is_rev = not read_is_rev
#     else:
#         frag_is_rev = read_is_rev

#     # 3. Apply the Forward/Reverse logic to the Fragment orientation
#     if strand_mode == "Forward":
#         read_strand = "-" if frag_is_rev else "+"
#     else: # Reverse (10x / dUTP)
#         read_strand = "+" if frag_is_rev else "-"
    
#     if read_strand != te_strand:
#         continue