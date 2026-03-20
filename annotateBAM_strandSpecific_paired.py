import pysam
import sys

# Parameters
# 1: input_bam, 2: te_bed, 3: output_bam, 4: min_overlap_pct (0.0-1.0), 
# 5: mapq_thr, 6: strand_mode (Unstranded, Forward, Reverse)

inputfile = sys.argv[1]
te_bed_file = sys.argv[2]
outfilename = sys.argv[3]
overlap_threshold = float(sys.argv[4])
mapq_thr = int(sys.argv[5])
strand_mode = sys.argv[6] # 'Unstranded', 'Forward', or 'Reverse'

samfile = pysam.AlignmentFile(inputfile, "rb")
te_bed_iterator = pysam.tabix_iterator(open(te_bed_file, "r"), pysam.asBed())
outfile = pysam.AlignmentFile(outfilename, "wb", template=samfile)

for te in te_bed_iterator:
    te_chrom = te[0]
    te_start = int(te[1])
    te_end = int(te[2])
    te_locusname = te[3]
    te_strand = te[5]
    te_name = te_locusname.split("|")[3] if "|" in te_locusname else te_locusname

    # Fetch reads overlapping the TE coordinates
    for sam_record in samfile.fetch(te_chrom, te_start, te_end):
        
        # --- 1. STRANDEDNESS LOGIC (Fragment-aware) ---
        if strand_mode != "Unstranded":
            read_is_rev = sam_record.is_reverse
            
            # Reconcile R1/R2 for Paired-End to find original RNA orientation
            if sam_record.is_paired:
                if sam_record.is_read1:
                    frag_is_rev = read_is_rev
                else:
                    frag_is_rev = not read_is_rev
            else:
                frag_is_rev = read_is_rev

            # Map orientation to Strand
            if strand_mode == "Forward":
                read_strand = "-" if frag_is_rev else "+"
            else: # Reverse (e.g., 10x Genomics)
                read_strand = "+" if frag_is_rev else "-"
            
            if read_strand != te_strand:
                continue

        # --- 2. PROPORTIONAL INTERSECTION (Splicing-aware) ---
        # get_blocks() returns [(start, end)] of mapped segments, skipping "N"
        read_blocks = sam_record.get_blocks()
        total_overlap_bp = 0
        
        for b_start, b_end in read_blocks:
            inter_start = max(te_start, b_start)
            inter_end = min(te_end, b_end)
            if inter_start < inter_end:
                total_overlap_bp += (inter_end - inter_start)

        # Calculate if overlap meets the user-defined percentage of the mapped read
        # Using query_alignment_length (excludes soft-clips)
        if total_overlap_bp < (sam_record.query_alignment_length * overlap_threshold):
            continue

        # --- 3. TAGGING AND WRITING ---
        sam_record.set_tag("GX", te_locusname)
        if sam_record.mapping_quality >= mapq_thr:
            sam_record.set_tag("GN", "SoloTE|" + te_locusname)
        else:
            sam_record.set_tag("GN", "SoloTE|" + te_name)
            
        outfile.write(sam_record)

samfile.close()
outfile.close()