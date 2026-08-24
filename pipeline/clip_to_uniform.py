#!/usr/bin/env python3
"""
Quick helper: read already-trimmed FASTAs from 04_trimmed/,
clip all seqs to the mode length (from 3' end), drop any too short,
write uniform-length files back, and report drops per gene.
"""

import glob
import os
from collections import Counter

IN_DIR = "data/04_trimmed"
NOTES_DIR = os.path.join(IN_DIR, "notes")

os.makedirs(NOTES_DIR, exist_ok=True)
log_path = os.path.join(NOTES_DIR, "uniform_clip_log.tsv")
log = open(log_path, "w")
log.write("gene\ttotal\tkept\tdropped\tdrop_pct\ttoo_short\tclipped_3prime\tat_mode\tuniform_len\n")

for fasta in sorted(glob.glob(os.path.join(IN_DIR, "*.fasta"))):
    gene = os.path.basename(fasta).replace(".fasta", "")

    # Read all sequences
    entries = []
    header, parts = None, []
    with open(fasta) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
        if header is not None:
            entries.append((header, "".join(parts)))

    # Find mode length
    lens = Counter(len(s) for _, s in entries)
    mode_len = lens.most_common(1)[0][0]

    # Clip and count
    kept = []
    too_short = 0
    clipped = 0
    at_mode = 0

    for h, s in entries:
        if len(s) < mode_len:
            too_short += 1
        elif len(s) > mode_len:
            clipped += 1
            kept.append((h, s[:mode_len]))
        else:
            at_mode += 1
            kept.append((h, s))

    total = len(entries)
    dropped = total - len(kept)
    drop_pct = 100.0 * dropped / total if total > 0 else 0

    # Write back
    with open(fasta, "w") as out:
        for h, s in kept:
            out.write(f">{h}\n{s}\n")

    print(f"  {gene:12s}  total={total:>10,}  kept={len(kept):>10,}  "
          f"dropped={dropped:>6} ({drop_pct:.3f}%)  "
          f"[too_short={too_short} clipped_3prime={clipped} already_ok={at_mode}]  "
          f"uniform={mode_len}bp")

    log.write(f"{gene}\t{total}\t{len(kept)}\t{dropped}\t{drop_pct:.3f}\t"
              f"{too_short}\t{clipped}\t{at_mode}\t{mode_len}\n")

log.close()
print(f"\nLog written to {log_path}")
