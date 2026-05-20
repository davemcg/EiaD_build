import os
import json
import pyBigWig
import tiledb
import numpy as np

file_path = snakemake.params.file_path
db_uri = snakemake.params.db_uri
sample_map_file = snakemake.input.sample_map_file
output_flag = snakemake.output.flag
sample_name = snakemake.wildcards.sample
threads = snakemake.threads

with open(sample_map_file, 'r') as f:
    sample_map = json.load(f)
sample_id = sample_map[sample_name]

# Bind TileDB explicitly to the threads requested by Snakemake
config = tiledb.Config({
    "sm.compute_concurrency_level": str(threads),
    "sm.io_concurrency_level": str(threads),
    "vfs.num_threads": str(threads)
})
ctx = tiledb.Ctx(config)

bw = pyBigWig.open(file_path)
for chrom in bw.chroms().keys():
    intervals = bw.intervals(chrom)
    if not intervals:
        continue
        
    arr = np.array(intervals)
    starts = arr[:, 0].astype(np.int32) + 1
    ends = arr[:, 1].astype(np.int32)
    scores = np.round(arr[:, 2].astype(np.float64) * 1000).astype(np.int32)
    samples = np.full(len(starts), sample_id, dtype=np.uint16)
    
    chrom_uri = os.path.join(db_uri, chrom)
    with tiledb.SparseArray(chrom_uri, mode='w', ctx=ctx) as A:
        A[starts, ends, samples] = {"score": scores}
        
bw.close()

with open(output_flag, 'w') as f:
    f.write("done")
