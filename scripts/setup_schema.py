import os
import json
import tiledb
import numpy as np

db_uri = snakemake.params.db_uri
chroms = snakemake.params.chroms
samples = snakemake.params.samples
max_coord = snakemake.params.max_coord
output_map = snakemake.output.sample_map
output_flag = snakemake.output.flag

os.makedirs(os.path.dirname(output_flag), exist_ok=True)

sample_map = {name: i + 1 for i, name in enumerate(samples)}
with open(output_map, 'w') as f:
    json.dump(sample_map, f, indent=4)

for chrom in chroms:
    chrom_uri = os.path.join(db_uri, chrom)
    if not os.path.exists(chrom_uri):
        zstd_level = 19
        zstd_filter = tiledb.FilterList([tiledb.ZstdFilter(level=zstd_level)])
        coord_filter = tiledb.FilterList([tiledb.DoubleDeltaFilter(), tiledb.ZstdFilter(level=zstd_level)])
        
        dim_start = tiledb.Dim(name="start", domain=(1, max_coord), tile=1000000, dtype=np.int32, filters=coord_filter)
        dim_end = tiledb.Dim(name="end", domain=(1, max_coord), tile=1000000, dtype=np.int32, filters=coord_filter)
        dim_sample = tiledb.Dim(name="sample", domain=(1, 65535), dtype=np.uint16, filters=zstd_filter)
        
        dom = tiledb.Domain(dim_start, dim_end, dim_sample)
        score_filter = tiledb.FilterList([tiledb.DeltaFilter(), tiledb.ZstdFilter(level=zstd_level)])
        attr_score = tiledb.Attr(name="score", dtype=np.int32, filters=score_filter)
        
        schema = tiledb.ArraySchema(
            domain=dom, sparse=True, attrs=[attr_score],
            cell_order='row-major', tile_order='row-major', capacity=500000
        )
        tiledb.Array.create(chrom_uri, schema)
        
        with tiledb.SparseArray(chrom_uri, mode='w') as A:
            A.meta["sample_map"] = json.dumps(sample_map)

with open(output_flag, 'w') as f:
    f.write("done")
