import os
import tiledb

db_uri = snakemake.params.db_uri
chrom = snakemake.wildcards.chrom
output_flag = snakemake.output.flag
threads = snakemake.threads

chrom_uri = os.path.join(db_uri, chrom)

consolidation_config = tiledb.Config({
    "sm.mem.total_budget": "64424509440",
    "sm.consolidation.buffer_size": "500000000",
    "sm.consolidation.step_min_frags": "2",
    "sm.consolidation.step_max_frags": "20",
    "sm.compute_concurrency_level": str(threads),
    "sm.io_concurrency_level": str(threads),
    "vfs.num_threads": str(threads)
})
main_ctx = tiledb.Ctx(consolidation_config)

tiledb.consolidate(chrom_uri, config=consolidation_config, ctx=main_ctx)
tiledb.vacuum(chrom_uri, config=consolidation_config, ctx=main_ctx)

with open(output_flag, 'w') as f:
    f.write("done")
