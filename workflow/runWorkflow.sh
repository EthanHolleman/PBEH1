rm -r logs/*
snakemake -j 1 --configfile config/config.yml --unlock
snakemake \
    --snakefile Snakefile \
    --profile profile/ \
    --configfile config/config.yml \
    --use-conda \
    --rerun-incomplete \
    -k