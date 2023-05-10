rm -r logs/*
snakemake -j 1 --configfile config/config.yml --unlock
snakemake --profile profile --snakefile mainSnake --use-conda
    
