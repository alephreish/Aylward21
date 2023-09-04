cfg=
if [ -s config.yaml ]; then
    cfg="--configfile config.yaml"
fi
snakemake -c60 --use-conda --resources ncbi=1 $cfg "$@"
