#!/bin/bash -l
set -eu

src=$(realpath ../src/kegg-plots.R)
#src=$(realpath ../src/biogrid-plots.R)
out=$(realpath ../analysis)

sbatch -A facility -w franklin -n 1 -p core --mem=400G --mail-type=ALL \
--mail-user=nicolas.delhomme@umu.se \
-o seidr_plots_kegg.out -e seidr_plots_kegg.err \
$(realpath ../UPSCb-common/pipeline/runKnitR.sh) $src
