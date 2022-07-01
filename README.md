# seidr-manuscript
Public code related to the seidr manuscript


## Index

* `src/srasalmon.py` - Python script used to preprocess one sample from the SRA using SortMeRNA -> Trimmomatic -> Salmon
* `src/{at,dm,sc}-biogrid-plots.R` - R plotting script to create BioGRID plots in supplementary S2.
* `src/yeastract-plots.R` - R plotting script to create Yeastract plots in supplementary S2.
* `src/kegg-plots.R` - R plotting script to create KEGG plots in supplementary S2.
* `src/enr-script.R` - R script used to run GO and MapMan enrichment.
* `src/gsea-drought.R` - R script used for drought stress GSEA analysis and visualization.
* `src/{athal,dmel,scere}_gs.R` - R scripts used to create KEGG "gold" standards for benchmarking.
* `src/normalize_exp.R` - R script used to import salmon counts and variance stabilize / median center.
