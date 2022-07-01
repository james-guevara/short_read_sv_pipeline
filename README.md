# Illumina short-read structural variation pipeline (using Nextflow)

Example run command:
```
nextflow run main.nf --reference_fasta /expanse/lustre/projects/ddp195/eiovino/fasta/Homo_sapiens_assembly38.fasta \
		     --sample_alignments_tsv input_REACH2.csv --outdir /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/SV_result_REACH/REACH_cohort \
		     --delly_exclude_regions_bed data/exlude.regions.delly.human.hg38.excl.tsv \
		     --smoove_exclude_regions_bed data/exclude.smoove.cnvnator_100bp.GRCh38.20170403.bed \
		     --expansion_hunter_variant_catalog_json Hg38/ExpansioHunter_variant_catalog.json \
		     --recode_delly_python_script data/recode_delly.py \
		     --mosdepth_segmental_duplications_bed data/Segmental_dups_hg38_frt_srt.bed \
		     --bind_path /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf,/expanse/lustre/projects/ddp195/eiovino/fasta,/expanse/lustre/projects/ddp195/eiovino/test_pipeline/Sv_pipeline/data,/expanse/lustre/projects/ddp195/eiovino/test_pipeline/Sv_pipeline/Hg38,/expanse/projects/sebat1/genomicsdataanalysis/REACH_MSSNG/request_202001/ \
		     -profile slurm \
		     -resume 
```
