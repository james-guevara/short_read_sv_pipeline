nextflow.enable.dsl = 2

data = Channel
        .fromPath("mycohort.csv", type: "file", checkIfExists: true)
        .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
        .map { row -> tuple(row.sample_id, row.alignment_file) }
params.reference_fasta = file("/expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", type: "file", checkIfExists: true)
params.reference_fasta_fai = file("${params.reference_fasta}.fai", checkIfExists: true)
params.sv_bed = file("/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/results/HG002/raw/chr22.bed", type: "file", checkIfExists: true)
params.snv_vcf = file("/expanse/lustre/projects/ddp195/j3guevar/sv_testing/resources/snv_indel_vcfs_for_cerds/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", type: "file", checkIfExists: true)
params.pedigree_file = file("/home/j3guevar/sv2_2/sv2_old_version/sv2/f01.ped", type: "file", checkIfExists: true)
params.clf_folder = file("/home/j3guevar/sv2_2/sv2_old_version/sv2/data/trained_classifiers", type: "file", checkIfExists: true)
params.exclude_regions_bed = file("/home/j3guevar/sv2_2/data/excluded_regions_bed_files/hg38_excluded.bed.gz", type: "file", checkIfExists: true)

// Run the SV2 genotyper:
process SV2 {
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("sv2_output/${sample_id}_*.vcf.gz")
    script:
    """
    python run_sv2.py --alignment_file $alignment_file --reference_fasta $reference_fasta --snv_vcf_file $snv_vcf --sv_bed_file $sv_bed --sample_name $sample_id --ped_file $pedigree_file --exclude_regions_bed $exclude_regions_bed --output_folder sv2_output --clf_folder $clf_folder -M
    """
    stub:
    """
    mkdir sv2_output
    touch sv2_output/${sample_id}.vcf.gz
    """
}

workflow {
    SV2(data)
}
