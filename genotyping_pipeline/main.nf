nextflow.enable.dsl = 2

data = Channel
        .fromPath(params.sample_tsv, type: "file", checkIfExists: true)
        .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
        .map { row -> tuple(row.sample_id, row.alignment_file) }
params.reference_fasta = file(params.reference_fasta, type: "file", checkIfExists: true)
params.reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)
params.sv_bed = file(params.sv_bed, type: "file", checkIfExists: true)
params.snv_vcf = file(params.snv_vcf, type: "file", checkIfExists: true)
params.pedigree_file = file(params.pedigree_file, type: "file", checkIfExists: true)
params.clf_folder = file(params.clf_folder, type: "file", checkIfExists: true)
params.exclude_regions_bed = file(params.exclude_regions_bed, type: "file", checkIfExists: true)

// Run the SV2 genotyper:
process SV2 {
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("sv2_output/$sample_id*.vcf.gz")
    script:
    """
    python run_sv2.py --alignment_file $alignment_file --reference_fasta $reference_fasta --snv_vcf_file $snv_vcf --sv_bed_file $sv_bed --sample_name $sample_id --ped_file $pedigree_file --exclude_regions_bed $exclude_regions_bed --output_folder sv2_output --clf_folder $clf_folder
    """
    stub:
    """
    touch $sample_id*.vcf.gz
    """
}

workflow {
    SV2(data)
}
