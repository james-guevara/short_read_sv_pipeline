import pysam
import sys

# Read the input VCF file:
delly = sys.argv[1]
out = sys.argv[2]
mybcf = pysam.VariantFile(delly, "r")

# New VCF file to be written:
outvcf = pysam.VariantFile(out, "w", header = mybcf.header)
for record in mybcf:
    svlen = record.stop - record.start
    info = record.info
    if info["SVTYPE"] != "INS":
        record.info.__setitem__("SVLEN", svlen)
    outvcf.write(record)

