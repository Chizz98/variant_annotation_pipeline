"""

"""
import os
import subprocess
import gzip
import re


def parse_gene_tsv(gene_vcf_fn: str) -> list:
    """ Reads gene name, position, and start and end positions from a csv

    :param gene_vcf_fn: Filename of the input csv
    :return: List of gene names
    """
    out_list = []

    with open(gene_vcf_fn, "r") as infile:
        header = infile.readline().strip("\n").split("\t")
        for line in infile:
            line = line.strip("\n").split("\t")
            gene = line[2]
            if gene.startswith("LOC"):
                out_list.append(gene)
    return out_list


def query_feature_table(
        feature_table_fn: str,
        genes: list
) -> tuple[dict, dict]:
    """ Queries an NCBI feature table to extract gene information

    :param feature_table_fn: Filename of the NCBI feature table
    :param genes: List of genes to extract from the feature table
    :return: A dictionary with gene names input and the chr, start and stop as
        values. A dictionary with protein ids as the key and gene, mrna, chr,
        start and stop as values
    """
    gene_dict = {}
    protein_dict = {}

    with open(feature_table_fn, "r") as infile:
        for line in infile:
            line = line.strip().split("\t")
            gene_name = line[14]
            feature = line[0]
            if gene_name in genes:
                gene = list(set(line) & set(genes))[0]
                chr, start_pos, stop_pos = line[6:9]
                if feature == "gene":
                    gene_dict[gene] = [chr, start_pos, stop_pos]
                if feature == "mRNA":
                    protein_id = line[12]
                    mrna_id = line[10]
                    protein_dict[protein_id] = [
                        gene, mrna_id, chr, start_pos, stop_pos
                    ]
    return gene_dict, protein_dict


def extract_gene(
        chromosome: str,
        start_pos: int,
        end_pos: int,
        vcf_fn: str,
        out_fn: str,
        run_dry: bool = False
):
    """ Runs a bcftools command to extract given position from input vcf

    :param chromosome: Contig name
    :param start_pos: Start of region to extract
    :param end_pos: End of region to extract
    :param vcf_fn: Filename of the input vcf
    :param out_fn: Filename of the out vcf
    :param run_dry: If True, prints command to std out instead of excecuting it
    :return: None, runs bcftools to write an output file
    """
    command = f"bcftools view -r {chromosome}:{start_pos}-{end_pos} {vcf_fn} " \
              f"-Oz -o {out_fn}"

    if run_dry:
        print(command)
    else:
        subprocess.run(command, shell=True)


def parse_snpeff_to_provean(snpeff_vcf_fn: str) -> dict:
    """ Takes a snpeff annotated vcf and writes provean compatible variant file

    :param snpeff_vcf_fn: snpeff annotated vcf filename
    :return: dictionary with transcript ids as keys and provean formatted
    variants as values
    """
    open_vcf = gzip.open if snpeff_vcf_fn.endswith(".gz") else open

    out_dict = {}

    with open_vcf(snpeff_vcf_fn, "r") as vcf:
        for line in vcf:
            # Skip vcf header
            if line.startswith("#"):
                continue

            line = line.split()
            info_field = line[7]
            re_match = re.match(r".*(ANN[^;]+)", info_field)
            annotations = re_match.groups()[0].split(",")

            for annotation in annotations:
                annotation = annotation.split("|")
                # Skip non moderate annotations
                if not annotation[2] == "MODERATE":
                    continue
                rna_id = annotation[6].replace("rna-", "")
                aa_change = annotation[10]

                variant = var_snpeff_to_provean(aa_change)

                if rna_id in out_dict:
                    out_dict[rna_id] += [variant]
                else:
                    out_dict[rna_id] = [variant]
    return out_dict


def var_snpeff_to_provean(snpeff_var: str) -> str:
    """ Changes snpeff variant format to provean variant format

    :param snpeff_var: variant in snpeff format
    :return: variant in provean format
    """
    aa_lut = {
        'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
        'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
        'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
        'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
    }

    variant = None

    # Single substitution
    m_sub = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$",
                     snpeff_var)
    if m_sub:
        ref = aa_lut[m_sub.group(1)]
        pos = m_sub.group(2)
        alt = aa_lut[m_sub.group(3)]
        variant = f"{ref}{pos}{alt}"

    # Deletion
    m_del = re.match(
        r"p\.([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?del$",
        snpeff_var)
    if m_del:
        start = m_del.group(2)
        end = m_del.group(4) if m_del.group(4) else start
        variant = f"DEL{start}" if start == end \
            else f"DEL{start}_{end}"

    # Insertion
    m_ins = re.match(
        r"p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins"
        r"([A-Za-z]+)$",
        snpeff_var)
    if m_ins:
        start = m_ins.group(2)
        end = m_ins.group(4)
        inserted_aa_3 = re.findall(r"[A-Z][a-z]{2}", m_ins.group(5))
        inserted_aa_1 = ''.join(
            [aa_lut[aa] for aa in inserted_aa_3])
        variant = f"INS{start}_{end}:{inserted_aa_1}"
    return variant


def extract_protein_sequences(
        protein_fa_fn: str,
        protein_ids: list,
        out_fn: str
) -> None:
    """ Extracts protein sequences from a protein fasta

    :param protein_fa_fn: input fasta with protein sequences.
    :param protein_ids: list of proteins to extract from the input fasta.
    :param out_fn: filename to write the output fasta to.
    :return: None, writes outfile
    """
    in_prot_of_interest = False
    out_lines = []
    out_dict = {}
    with open(protein_fa_fn, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                fasta_header = line.strip().split(" ")
                protein_id = fasta_header[0].replace(">", "")
                if protein_id in protein_ids:
                    in_prot_of_interest = True
                    out_dict[protein_id] = line.strip()
                    line = f">{protein_id}\n"
                else:
                    in_prot_of_interest = False
            if in_prot_of_interest:
                out_lines.append(line)
    with open(out_fn, "w") as outfile:
        for line in out_lines:
            outfile.write(line)
    return out_dict


def write_regions_file(gene_dict: dict, out_fn: str) -> None:
    """ Writes a regions file to be parsed by bcftools

    :param gene_dict: first dictionary output by query_feature_table
    :param out_fn: name of the output file
    :return: None, writes an output file
    """
    with open(out_fn, "w") as outfile:
        for value in gene_dict.values():
            outfile.write("\t".join(value) + "\n")


def main():
    # Get cmd arguments

    gene_dir = "genes"
    protein_dir = "proteins"
    if not os.path.exists(gene_dir):
        os.mkdir(gene_dir)
    if not os.path.exists(protein_dir):
        os.mkdir(protein_dir)

    # Parse input file
    gene_list = parse_gene_tsv("genes_test.txt")
    
    # Parse feature table
    gene_dict, protein_dict = query_feature_table(
        "../GCF_002870075.3_Lsat_Salinas_v8_feature_table.txt",
        gene_list
    )
    
    # Extract protein sequences
    protein_fa = "proteins/protein.fa"
    
    extract_protein_sequences(
        "../GCF_002870075.3_Lsat_Salinas_v8_protein.faa",
        list(protein_dict.keys()),
        protein_fa
    )
    
    # Extract regions of interest from annotated vcf
    regions_file = "genes/regions.txt"
    input_vcf = "../Lser_200_filtered_missing_lt100_het_lt20.ann.vcf.gz"
    out_vcf = "genes/test.ann.vcf"

    write_regions_file(gene_dict, regions_file)
    bcftools_cmd = f"bcftools view -R {regions_file} {input_vcf} -o {out_vcf}"
    subprocess.run(bcftools_cmd, shell=True)
    
    # Write provean variant file
    prov_variant_file = "proteins/provean_vars.txt"
    
    out_dict = parse_snpeff_to_provean(out_vcf)
    with open(prov_variant_file, "w") as outfile:
        for mrna_id, prov_variants in out_dict.items():
            protein_id_match = None
            for protein_id, prot_vals in protein_dict.items():
                if mrna_id in prot_vals:
                    protein_id_match = protein_id
            if protein_id_match is not None:
                for variant in prov_variants:
                    if variant is not None:
                        outfile.write(f"{protein_id_match}\t{variant}\n")
    
    # Run interpro query
    interpro_cmd = f"interproscan.sh -i {protein_fa} -f tsv -dp -appl Pfam"
    subprocess.run(interpro_cmd, shell=True)

    """
    # Extract genes
    extract_gene("CHR", 12, 13, "infile", "outfile", True)

    # Write provean variant file
    parse_snpeff_to_provean("test_gene.ann.vcf")
    """


if __name__ == "__main__":
    main()
