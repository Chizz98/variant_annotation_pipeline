"""

"""
import os
import subprocess
import gzip
import re


def parse_gene_tsv(gene_vcf_fn: str) -> list:
    """ Reads gene name, position, and start and end positions from a csv

    :param gene_vcf_fn: Filename of the input csv
    :return: list of gene names
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


def query_feature_table(feature_table_fn: str, genes: list):
    """

    :param feature_table_fn:
    :param genes:
    :return:
    """
    out_dict = {}

    with open(feature_table_fn, "r") as infile:
        for line in infile:
            line = line.strip().split("\t")
            print(line)


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
            else f"DEL{start}-{end}"

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
        variant = f"INS{start}-{end}:{inserted_aa_1}"

    return variant


def main():
    # Get cmd arguments

    gene_dir = "genes"
    protein_dir = "proteins"
    if not os.path.exists(gene_dir):
        os.mkdir(gene_dir)
    if not os.path.exists(protein_dir):
        os.mkdir(protein_dir)

    # Parse input file
    parse_gene_tsv("genes.txt")

    """
    # Extract genes
    extract_gene("CHR", 12, 13, "infile", "outfile", True)

    # Write provean variant file
    parse_snpeff_to_provean("test_gene.ann.vcf")
    """


if __name__ == "__main__":
    main()
