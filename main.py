"""

"""
import os
import subprocess
import gzip
import re
import argparse as arg
import pathlib


def arg_reader():
    """ Reads arguments from command line

    :return class containing the arguments
    """
    arg_parser = arg.ArgumentParser(
        description="Takes a list of genes and isolates the regions spanning "
                    "the mRNA from the input VCF. Also isolates protein "
                    "sequences and optionally runs those through interproscan."
    )
    arg_parser.add_argument(
        "vcf",
        help="The input vcf file."
    )

    database_group_top = arg_parser.add_argument_group(
        title="Database method",
        description="Either supply a database which is in the database.config "
                    "in the main project folder with -d, or supply a feature "
                    "table and protein fasta yourself using -f and -p."
    )
    database_group = database_group_top.add_mutually_exclusive_group(
        required=True
    )
    database_group.add_argument(
        "-f", "--feature_table",
        help="Filename of feature table for your organism, downloadable from "
             "NCBI. Mutually exclusive with -d. Requires protein fasta to be "
             "supplied through -p."
    )
    database_group.add_argument(
        "-p", "--protein_fasta",
        help="Fasta file containing the protein sequences of your organism."
             "Mutually exclusive with -d. Requires feature table to be supplied"
             " trough -f."
    )
    database_group.add_argument(
        "-d",
        "--database",
        help="Internal database to query, mutually exclusive with -f and -p."
    )

    gene_input_group_top = arg_parser.add_argument_group(
        title="Gene input",
        description="Mode to provide genes."
    )
    gene_input_group = gene_input_group_top.add_mutually_exclusive_group(
        required=True
    )
    gene_input_group.add_argument(
        "-g",
        "-genes",
        help="Comma separated list of genes (or a single gene). No spaces "
             "between gene names."
    )
    gene_input_group.add_argument(
        "-G",
        "-gene_file",
        help="Name of a file containing the gene names, with each line "
             "containing one gene name."
    )

    arg_parser.add_argument(
        "-o",
        "--out_dir",
        default="./out",
        help="The output directory to write files to, ./out by default."
    )
    arg_parser.add_argument(
        "-i",
        "--interproscan",
        action="store_true",
        help="If set, runs interproscan on all proteins (just pfam for now)."
    )
    return arg_parser


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


def parse_db_config(config_fn: str) -> dict:
    """ Parser for database.config file

    :param config_fn: name of the config file
    :return: dictionary containing the database name as keys and the related
        urls as values.
    """
    out_dict = {}
    with open(config_fn) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            db_name = line[0]
            protein_url = line[1]
            feature_url = line[2]
            out_dict[db_name] = {
                "protein_url": protein_url,
                "feature_url": feature_url
            }
    return out_dict


def database_handler(script_dir: str, curr_db_dir: str, database: str):
    """ Handles the logic for checking for and retrieving databases

    :param database: Database name as in config / databases folder in project
    :param script_dir: The dir the project is in
    :param curr_db_dir: The dir of the provided database
    :return: filenames of the protein fasta and feature table in the provided
        database.
    """
    all_protein_fa = os.path.join(curr_db_dir, 'proteins.fa')
    feature_table = os.path.join(curr_db_dir, 'features.txt')
    if not (os.path.exists(all_protein_fa) & os.path.exists(feature_table)):
        config_dbs = parse_db_config(
            os.path.join(script_dir, "database.config")
        )
        if database in config_dbs:
            if not os.path.exists(curr_db_dir):
                os.mkdir(curr_db_dir)
            protein_url = config_dbs[database]["protein_url"]
            feature_url = config_dbs[database]["feature_url"]

            # Download data
            prot_cmd = f"wget -O {all_protein_fa + '.gz'} {protein_url}"
            feat_cmd = f"wget -O {feature_table + '.gz'} {feature_url}"
            subprocess.run(prot_cmd, shell=True)
            subprocess.run(feat_cmd, shell=True)

            # Unzip data
            subprocess.run(f"gunzip {all_protein_fa + '.gz'}", shell=True)
            subprocess.run(f"gunzip {feature_table + '.gz'}", shell=True)

        else:
            raise Exception(f"Database {database} not found in database.config")
    return all_protein_fa, feature_table



def main():
    script_dir = pathlib.Path(__file__).resolve().parent
    db_dir = os.path.join(script_dir, "databases")

    # Get cmd arguments
    arg_parser = arg_reader()
    args = arg_parser.parse_args()

    input_vcf = args.vcf
    database = args.database
    run_interpro = args.interproscan
    out_dir = args.out_dir
    genes = args.genes
    gene_list = args.gene_file

    if bool(args.protein_fasta) != bool(args.feature_table):
        arg_parser.error("-f and -p have to be submitted together")
    if args.protein_fasta is not None and args.feature_table is not None:
        all_protein_fa = args.protein_fasta
        feature_table = args.feature_table
        if database:
            arg_parser.error("-f and -p options are incompatible with -d")
    else:
        if not os.path.exists(db_dir):
            os.mkdir(db_dir)
        curr_db_dir = os.path.join(db_dir, database)
        all_protein_fa, feature_table = database_handler(
            script_dir, curr_db_dir, database
        )

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    regions_file = os.path.join(out_dir, "regions.txt")
    out_vcf = os.path.join(out_dir, "genes.vcf")

    # Parse input file
    if genes:
        gene_list = genes.split(",")
    else:
        # To do rework for described file format in argparser
        gene_list = parse_gene_tsv("genes_test.txt")
    
    # Parse feature table
    gene_dict, protein_dict = query_feature_table(
        feature_table,
        gene_list
    )

    protein_fa = os.path.join(out_dir, "proteins.fa")
    extract_protein_sequences(
        all_protein_fa,
        list(protein_dict.keys()),
        protein_fa
    )
    
    # Extract regions of interest from annotated vcf
    write_regions_file(gene_dict, regions_file)
    bcftools_cmd = f"bcftools view -R {regions_file} {input_vcf} -o {out_vcf}"
    subprocess.run(bcftools_cmd, shell=True)
    
    # Run interpro query
    if run_interpro:
        interpro_cmd = f"interproscan.sh -i {protein_fa} -f tsv -dp " \
                       f"-appl Pfam -o " \
                       f"{os.path.join(out_dir, 'interpro_results.tsv')}"
        subprocess.run(interpro_cmd, shell=True)


if __name__ == "__main__":
    main()
