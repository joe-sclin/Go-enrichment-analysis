import os
from argparse import ArgumentParser
from GO_enrichment import enrichment_analysis


def get_argument_parser():
    p = ArgumentParser(description="Perform GO enrichment analysis with optional background genes.")
    p.add_argument("gene", help="gene list txt file")
    p.add_argument("--bg", default="", help="background txt file")
    p.add_argument("--namespace", default="BP", choices=("BP", "CC", "MF"),
        help="namespace applied for analysis")
    p.add_argument("--source", default='NCBI', choices=("NCBI", "GO"),
        help="source of gene annotation file")
    return p

def main(args):
    if os.path.exists(args.gene):
        with open(args.gene, 'r') as file:
            gene_list = [line.strip() for line in file]
    else:
        raise FileExistsError('Input gene list file not exists.')
    if args.bg != "":
        if os.path.exists(args.bg):
            with open(args.bg, 'r') as file:
                background = [line.strip() for line in file]
        else:
            raise FileExistsError('Input background file not exists.')
    else:
        background = []

    enrichment_analysis(gene_list, background, name=args.namespace, source=args.source)
if __name__ == "__main__":
    p = get_argument_parser()
    args = p.parse_args()
    main(args)

