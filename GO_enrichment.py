import pandas as pd
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.gaf_reader import GafReader
from goatools.base import download_ncbi_associations, download_go_basic_obo, dnld_gaf
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
import networkx as nx
from networkx.drawing.nx_agraph import to_agraph
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from time import time


def enrichment_analysis(gene_list: list, background: list, name: str = 'BP', source: str = 'NCBI'):
    """
    Perform GO terms enrichment analysis using input gene list and optional backgroun gene
    Output: csv containing information about significantly enriched / purified GO terms and network diagram dot for each of
    analysis result if significantly GO terms were found (Enriched: red and Purified: green)
    """
    start = time()

    download_obo = download_go_basic_obo()  # GOATOOLs function to download Ontologies (obo)
    obo = GODag(download_obo)
    if source == 'NCBI':
        # Section: preparation of necessary files
        human_gene = pd.read_csv('human_gene.txt', sep='\t', index_col=None, usecols=['GeneID', 'Symbol'])
        genome = dict(zip(human_gene['Symbol'], human_gene['GeneID']))
        # GOATOOLS function to download associations ('gene2go')
        _ = download_ncbi_associations()
        ns2assoc = Gene2GoReader('gene2go', taxids=[9606]).get_ns2assc()  # TaxID of human: 9606
    else:
        dnld_gaf(species_txt='goa_human')
        gaf = GafReader('goa_human.gaf')
        ns2assoc = gaf.get_ns2assc(taxid=[9096])
        genome = {}
        # If statement to skip duplicated entries (Reduce runtime)
        genome = {nt.DB_Symbol: nt.DB_ID for nt in gaf.associations if nt.DB_Symbol not in genome}

    if background != []:  # mapping of background gene list to source data
        background_mapper = {k: v for k, v in genome.items() if
                             k in background}  # For third analysis (Disease genes as population)
        background_inv_map = {v: k for k, v in background_mapper.items()}
    else:
        background_mapper = genome
        background_inv_map = {v: k for k, v in genome.items()}

    goeaobj = GOEnrichmentStudy(  # Create goea object for analysis considering selected namespaces
        background_inv_map.keys(),  # Population
        ns2assoc[name],  # geneid / GO associations
        obo,  # Ontologies
        propagate_counts=False,
        methods=['fdr_bh'])  # defult multipletest correction method
    GO_items = set()
    for item in ns2assoc[name].values():
        GO_items.update(item)
    GO_items = list(GO_items)

    def go_analysis(test_genes):
        print('Input genes:', len(test_genes))
        mapped_genes = [background_mapper.get(gene) for gene in test_genes if gene in background_mapper]
        print('mapped genes:', len(mapped_genes))
        print('GO enrichment in namespace', name, 'using', source, 'data')
        goea_results_all = goeaobj.run_study(mapped_genes)
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

        rows = [(x.GO, x.goterm.name, x.goterm.namespace, x.enrichment, x.p_uncorrected, x.p_fdr_bh,
                 x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),
                 [background_inv_map[y] for y in x.study_items]) for x in goea_results_sig]
        columns = ['GO', 'term', 'class', 'e_or_p', 'p', 'p_corr', 'n_genes', 'n_study', 'n_go', 'study_genes']
        GO = pd.DataFrame(rows, columns=columns)  # Res df for return
        GO = GO[GO.n_genes > 1]  # Exclude GO terms with only 1 gene (Not representative)
        enriched = []
        purified = []
        for x in goea_results_sig:  # Generate enriched and purified list for GO network drawing
            if x.enrichment == 'e':
                enriched.append(x.GO)
            else:
                purified.append(x.GO)
        return GO, enriched, purified

    def drawing_GO_network(purified=None, enriched=None):
        if purified is None: purified = []
        if enriched is None: enriched = []
        root = ['GO:0008150', 'GO:0005575', 'GO:0003674']  # GO terms for 3 namespaces
        go_terms = purified + enriched
        GO_network = nx.DiGraph()

        # Recursive function to compute network until reaching namespaces
        def adding_term_to_network(current_term, network):
            if current_term.id not in network.nodes():
                network.add_node(current_term.id, label=current_term.name)
            while current_term.level > 0:  # level = 0 -> namespace
                parents = current_term.parents
                for parent in parents:
                    if parent.id not in network.nodes():
                        network.add_node(parent.id, label=parent.name)
                        network.add_edge(current_term.id, parent.id)
                        adding_term_to_network(parent, network)
                    else:
                        network.add_edge(current_term.id, parent.id)
                current_term = parent

        for term in go_terms:
            adding_term_to_network(obo[term], GO_network)
        A = to_agraph(GO_network)  # Agraph object for dot format drawing
        for n in GO_network.nodes:
            node = A.get_node(str(n))
            node_label = f"{n}\n{node.attr['label']}"
            node.attr['label'] = node_label
            if n in root:  # Namespace: lightblue
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = 'lightblue'
            if n in purified:  # Purified: green
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = '#99FF99'
            if n in enriched:  # Enrichment: red
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = '#FF7C80'
        A.graph_attr['rankdir'] = 'BT'  # Namespace on top
        with open(''.join(['GO_network_', name,'_', source, '.dot']), 'w') as file:
            file.write(A.to_string())

    df, e_list, p_list = go_analysis(gene_list)
    if df.shape[0] > 0:  # Only genrerate csv and plotting bar chart when significant terms found
        df.to_csv('_'.join(['GO_terms_', name, source, '.csv']),
                  index=False)
        class_colors = {'Overrepresented': 'red', 'Underrepresented': 'green'}
        label_map = {'e': 'Overrepresented', 'p': 'Underrepresented'}
        namespace_map = {'BP': 'Biological processes', 'CC': 'Cellular components', 'MF': 'Metabolic functions'}
        class_palette = [class_colors[label_map[label]] for label in df['e_or_p'].unique()]
        fig, ax = plt.subplots(figsize=(12, max(df.shape[0], 4)))
        plt.subplots_adjust(left=0.2)
        sns.barplot(data=df, x='n_genes', y='term', width=0.2, ax=ax, palette=class_palette, hue='e_or_p',
                    dodge=False, alpha=0.7)
        ax.set_yticklabels([textwrap.fill(e, 22) for e in df['term']])  # Max 22 char per line (for long GO Term)
        plt.subplots_adjust(wspace=0.05)
        ax.set_title(' '.join(['GO enrichment chart for namespace', namespace_map[name]]))  # add title to the plot
        handles, labels = ax.get_legend_handles_labels()
        labels = [label_map[label] for label in labels]
        ax.legend(handles, labels, loc=0, borderaxespad=0)
        ax.set_xlabel('Number of genes')
        ax.set_ylabel('GO terms')
        fig.savefig('_'.join(['GO_chart', name, source, '.png']))
        plt.close()
        drawing_GO_network(p_list, e_list)
    print('GO terms enrichment analysis completed in', round(time() - start, 3), 's.')
