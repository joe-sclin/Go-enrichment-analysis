# Go enrichment analysis
A convenient Python script to perform GO enrichment analysis using multitest correlation method: Benjamini/Hochberg FDR

### Features of script
* Automatic download Gene Ontologies file (Default to use go-basic.obo)
* Two different source of gene annotation files are available
* Gene Annotation file
  * source option: GO (Gene Ontology Resource): Automatic download (Default to use goa_human.gaf)
* source option: NCBI: provided
```
Download human_gene.txt from NCBI with following search details to include all possible cases for miRNAs
https://www.ncbi.nlm.nih.gov/gene
("Homo sapiens"[Organism] AND alive[property]) AND ((mirna[All Fields] AND (("genetype miscrna"[Properties] OR
"genetype ncrna"[Properties] OR "genetype scrna"[Properties] OR "genetype snrna"[Properties] OR
"genetype snorna"[Properties]) NOT "genetype protein coding"[Properties] AND "annotated genes"[Filter]))
OR genetype protein coding[Properties] OR summary:microRNA)
```
* Optional: input background gene list
* Users can select which namespace applied for analysis to reduce runtime (analysis result for a specific namespace only was common for many researcher)

* Number of input and mapped genes --> evaluate quality of input gene list
* 3 output files for different needs of users
  * CSV file contained information about each over- or under- represented GO terms
  * Bar chart to show number of genes included in each over- or under- represented GO terms
  * DOT file to draw the network about relationship of each over- or under- represented GO terms (providing DOT file is more flexible than directly providing the graph)

### Requirement

* The list of package version requirements is available here:
  * `argparse` is used for program argument and command line option management
  * `goatools` is the main tools for GO enrichment analysis
  * `pandas` is used for csv formatting and exporting
  * `networkx` is used for network graph DOT
  * `matplotlib` and `seaborn` is used for bar chart
```
argparse 1.4.0
pandas 1.5.3
goatools 1.3.1
networkx 3.0
matplotlib 3.6.1
seaborn 0.12.2
```

### User guide:
* Basic usage: input.txt is a file each line containing the gene in your input gene list

  `python main.py input.txt` 
* Applying background gene list with option `--bg`:
  
  `python main.py input.txt --bg background.txt` 
* Changing the namespace (Default BP) with option `--namespace`:
  
  `python main.py input.txt --bg background.txt --namespace "CC"`
* Changing the source (Default NCBI) with option `--source`:
  
  `python main.py input.txt --bg background.txt --namespace "CC" --source "GO"`

### Format of input:
* Strongly recommend to use official gene symbol that approved by HGNC to improve the mapping rate

### Example data input:
* `BRCA.txt` is a list of known cancer-associated genes for human neoplasms: Breast Carcinoma
* Run example data: `python main.py example_data\BRCA.txt`
* Source of data: DisGeNET using UMLS: `C0678222` with selection filter `DSI >= 0.5 and DPI <= 0.5`

### Development:
* Adding the multiple symbol or ID conversion to HGNC official symbol as a separated script