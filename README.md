## Scripts used in synthaser manuscript

| File | Description |
| ---- | ----------- |
| ``sum_bitscores.py`` | Script to sum bitscores of identical query/target sequence pairs in BLAST results. |
| ``extract_PKS.py`` | Script to extract PKS/NRPS sequences from MIBiG GenBank/JSON files. |
| ``mibig/*`` | synthaser output for MIBiG synthases |
| ``ks_domains/*`` | Files generated during KS domain network construction |

## Methods
### Comparison to MIBiG domain architectures
Download MIBiG database GenBank and JSON dumps and extract contents:

	wget https://dl.secondarymetabolites.org/mibig/mibig_json_2.0.tar.gz
	wget https://dl.secondarymetabolites.org/mibig/mibig_gbk_2.0.tar.gz
	tar xzvf mibig_json_2.0.tar.gz
	tar xzvf mibig_gbk_2.0.tar.gz

Retrieve all annotated PKS sequences using ``extract_PKS.py``:

	python3 extract_PKS.py \
		mibig_gbk_2.0/ \               # GenBank folder
		mibig_json_2.0/ \              # JSON folder
		mibig_table.tsv \              # Output, table with MIBiG metadata
		--fasta mibig_synthases.fasta  # Output, FASTA file with PKS 

Setup synthaser:

	pip install --user synthaser

Run synthaser on PKS sequences, saving HTML plot and search session:

	synthaser search \
		--query_file mibig_synthases.fasta \
		--json_file mibig_synthases.json \
		--output mibig_predictions.txt \
		--plot mibig.html \
		--long_form

The MIBiG metadata table (``mibig_table.tsv``) was then merged with
the synthaser predictions table (``mibig_predictions.txt``). MIBiG domain
architectures were copied from the 'NRPS/PKS domains' tab of each
MIBiG entry, added to the table and compared to the predictions in the
synthaser output.

### Creation of the Aspergillus KS network
Retrieve sequences from NCBI containing the ``cond_enzymes`` conserved
domain family, removing any unnecessary information from FASTA description
lines:

	esearch -db cdd -query 238201 |\
		elink -target protein |\
		efilter -query "Aspergillus"[ORGN] -source genbank |\
		efetch -format fasta |\
		sed 's/ .*$//g' - > synthases.faa

Analyse sequences using synthaser:

	synthaser search \
		--query_file synthases.faa \
		--json_file synthases.json \
		--output architectures.txt \
		--long_form

Extract KS domains from the search session:

	synthaser extract \
		synthases.json \  # Session file
		synthases_  \     # Output file prefix, e.g. synthases_KS.faa
		--mode domain \   # Specify domain extraction
		--types KS        # Specify KS domains

Build DIAMOND database from extracted KS domains:

	diamond makedb --in synthases_KS.faa --db KS

Perform all vs all alignments:

	diamond blastp --query domains.faa \
		--db KS.dmnd \
		--more-sensitive \
		--outfmt "6 qseqid sseqid bitscore" \
		--out KS_alignments.tsv

Sum bitscores of all non-overlapping high-scoring segment pairs (HSPs):

	python3 sum_bitscores.py KS_alignments.tsv summed.tsv

The summed alignment table (``summed.tsv``) was then imported into CytoScape
v3.7.2 to build a similarity network. Domain architecture predictions from
synthaser (``architectures.txt``) were imported and connected to their
corresponding nodes, which were then coloured based on an alphabetical ordering
of architectures.
