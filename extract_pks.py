"""
Extract PKS/NRPS sequences from MIBiG JSON/GBK files.

Generates a table with each sequence and their corresponding MIBiG entries.
"""

import re
import json
import argparse

from pathlib import Path


def parse_gbk(fp):
    """Parse a MIBIG GenBank file.

    Gets the Protein ID, amino acid sequence and if present, the gene name.
    """
    pattern = re.compile(
        r'     CDS    (.+?)/protein_id="(.+?)".+?/translation="(.+?)"',
        re.DOTALL,
    )
    raw = fp.read()
    results = []
    for gene in pattern.finditer(raw):
        group = gene.group(1)
        conditions = [
            'NRPS_PKS="type: NRPS' in group,
            'NRPS_PKS="type: Hybrid PKS-NRPS' in group,
            'NRPS_PKS="type: PKS' in group,
            'PKS_KS' in group,
            'fatty acid synthase' in group.lower(),
        ]
        if not any(conditions):
            continue
        entry = {
            "pid": gene.group(2),
            "seq": "".join(gene.group(3).split()),
        }
        for line in gene.group(1).split("\n"):
            if "/gene=" in line:
                entry["gene"] = line.split('"')[1]
                break
        results.append(entry)
    return results


def parse_json(fp):
    """Parse MIBIG JSON file.

    Gets organism name, list of compounds encoded by the cluster, as well as any
    pubmed/DOIs corresponding to the entry.
    """
    dic = json.load(fp)
    organism = dic["cluster"]["organism_name"]
    compounds = [entry["compound"] for entry in dic["cluster"]["compounds"]]
    citations = dic["cluster"]["publications"]
    return organism, compounds, citations


def main(gbk_folder, json_folder, output=None, fasta=None, accessions=None):
    # Parse corresponding files
    if not Path(gbk_folder).is_dir() or not Path(json_folder).is_dir():
        raise IOError("Given paths are not directories")
    
    clusters = {}
    
    for gbk in Path(gbk_folder).iterdir():
        accession = gbk.stem
        
        if accessions and accession not in accessions:
            print(f"{accession} not in specified accessions, skipping")
            continue
       
        print(f"Parsing {accession}")
        
        js = Path(json_folder) / f"{accession}.json"
        
        clusters[accession] = {}
        
        with gbk.open() as fp:
            clusters[accession]["synthases"] = parse_gbk(fp)
        
        try:
            with js.open() as fp:
                organism, compounds, citations = parse_json(fp)
                clusters[accession]["organism"] = organism
                clusters[accession]["compounds"] = compounds
                clusters[accession]["citations"] = citations            
        except FileNotFoundError:
            print(f"{accession} has no JSON file")
            
    # Print table
    handle = open(output, mode="w")
    print("Organism\tMIBIG\tAccession\tGene\tLength\tCompounds\tCitations", file=handle)
    for mibig, cluster in clusters.items():
        compounds = ", ".join(cluster["compounds"])
        citations = ", ".join(cluster["citations"])
        for synthase in cluster["synthases"]:
            print(
                cluster["organism"],
                mibig,
                synthase["pid"],
                synthase["gene"] if "gene" in synthase else "",
                len(synthase["seq"]),
                compounds,
                citations,
                sep="|",
                file=handle
            )
    handle.close()

    # Write a FASTA file
    if fasta:
        with open(fasta, "w") as fp:
            for cluster in clusters.values():
                for synthase in cluster["synthases"]:
                    fp.write(f">{synthase['pid']}\n{synthase['seq']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gbk", help="GenBank folder")
    parser.add_argument("json", help="JSON folder")
    parser.add_argument("out", help="Output file")
    parser.add_argument("--fasta", help="Save sequences to FASTA")
    parser.add_argument("--accessions", nargs="+", help="Specific MIBiG accessions")
    args = parser.parse_args()
    main(args.gbk, args.json, output=args.out, fasta=args.fasta, accessions=args.accessions)
