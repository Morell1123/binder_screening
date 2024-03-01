import os
import pandas as pd
import json
import numpy as np
from Bio import SeqIO
import argparse
import shutil

"""
This script extracts the AlphaFold (AF) metrics of ipAE and pLDDT 
from the output of an AF run on binders and targets (multimers)

The Script need to input fasta files used for alpha fold. 
For this script to work the binding complexes (multimers) in the .fasta file input to alphafold must be of the form:
Fasta file format:
-------------------------------------------------------------
>binder_{NAME_OF_BINDER}_target_{NAME_OF_target}
{BINDER_SEQUENCE}:{TARGET_SEQUENCE}
--------------------------------------------------------------
The binder sequence must come first in the fasta file entries, as the script uses the length of the sequences 
(or the location of the split between binder and target) to compute the pLDDT binder and ipAE from the AF output
"""


def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]

def get_pae_interaction(data, complex_seq):
    pae = data.get('pae', [])
    target_index = complex_seq.find(':') # FIXME: this assumes that the binder is never a multimer, and is always written first in the fasta file
    pae_binder_to_target = [item for sublist in pae[0:target_index] for item in sublist[target_index:]]
    pae_target_to_binder = [item for sublist in pae[target_index:] for item in sublist[0:target_index]]
    return np.mean(pae_binder_to_target + pae_target_to_binder)

def get_mean_pae_overall(data):
    pae = data.get('pae', [])
    flat_pae = flatten_list(pae)
    return np.mean(flat_pae)


def get_mean_plddt(data):
    plddt = data.get('plddt', [])
    return np.mean(plddt)

def get_mean_plddt_binder(data, complex_seq):
    plddt = data.get('plddt', [])
    target_index = complex_seq.find(':')
    plddt_binder = plddt[:target_index]
    return np.mean(plddt_binder)


def Extract_AF_metrics(output_folder_path, sequences):
    if not os.path.isdir(output_folder_path):
        os.mkdir(output_folder_path)
    data = []
    for filename in os.listdir(output_folder_path):
        if 'rank_001' in filename and os.path.splitext(filename)[1] == '.json':
            complex_name = filename.split('_scores')[0]
            file_path = os.path.join(output_folder_path, filename)
            with open(file_path, 'r') as file:
                json_data = json.load(file)
                # mean_plddt = get_mean_plddt(json_data)
                complex_seq = sequences[complex_name]
                mean_pae_inter = get_pae_interaction(json_data, complex_seq)
                # mean_pae_overall = get_mean_pae_overall(json_data)
                mean_plddt_binder = get_mean_plddt_binder(json_data, complex_seq)


            data.append({
                "Filename": filename[: filename.index('_scores')],
                "IpAE": mean_pae_inter,
                "binder_pLDDT": mean_plddt_binder
            })

    df_metrics = pd.DataFrame(data)

    return df_metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process inputs')

    parser.add_argument('--colab_input', type=str, help='The original fasta file inputted to ColabFold')
    parser.add_argument('--colab_output_folder', type=str, help='The folder in which the ColabFold output is located')
    parser.add_argument('--output_csv', type=str, help='File path of output (the extracted AF2 metrics)')
    parser.add_argument('--save_structures', type=bool, default=True, help='Save structures or keep delete them to be more memory efficient')
    parser.add_argument('--save_raw_metrics_data', type=bool, default=True, help='Save structures or keep delete them to be more memory efficient')
    parser.add_argument('--save_all_outdata', type=bool, default=False, help='Save structures or keep delete them to be more memory efficient')

    args = parser.parse_args()
    output_dir = args.colab_output_folder

    complexes_path = args.colab_input
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(complexes_path, "fasta")}

    metrics_df = Extract_AF_metrics(output_dir, sequences)
    if not args.save_all_outdata:
        for filename in os.listdir(output_dir):
            file = os.path.join(output_dir,filename)
            if os.path.isdir(file):
                shutil.rmtree(file)
                continue
            elif os.path.splitext(file)[1] == '.pdb' and args.save_structures and 'rank_001' in file:
                continue
            elif os.path.splitext(file)[1] == '.json' and args.save_raw_metrics_data and 'rank_001' in file:
                continue
            os.remove(file)
    try:
        os.rmdir(output_dir)
    except:
        pass

    assert len(metrics_df) != 0, "Your metrics dataframe is empty. Please check you have both an input fasta file and files in the output, which matchES"
    metrics_df = metrics_df[metrics_df['Filename'].str.contains('_target_')]
    assert metrics_df['Filename'].str.contains('_target_').all(), "Error: Your input format is wrong. Not all your filenames in the .fasta input file contains '_target_'"
    assert (metrics_df['Filename'].str.find('binder_') == 0).all(), "Error: Your input format is wrong. Not all your filenames in the .fasta input file start with 'binder_' "

    metrics_df['binder'] = metrics_df['Filename'].apply(lambda x: x[x.index('binder')+7: x.index('_target')])
    metrics_df['target'] = metrics_df['Filename'].apply(lambda x: x[x.index('target')+7:])

    metrics_df = metrics_df.pivot(index='binder', columns='target', values=['IpAE', 'binder_pLDDT'])
    metrics_df.columns = ['_'.join(col).strip() for col in metrics_df.columns.values]

    metrics_df.to_csv(args.output_csv, sep=';')