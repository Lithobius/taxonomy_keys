#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
based on https://github.com/gjeunen/reference_database_creator/blob/main/function/crabs_functions.py
Copyright for original functions to gjeunen

I'm having trouble with this figure and in the future I want to update taxonomy for it.

Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev

"""

# module import
# from original; I don't know if I need all these?
import requests, tarfile, rich, os, zipfile, shutil, ftplib, fnmatch, collections, tempfile, time
import rich.progress
import rich_click as click
import subprocess as sp
from rich.progress import Progress, BarColumn, TextColumn
from matplotlib import pyplot as plt
import numpy as np

from pyprojroot import here # like here in R

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'primerefficiency-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


def amplicon_import(amplicons, tax_group):
    '''
    takes an input file and reads it into memory
    '''
    amplicons_dict = {}
    with open(amplicons, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            lineparts = line.split('\t')
            if tax_group:
                if tax_group in line:
                    amplicons_dict[lineparts[0]] = lineparts[-1]
            else:
                amplicons_dict[lineparts[0]] = lineparts[-1]
    return amplicons_dict

def raw_import(input, amplicons_dict):
    '''
    takes an input file and reads it into memory
    '''
    raw_dict = {}
    with open(input, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            lineparts = line.split('\t')
            if lineparts[0] in amplicons_dict:
                raw_dict[lineparts[0]] = lineparts[-1]
    return raw_dict

def rev_comp(oligo):
    '''
    takes in an oligo nucleotide string and returns the reverse complement
    '''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
    return ''.join(complement[base] for base in oligo[::-1])

def extract_primer_regions(amplicons_dict, raw_dict, forward, reverse):
    '''
    extract primer binding regions from barcodes and returns a dict
    '''
    primer_binding_region_dict = collections.defaultdict(dict)
    for item in amplicons_dict:
        forward_region = ''
        reverse_region = ''
        try:
            raw_seq = raw_dict[item]
            barcode = amplicons_dict[item]
            if raw_seq.find(barcode) != -1:
                logging.info(raw_seq.find(barcode))
                if raw_seq.find(barcode) == 0:
                    logging.info(raw_seq.find(forward))
                    forward_region = raw_seq[raw_seq.find(barcode) : raw_seq.find(barcode) + len(forward)]
                    reverse_region = rev_comp(raw_seq[len(barcode) - len(reverse) : len(barcode)])
                    #logging.info(raw_seq[len(barcode) : len(barcode) + len(reverse)])
                    #logging.info(reverse_region)
                    #logging.info(rev_comp(forward_region))
                else:
                    forward_region = raw_seq[raw_seq.find(barcode) - len(forward) : raw_seq.find(barcode)]
                    reverse_region = rev_comp(raw_seq[raw_seq.find(barcode) + len(barcode) : raw_seq.find(barcode) + len(barcode) + len(reverse)])
                #logging.info(forward_region + ' ' + str(raw_seq.find(barcode) + len(barcode) + len(reverse)))
            else:
                rev_raw_seq = rev_comp(raw_seq)
                if rev_raw_seq.find(barcode) != -1:
                    logging.info('reverse complement')
                    if rev_raw_seq.find(barcode) == 0:
                        forward_region = raw_seq[rev_raw_seq.find(barcode) : rev_raw_seq.find(barcode) + len(forward)]
                        reverse_region = rev_comp(rev_raw_seq[rev_raw_seq.find(barcode) +len(forward) + len(barcode) : rev_raw_seq.find(barcode) + len(barcode) +len(forward) + len(reverse)])
                    else:
                        forward_region = rev_raw_seq[rev_raw_seq.find(barcode) - len(forward) : rev_raw_seq.find(barcode)]
                        reverse_region = rev_comp(rev_raw_seq[rev_raw_seq.find(barcode) + len(barcode) : rev_raw_seq.find(barcode) + len(barcode) + len(reverse)])
                    #logging.info(reverse_region)
        except KeyError:
            continue
        if len(forward_region) == len(forward) and len(reverse_region) == len(reverse):
            primer_binding_region_dict[item]['forward'] = forward_region
            primer_binding_region_dict[item]['reverse'] = reverse_region
        logging.info(forward_region + ' ' + reverse_region)
    return primer_binding_region_dict

def deconstruct_primer_regions(primer_dict, key):
    '''
    deconstructs string and places them in a dictionary of lists
    '''
    position_dict = collections.defaultdict(list)
    for item in primer_dict:
        for i in range(len(primer_dict[item][key])):
            position_dict[i].append(primer_dict[item][key][i])
    return position_dict

def dict_to_array(position_dict):
    '''
    takes a dict and returns an np.array
    '''
    positions = []
    ordered_counts = []
    for position in position_dict:
        sequence = position_dict[position]
        counts = {
            'A': sequence.count('A') / len(sequence) * 100,
            'C': sequence.count('C') / len(sequence) * 100,
            'G': sequence.count('G') / len(sequence) * 100,
            'T': sequence.count('T') / len(sequence) * 100,
        }
        counts['Other'] = 100 - sum(counts.values())
        sorted_counts = sorted(counts.items(), key = lambda x: x[1], reverse = True)
        positions.append(position)
        ordered_counts.append(sorted_counts)
    positions = np.array(positions)
    bottoms = np.zeros(len(positions))
    return positions, ordered_counts, bottoms


def parse_primer(primer):
    '''
    parses a primer sequence for plotting
    '''
    ordered_counts = []
    for sequence in primer:
        counts = {
            'A': sequence.count('A') / len(sequence) * 100,
            'C': sequence.count('C') / len(sequence) * 100,
            'G': sequence.count('G') / len(sequence) * 100,
            'T': sequence.count('T') / len(sequence) * 100,
        }
        counts['Other'] = 100 - sum(counts.values())
        sorted_counts = sorted(counts.items(), key = lambda x: x[1], reverse = True)
        ordered_counts.append(sorted_counts)
    logging.info('parsed ' + primer)
    return ordered_counts

# make plot
def efficiency_barplot(forward_positions, forward_ordered_counts, forward_bottoms, reverse_positions, reverse_ordered_counts, reverse_bottoms, forward_primer_info, reverse_primer_info, forward, reverse, output):
    '''
    generates a bar plot with base compositions indicating amplification efficiency
    '''
    width = 0.8
    fig, axs = plt.subplots(2, 2, gridspec_kw = {'height_ratios': [20, 1]})
    # forward primer-binding region subfigure
    for i in range(5):
        forward_counts = [forward_ordered_counts[j][i][1] for j in range(len(forward_ordered_counts))]
        forward_labels = [forward_ordered_counts[j][i][0] for j in range(len(forward_ordered_counts))]
        reverse_counts = [reverse_ordered_counts[j][i][1] for j in range(len(reverse_ordered_counts))]
        reverse_labels = [reverse_ordered_counts[j][i][0] for j in range(len(reverse_ordered_counts))]
        fprimer_counts = [forward_primer_info[j][i][1] for j in range(len(forward_primer_info))]
        fprimer_labels = [forward_primer_info[j][i][0] for j in range(len(forward_primer_info))]
        rprimer_counts = [reverse_primer_info[j][i][1] for j in range(len(reverse_primer_info))]
        rprimer_labels = [reverse_primer_info[j][i][0] for j in range(len(reverse_primer_info))]
        colors = {'A': '#e09f3e', 'C': '#335c67', 'G': '#fff3b0', 'T': '#9e2a2b', 'Other': 'gray'}
        forward_color = [colors[label] for label in forward_labels]
        reverse_color = [colors[label] for label in reverse_labels]
        fprimer_color = [colors[label] for label in fprimer_labels]
        rprimer_color = [colors[label] for label in rprimer_labels]
        logging.info(forward_counts)
        axs[0, 0].bar(forward_positions, forward_counts, width = width, bottom = forward_bottoms, color = forward_color, label = forward_labels[0] if i == 0 else "_nolegend_")
        axs[0, 1].bar(reverse_positions, reverse_counts, width = width, bottom = reverse_bottoms, color = reverse_color, label = reverse_labels[0] if i == 0 else "_nolegend_")
        axs[1, 0].bar(forward_positions, fprimer_counts, width = width, bottom = forward_bottoms, color = fprimer_color, label = fprimer_labels[0] if i == 0 else "_nolegend_")
        axs[1, 1].bar(reverse_positions, rprimer_counts, width = width, bottom = reverse_bottoms, color = rprimer_color, label = rprimer_labels[0] if i == 0 else "_nolegend_")
        forward_bottoms += forward_counts
        reverse_bottoms += reverse_counts
    handles = [plt.Line2D([0], [0], color = '#e09f3e', lw = 4, label = 'A'),
               plt.Line2D([0], [0], color = '#335c67', lw = 4, label = 'C'),
               plt.Line2D([0], [0], color = '#fff3b0', lw = 4, label = 'G'),
               plt.Line2D([0], [0], color = '#9e2a2b', lw = 4, label = 'T'),
               plt.Line2D([0], [0], color = 'gray', lw = 4, label = 'Other'),]
    axs[0,0].set_ylabel('Proportion of bp occurrences')
    axs[0,0].tick_params(bottom=False)
    axs[1,0].tick_params(left=False)
    axs[1,0].tick_params(bottom=False)
    axs[1,0].set(yticklabels=[])
    axs[1,0].set(xticklabels=[])
    axs[0,0].set(xticklabels=[])
    axs[1,0].margins(y=0)
    axs[0,0].margins(y=0)
    axs[0,0].set_title('Forward primer')
    axs[0,1].legend(bbox_to_anchor=(1.0, 1.0), handles = handles, title = 'Nucleotide')

    axs[0,1].tick_params(left=False)
    axs[0,1].tick_params(bottom=False)
    axs[0,1].set(yticklabels=[])
    axs[0,1].set(xticklabels=[])
    axs[1,1].margins(y=0)
    axs[0,1].margins(y=0)
    axs[0,1].set_title('Reverse primer')
    axs[1,1].set(yticklabels=[])
    axs[1,1].tick_params(left=False)
    axs[1,1].tick_params(bottom=False)
    axs[1,1].set(xticklabels=[])
    for x, p in zip(reverse_positions, list(reverse)):
        axs[1, 1].text(x, 50, p, color = 'black', ha = 'center', va = 'center')
    for x, p in zip(forward_positions, list(forward)):
        axs[1, 0].text(x, 50, p, color = 'black', ha = 'center', va = 'center')
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.05, wspace = 0.05)
    plt.savefig(output)

# basic-body

if __name__ == "__main__":

    crabs_pcr_out = here("./database12s/crabs/20241015_crabs3_pcrout.txt")
    #selected_group  = 'Salmonidae'
    selected_group = 'Clupeidae'
    # make a dictionary of amplicons from the in silico PCR, filtered to a taxonomy level
    amplicons_dict = amplicon_import(crabs_pcr_out, selected_group)
    merged_crabs = here("./database12s/crabs/20241015_12s_crabs-merge.txt")
    # dictionary with all the sequences
    raw_dict = raw_import(merged_crabs, amplicons_dict)
    # Primer info
    mifish_fw = 'GTCGGTAAAACTCGTGCCAGC'
    mifish_rv = 'CATAGTGGGGTATCTAATCCCAGTTTG'
    # extract primer binding regions
    primer_binding_region_dictionary = extract_primer_regions(amplicons_dict, raw_dict, mifish_fw, mifish_rv)
    # calculate base proportions
    forward_position_dict = deconstruct_primer_regions(primer_binding_region_dictionary, 'forward')
    reverse_position_dict = deconstruct_primer_regions(primer_binding_region_dictionary, 'reverse')
    # transform dict to np.array
    forward_positions, forward_ordered_counts, forward_bottoms = dict_to_array(forward_position_dict)
    reverse_positions, reverse_ordered_counts, reverse_bottoms = dict_to_array(reverse_position_dict)
    # parse primers to plot
    forward_primer_info = parse_primer(mifish_fw)
    reverse_primer_info = parse_primer(mifish_rv)
    # figure
    exportfile = here("./database12s/crabs/plots/20241015_primerefficiency.png")
    efficiency_barplot(forward_positions, forward_ordered_counts, forward_bottoms, reverse_positions, reverse_ordered_counts, reverse_bottoms, forward_primer_info, reverse_primer_info, mifish_fw, mifish_rv, exportfile)
    #logging.info(primer_binding_region_dictionary)

# write code here
# don't forget logging.debug('message' or object) periodically


# end
logging.info('Fin!')
