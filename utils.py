## Standard libraries
import os
import re
import argparse
import json
import datetime
import itertools
import warnings

##
import uproot
import numpy as np
import pandas as pd
import awkward as ak
import requests
import urllib3
#import tsgauth
import awkward as ak
import mplhep as hep
from prettytable import PrettyTable

# SciPy 
from scipy.ndimage import median_filter, uniform_filter
from scipy.signal import savgol_filter

# Matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.dates as mdates
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

mplstyle.use('fast')
hep.style.use('CMS')
pd.set_option('display.max_rows', 800)

warnings.simplefilter("ignore")

# Suppress warnings related to numpy subnormal value issues
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in divide")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in divide")
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float32'> type is zero.")

LS_seconds = 2**18 / 11245.5 # copied from Silvios OMS rates ntuples 

def multipage(filename, figs=None, dpi=200):
	"""Creates a pdf with one page per plot"""
	pp = PdfPages(filename)
	if figs is None:
		figs = [plt.figure(n) for n in plt.get_fignums()]
		print("No figures handed")
	for fig in figs:
		plt.figure(fig).savefig(pp, format='pdf',bbox_inches='tight')
	pp.close()

def add_era_lines(ax, ax_upper, eras, x_range):
    min_x_, max_x_ = x_range
    if isinstance(min_x_, pd.Timestamp):
        min_x, max_x = min_x_.to_pydatetime(), max_x_.to_pydatetime()
    else:
        min_x, max_x = float(min_x_), float(max_x_)

    era_keys = list(eras.keys())
    tick_positions = []
    era_labels = []

    for idx, era in enumerate(era_keys):
        start_, end_ = eras[era][0], eras[era][1]
        try:
            if isinstance(start_, pd.Timestamp) and isinstance(end_, pd.Timestamp):
                start, end = start_, end_
            else:
                start, end = int(start_), int(end_)
        except ValueError:
            try:
                start, end = pd.to_datetime(start_), pd.to_datetime(end_)
            except Exception as e:
                print(f"Error converting start/end for era {era}: {e}")
                continue

        # Ensure both start and end are either int or datetime for comparison
        if isinstance(start, str) or isinstance(end, str):
            print(f"Skipping era {era} due to incompatible types: start={start}, end={end}")
            continue

        # Ignore eras that are completely outside the range of the data
        if isinstance(start, (int, float)) and isinstance(end, (int, float)):
            if start > max_x or end < min_x:
                continue
        elif isinstance(start, datetime.datetime) and isinstance(end, datetime.datetime):
            if start > max_x or end < min_x:
                continue

        # Adjust the range to be within the data limits
        if isinstance(start, (int, float)) and start < min_x:
            start = min_x
        if isinstance(end, (int, float)) and end > max_x:
            end = max_x
        if isinstance(start, datetime.datetime) and start < min_x:
            start = min_x
        if isinstance(end, datetime.datetime) and end > max_x:
            end = max_x

        # Add a vertical line at the end of the era (except for the last era)
        if idx < len(era_keys) - 1:
            ax.axvline(x=end, color='k', linestyle='--', linewidth=0.8)

        # Set tick positions and labels for the upper x-axis
        if isinstance(start, datetime.datetime) and isinstance(end, datetime.datetime):
            tick_positions.append(start + (end - start) / 2)
        elif isinstance(start, (int, float)) and isinstance(end, (int, float)):
            tick_positions.append((start + end) / 2)
        era_labels.append(era)

    # Set ticks and labels only if there are any valid ones
    if tick_positions:
        ax_upper.set_xticks(tick_positions)
        ax_upper.set_xticklabels(
            era_labels,
            rotation=0,
            ha='center',
            va='bottom',  # Vertical alignment
            fontsize=9,   # Adjust font size if needed
            transform=ax_upper.get_xaxis_transform() + ax_upper.transAxes
        )

