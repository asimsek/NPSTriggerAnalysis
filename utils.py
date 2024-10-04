import pandas as pd
pd.set_option('display.max_rows', 800)
import re, os
import argparse
## import tsgauth
import requests
import urllib3
import json 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import awkward as ak

import uproot
import mplhep as hep
hep.style.use('CMS')
import matplotlib.style as mplstyle
mplstyle.use('fast')

from scipy.ndimage import median_filter, uniform_filter
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages

import itertools
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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
	"""
	Adds vertical lines at the end of each era (except the last one) and labels on the upper x-axis.

	Parameters:
	- ax: matplotlib axis object for plotting the vertical lines.
	- ax_upper: matplotlib secondary axis object for labeling the eras.
	- eras: dictionary containing the eras and their respective ranges.
	- x_range: tuple containing the min and max values for the x-axis data.
	"""
	min_x_, max_x_ = x_range
	min_x, max_x = int(min_x_), int(max_x_)

	era_keys = list(eras.keys())
	tick_positions = []
	era_labels = []

	for idx, era in enumerate(era_keys):
		start_, end_ = eras[era]
		start, end = int(start_), int(end_)

		# Ignore eras that are completely outside the range of the data
		if start > max_x or end < min_x:
			continue

		# Adjust the range to be within the data limits
		if start < min_x:
			start = min_x
		if end > max_x:
			end = max_x

		# Add a vertical line at the end of the era (except for the last era)
		if idx < len(era_keys) - 1:
			ax.axvline(x=end, color='k', linestyle='--', linewidth=0.8)

		# Set tick positions and labels for the upper x-axis
		tick_positions.append((start + end) / 2)
		era_labels.append(era)

	# Set ticks and labels only if there are any valid ones
	if tick_positions:
		ax_upper.set_xticks(tick_positions)
		ax_upper.set_xticklabels(era_labels, rotation=45, ha='center')

