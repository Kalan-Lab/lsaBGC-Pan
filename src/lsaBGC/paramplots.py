"""
paramplots: static (matplotlib) report on how GCF clustering parameters (MCL inflation and Jaccard similarity cutoffs)
influence GCF delineations - produced by lsaBGC-Cluster when run in parameter testing mode.
"""

import math
import random
from collections import defaultdict
from lsaBGC import util
from lsaBGC.phyloplot import MPL_LOCK, _mpl_text, _natural_key

CORE_COLORS = {'Not Applicable': '#737574', 'Core Exists': '#222422', 'Core DNE': '#ed665f'}
COPY_COLORS = {'Multi-Copy': '#000000', 'Single-Copy': '#bebebe'}


def _read_table(path):
	rows = []
	with open(path) as of:
		header = None
		for line in of:
			line = line.rstrip('\n')
			if not line:
				continue
			ls = line.split('\t')
			if header is None:
				header = ls
				continue
			rows.append(dict(zip(header, ls)))
	return rows


def _num(x):
	try:
		v = float(x)
		if math.isnan(v):
			return None
		return v
	except (TypeError, ValueError):
		return None


def _fmt(v):
	"""Format numeric category values the way they are typically written (e.g. 0.8, 20.0)."""
	f = _num(v)
	if f is None:
		return str(v)
	return ('%g' % f) if f != int(f) else ('%.1f' % f)


def _ramp(palette, count):
	return util.brewerColorRamp(palette, count) if count > 0 else []


def _gradient(low, high):
	from matplotlib.colors import LinearSegmentedColormap
	return LinearSegmentedColormap.from_list('grad', [low, high])


def _style(ax):
	for side in ['top', 'right']:
		ax.spines[side].set_visible(False)


def _title_page(fig, text, fontsize):
	ax = fig.add_axes([0, 0, 1, 1])
	ax.axis('off')
	ax.text(0.5, 0.5, _mpl_text(text), ha='center', va='center', fontsize=fontsize, fontweight='bold')
	return ax


def create_param_influence_report(plot_input_1_file, plot_annot_file, plot_input_2_file, plot_overview_file,
								  plot_input_3_file, plot_pdf_file):
	"""
	Create the multi-page parameter impact assessment PDF report for lsaBGC-Cluster.
	********************************************************************************************************************
	Parameters:
	- plot_input_1_file: Per-GCF statistics for each parameter combination tested.
	- plot_annot_file: The BGC annotation classes (one per line, the first being "Unknown").
	- plot_input_2_file: BGC annotation counts for the (up to) 50 largest GCFs for each parameter combination.
	- plot_overview_file: Overview statistics for each parameter combination tested.
	- plot_input_3_file: Counts of samples with single- or multi-copy instances of each GCF.
	- plot_pdf_file: The resulting path to the report PDF file.
	********************************************************************************************************************
	"""
	from matplotlib.figure import Figure
	from matplotlib.backends.backend_pdf import PdfPages
	from matplotlib.lines import Line2D
	import numpy as np

	data = _read_table(plot_input_1_file)
	largest = _read_table(plot_input_2_file)
	overview = _read_table(plot_overview_file)
	copy_counts = _read_table(plot_input_3_file)
	with open(plot_annot_file) as oaf:
		annot_classes = [x.strip() for x in oaf.read().split('\n') if x.strip()]

	# colors for BGC annotation classes - "Unknown" is grey, the rest use a ramp of the ColorBrewer Set2 palette.
	annot_colors = {}
	if annot_classes:
		annot_colors[annot_classes[0]] = '#808080'
		# de-duplicate and never re-color the "Unknown" class (guards against older plotting inputs listing it twice)
		other_classes = [c for c in dict.fromkeys(annot_classes[1:]) if c != annot_classes[0]]
		for c, col in zip(other_classes, _ramp('Set2', len(other_classes))):
			annot_colors[c] = col

	inflations = sorted(set(_fmt(r['Inflation']) for r in overview), key=lambda x: float(x))
	jaccards = sorted(set(_fmt(r['JaccardSim']) for r in overview), key=lambda x: float(x))
	inflation_colors = dict(zip(inflations, _ramp('Set2', len(inflations))))
	markers = ['o', '^', 's', '+', 'x', 'D', 'v', '*', 'P', 'X', 'h', '<', '>', 'p']
	jaccard_markers = {j: markers[i % len(markers)] for i, j in enumerate(jaccards)}

	with MPL_LOCK, PdfPages(plot_pdf_file) as pdf:
		# ---- page 1: title and tradeoff between cluster granularity and singletons
		fig = Figure(figsize=(11, 8))
		tax = fig.add_axes([0, 0.67, 1, 0.33])
		tax.axis('off')
		tax.text(0.5, 0.5, 'Parameter Impact on\nGCF Delineations', ha='center', va='center', fontsize=40,
				 fontweight='bold')
		ax = fig.add_axes([0.1, 0.27, 0.85, 0.36])
		_style(ax)
		rng = random.Random(42)
		xs_all = sorted(set(_num(r['NumberClusters']) for r in overview if _num(r['NumberClusters']) is not None))
		ys_all = sorted(set(_num(r['SingletonToClustered']) for r in overview if _num(r['SingletonToClustered']) is not None))

		def _jitter(vals):
			# jitter overlapping points slightly (at most 40% of the data resolution and 2.5% of the axis range).
			if len(vals) < 2 or vals[-1] - vals[0] <= 0:
				return 0.02 * max(abs(vals[0]), 1.0) if vals else 0.0
			diffs = [b - a for a, b in zip(vals, vals[1:]) if b - a > 0]
			return min(0.4 * min(diffs), 0.025 * (vals[-1] - vals[0]))

		jx, jy = _jitter(xs_all), _jitter(ys_all)
		for r in overview:
			x, y = _num(r['NumberClusters']), _num(r['SingletonToClustered'])
			if x is None or y is None:
				continue
			infl, jac = _fmt(r['Inflation']), _fmt(r['JaccardSim'])
			ax.scatter([x + rng.uniform(-jx, jx)], [y + rng.uniform(-jy, jy)], color=inflation_colors[infl],
					   marker=jaccard_markers[jac], alpha=0.7, s=25)
		ax.set_title('Tradeoff between Cluster Granularity vs. Singletons Incurred', fontsize=11)
		ax.set_xlabel('Number of GCFs')
		ax.set_ylabel('Ratio of BGCs which are\nsingletons to those in GCFs')
		color_handles = [Line2D([0], [0], marker='o', color='none', markerfacecolor=inflation_colors[i],
								markeredgecolor='none', markersize=7) for i in inflations]
		marker_handles = [Line2D([0], [0], marker=jaccard_markers[j], color='#555555', linestyle='none',
								 markersize=7) for j in jaccards]
		leg1 = ax.legend(color_handles, inflations, title='MCL Inflation', loc='upper center',
						 bbox_to_anchor=(0.5, -0.14), ncol=min(len(inflations), 14), frameon=False, fontsize=8,
						 title_fontsize=8)
		ax.add_artist(leg1)
		ax.legend(marker_handles, jaccards, title='Jaccard Similarity Cutoff', loc='upper center',
				  bbox_to_anchor=(0.5, -0.33), ncol=min(len(jaccards), 14), frameon=False, fontsize=8,
				  title_fontsize=8)
		pdf.savefig(fig)

		# ---- page 2: overview heatmaps
		fig = Figure(figsize=(11, 8))
		heat_specs = [['SingletonToClustered', 'Ratio of BGCs which are\nsingletons to those in GCFs',
					   '#add8e6', '#00008b'],
					  ['NumberClusters', 'Number of GCFs', '#90ee90', '#006400'],
					  ['MixedAnnotationProportion', 'Proportion of GCFs which feature\nmultiple BGC annotation '
					   'categories\n(including unannotated)', '#ffff00', '#ffd700'],
					  ['MixedAnnotationWoHypotheticalProportion', 'Proportion of GCFs which feature\nmultiple BGC '
					   'annotation categories\n(not including unannotated)', '#ffc0cb', '#8b0000']]
		for hi, (col, title, low, high) in enumerate(heat_specs):
			ax = fig.add_subplot(2, 2, hi + 1)
			grid = np.full((len(jaccards), len(inflations)), np.nan)
			for r in overview:
				v = _num(r[col])
				if v is None:
					continue
				grid[jaccards.index(_fmt(r['JaccardSim'])), inflations.index(_fmt(r['Inflation']))] = v
			im = ax.imshow(grid, cmap=_gradient(low, high), aspect='auto', origin='lower', interpolation='nearest')
			ax.set_xticks(range(len(inflations)))
			ax.set_xticklabels(inflations, fontsize=7)
			ax.set_yticks(range(len(jaccards)))
			ax.set_yticklabels(jaccards, fontsize=7)
			ax.set_xlabel('MCL Inflation', fontsize=8)
			ax.set_ylabel('Jaccard Similarity Cutoff', fontsize=8)
			ax.set_title(title, fontsize=9)
			cb = fig.colorbar(im, ax=ax)
			cb.ax.tick_params(labelsize=7)
		fig.tight_layout()
		pdf.savefig(fig)

		# ---- page 3: boxplots of per-GCF statistics faceted by Jaccard similarity cutoff
		fig = Figure(figsize=(11, 8))
		box_specs = [['CoreSize', 'Number of Core Homolog Groups per GCF'],
					 ['Samples', 'Number of Samples with GCF'],
					 ['AvgGeneCount', 'Avg. Gene Count per BGC in GCF'],
					 ['StdDevGeneCount', 'Std. Deviation of Gene Count for BGCs in GCF']]
		gs = fig.add_gridspec(len(box_specs), max(1, len(jaccards)))
		for bi, (col, title) in enumerate(box_specs):
			for ji, jac in enumerate(jaccards):
				ax = fig.add_subplot(gs[bi, ji])
				_style(ax)
				series = []
				for infl in inflations:
					vals = [_num(r[col]) for r in data if _fmt(r['JaccardSim']) == jac and _fmt(r['Inflation']) == infl]
					series.append([v for v in vals if v is not None])
				positions = [i for i, s in enumerate(series) if s]
				if positions:
					bp = ax.boxplot([series[i] for i in positions], positions=positions, widths=0.6,
									patch_artist=True, flierprops=dict(markersize=2))
					for patch, i in zip(bp['boxes'], positions):
						patch.set_facecolor(inflation_colors[inflations[i]])
					for med in bp['medians']:
						med.set_color('#000000')
				ax.set_xticks(range(len(inflations)))
				ax.set_xticklabels(inflations if bi == len(box_specs) - 1 else [''] * len(inflations), fontsize=6)
				ax.set_xlim(-0.6, len(inflations) - 0.4)
				ax.tick_params(axis='y', labelsize=6)
				if bi == 0:
					ax.set_title('Jaccard Similarity Cutoff: %s' % jac, fontsize=7)
				if ji == 0:
					ax.set_ylabel(title, fontsize=7)
				if bi == len(box_specs) - 1:
					ax.set_xlabel('MCL Inflation', fontsize=7)
		fig.tight_layout()
		pdf.savefig(fig)

		# ---- one page per parameter combination
		params = []
		for r in data:
			if r['Parameters'] not in params:
				params.append(r['Parameters'])
		for param in params:
			data_filt = [r for r in data if r['Parameters'] == param]
			largest_filt = [r for r in largest if r['Parameters'] == param]
			copy_filt = [r for r in copy_counts if r['Parameters'] == param]

			fig = Figure(figsize=(11, 8))
			gs = fig.add_gridspec(4, 2, height_ratios=[1, 2, 2, 2])
			tax = fig.add_subplot(gs[0, :])
			tax.axis('off')
			tax.text(0.5, 0.5, _mpl_text(param), ha='center', va='center', fontsize=18, fontweight='bold')

			# histogram of sample counts per GCF, stacked by whether a core exists
			ax = fig.add_subplot(gs[1, 0])
			_style(ax)
			cats = [c for c in CORE_COLORS if any(r['CoreExists'] == c for r in data_filt)]
			series = [[_num(r['Samples']) for r in data_filt if r['CoreExists'] == c and _num(r['Samples']) is not None]
					  for c in cats]
			if any(series):
				ax.hist(series, bins=30, stacked=True, color=[CORE_COLORS[c] for c in cats], edgecolor='black',
						linewidth=0.3, label=cats)
				ax.legend(frameon=False, fontsize=7, loc='best')
			ax.set_xlabel('Sample Count', fontsize=8)
			ax.set_ylabel('count', fontsize=8)
			ax.tick_params(labelsize=7)

			# histogram of standard deviation in gene count
			ax = fig.add_subplot(gs[1, 1])
			_style(ax)
			vals = [_num(r['StdDevGeneCount']) for r in data_filt]
			vals = [v for v in vals if v is not None]
			if vals:
				ax.hist(vals, bins=30, color='black', edgecolor='black')
				ax.set_yscale('log')
			ax.set_xlabel('Std. Deviation in Gene Count', fontsize=8)
			ax.set_ylabel('count', fontsize=8)
			ax.tick_params(labelsize=7)

			# stacked bars of samples with single/multi-copy instances of each GCF, ordered by total BGC count
			ax = fig.add_subplot(gs[2, :])
			_style(ax)
			gcf_totals, gcf_counts = {}, defaultdict(dict)
			for r in copy_filt:
				gcf_totals[r['GCF']] = _num(r['Total.BGCs']) or 0.0
				gcf_counts[r['GCF']][r['Category']] = _num(r['Sample_Count']) or 0.0
			order = sorted(gcf_totals, key=lambda g: (gcf_totals[g], _natural_key(g)))
			bottom = np.zeros(len(order))
			for cat in sorted(COPY_COLORS):
				heights = np.array([gcf_counts[g].get(cat, 0.0) for g in order])
				if len(order):
					ax.bar(range(len(order)), heights, bottom=bottom, color=COPY_COLORS[cat], edgecolor='black',
						   linewidth=0.3, width=0.9)
				bottom += heights
			ax.set_xticks([])
			ax.set_xlim(-0.6, max(len(order), 1) - 0.4)
			ax.set_xlabel('All GCFs (grey=single-copy, black=multi-copy)', fontsize=8)
			ax.set_ylabel('Sample Count', fontsize=8)
			ax.tick_params(labelsize=7)

			# stacked bars of BGC annotation classes for the largest GCFs
			ax = fig.add_subplot(gs[3, :])
			_style(ax)
			gcf_totals, gcf_annots = {}, defaultdict(dict)
			for r in largest_filt:
				gcf_totals[r['GCF']] = _num(r['Total.BGCs']) or 0.0
				# accumulate rather than overwrite, in case an annotation class appears more than once for a GCF
				gcf_annots[r['GCF']][r['Annotation']] = gcf_annots[r['GCF']].get(r['Annotation'], 0.0) + (_num(r['Count']) or 0.0)
			order = sorted(gcf_totals, key=lambda g: (gcf_totals[g], _natural_key(g)))
			bottom = np.zeros(len(order))
			for annot in sorted(set(a for g in gcf_annots for a in gcf_annots[g])):
				heights = np.array([gcf_annots[g].get(annot, 0.0) for g in order])
				if len(order):
					ax.bar(range(len(order)), heights, bottom=bottom, color=annot_colors.get(annot, '#808080'),
						   edgecolor='black', linewidth=0.3, width=0.9)
				bottom += heights
			ax.set_xticks([])
			ax.set_xlim(-0.6, max(len(order), 1) - 0.4)
			ax.set_xlabel('Top 50 GCFs with Most BGCs (color = BGC product type)', fontsize=8)
			ax.set_ylabel('BGC Count', fontsize=8)
			ax.tick_params(labelsize=7)
			fig.tight_layout()
			pdf.savefig(fig)
