"""
phyloplot: phylogenetic heatmaps for lsaBGC-Pan, rendered both as interactive HTML (Plotly) and static PDF (matplotlib).

Each figure shows a (midpoint-rooted) phylogeny on the left with aligned leaf labels, optionally colored tip markers,
followed by one or more categorical heatmaps and an optional gene-architecture (arrow) panel, all sharing the same rows.
Figure dimensions scale with the number of leaves and columns, so large trees are not squeezed onto a single page; in
the browser, users can zoom/pan, hover for details, and download an SVG via the camera icon.

The tool-specific functions (e.g. see_html / see_pdf) build one figure specification from the track files produced by
the lsaBGC programs, which is then rendered by either backend, so the HTML and PDF versions always show the same data.

plotly.js is loaded from a CDN, so an internet connection is needed when opening the resulting HTML files.
"""

import os
import sys
import math
import threading
import traceback
from collections import defaultdict
from Bio import Phylo

DEFAULT_PALETTE = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
				   "#9c755f", "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf",
				   "#bcbd22", "#aec7e8", "#ffbb78", "#98df8a"]

HELP_TEXT = ("Drag to zoom, double-click to reset, hover for details; click legend entries to hide/show them; "
			 "use the camera icon to download an SVG.")

# matplotlib's text/font caches are not guaranteed to be thread-safe and several lsaBGC programs create figures from
# thread pools, so static figure rendering is serialized.
MPL_LOCK = threading.Lock()

# Most PDF viewers cannot display pages larger than 200 x 200 inches.
MAX_PDF_INCHES = 199.0


def _assign_colors(categories, fixed=None):
	"""Assign colors to categories (sorted naturally), respecting any fixed color assignments."""
	fixed = fixed or {}
	colors = {}
	i = 0
	for c in sorted(set(categories), key=_natural_key):
		if c in fixed:
			colors[c] = fixed[c]
		else:
			colors[c] = DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)]
			i += 1
	return colors


def _natural_key(s):
	import re
	return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]


def _nice_ticks(max_value, n=4):
	if max_value <= 0:
		return [0]
	return [float('%.3g' % (max_value * i / n)) for i in range(n + 1)]


def tree_layout(tree_file, midpoint_root=True):
	"""
	Compute a rectangular layout for a phylogeny.
	********************************************************************************************************************
	Returns:
	- names: leaf names ordered top to bottom (row i corresponds to y = i).
	- seg_x, seg_y: coordinates (with None separators) for drawing tree branches.
	- leaf_x: dictionary of leaf name -> x coordinate (distance from root).
	- max_depth: maximum distance from root.
	********************************************************************************************************************
	"""
	tree = Phylo.read(tree_file, 'newick')
	if midpoint_root and len(tree.get_terminals()) > 2:
		try:
			tree.root_at_midpoint()
		except Exception:
			pass
	tree.ladderize()

	depths = tree.depths()
	if max(depths.values()) <= 0:
		depths = tree.depths(unit_branch_lengths=True)

	leaves = tree.get_terminals()
	y = {}
	for i, leaf in enumerate(leaves):
		y[leaf] = float(i)
	for clade in tree.find_clades(order='postorder'):
		if not clade.is_terminal():
			child_ys = [y[c] for c in clade.clades]
			y[clade] = (min(child_ys) + max(child_ys)) / 2.0

	seg_x, seg_y = [], []
	for clade in tree.find_clades(order='preorder'):
		if clade.is_terminal():
			continue
		x_parent = depths[clade]
		child_ys = [y[c] for c in clade.clades]
		seg_x += [x_parent, x_parent, None]
		seg_y += [min(child_ys), max(child_ys), None]
		for c in clade.clades:
			seg_x += [x_parent, depths[c], None]
			seg_y += [y[c], y[c], None]

	names = [str(leaf.name) for leaf in leaves]
	leaf_x = {str(leaf.name): depths[leaf] for leaf in leaves}
	return names, seg_x, seg_y, leaf_x, max(depths.values())


def _heatmap_categories(h, names):
	"""Determine the (naturally sorted) categories present in a heatmap and their colors."""
	categories = set()
	for leaf in names:
		for col in h['columns']:
			v = h['values'].get(leaf, {}).get(col)
			if v is not None:
				categories.add(str(v))
	colors = _assign_colors(categories, h.get('colors'))
	return sorted(colors, key=_natural_key), colors


def create_phylo_heatmap_html(tree_file, output_html, title, heatmaps=None, gene_panel=None, tip_annotation=None,
							  midpoint_root=True, row_height=16):
	"""
	Create an interactive HTML phylogenetic heatmap.
	********************************************************************************************************************
	Parameters:
	- tree_file: Path to phylogeny in Newick format.
	- output_html: Path to the resulting HTML file.
	- title: Title of the figure.
	- heatmaps: List of dictionaries, each describing a categorical heatmap with keys:
	    - 'title': panel title.
	    - 'columns': ordered list of column names.
	    - 'values': dictionary of leaf -> {column: category}; missing cells are left blank.
	    - 'colors' (optional): dictionary of category -> color (others assigned automatically).
	    - 'show_legend' (optional, default True): whether to list categories in the legend.
	    - 'category_labels' (optional): dictionary of category -> label used in legend/hover.
	- gene_panel: Optional dictionary with keys 'title' and 'genes': leaf -> list of dictionaries with keys 'start',
	              'end', 'forward' (bool), 'name', 'group' and 'color'.
	- tip_annotation: Optional dictionary with keys 'title', 'values' (leaf -> category) and optionally 'colors'.
	- midpoint_root: Whether to midpoint root the phylogeny.
	- row_height: Height in pixels per leaf/row.
	********************************************************************************************************************
	"""
	import plotly.graph_objects as go
	from plotly.subplots import make_subplots

	heatmaps = [h for h in (heatmaps or []) if len(h.get('columns', [])) > 0]
	names, seg_x, seg_y, leaf_x, max_depth = tree_layout(tree_file, midpoint_root=midpoint_root)
	n = len(names)
	if n == 0:
		raise RuntimeError('No leaves found in phylogeny: %s' % tree_file)
	if max_depth <= 0:
		max_depth = 1.0

	# determine panel widths in pixels so that the figure grows with its content
	max_label_len = max(len(x) for x in names)
	label_px = min(450, 7 * max_label_len + 30)
	tree_core_px = 320
	panel_px = [tree_core_px + label_px]
	for h in heatmaps:
		panel_px.append(min(1600, max(140, 14 * len(h['columns']))))
	if gene_panel:
		panel_px.append(800)
	gap_px = 25
	plot_width = sum(panel_px) + gap_px * (len(panel_px) - 1)
	margin = dict(l=30, r=260, t=110, b=170)
	fig_width = plot_width + margin['l'] + margin['r']
	fig_height = max(450, n * row_height + margin['t'] + margin['b'])

	subplot_titles = ['Phylogeny'] + [h.get('title', '') for h in heatmaps]
	if gene_panel:
		subplot_titles.append(gene_panel.get('title', 'Gene architecture'))
	fig = make_subplots(rows=1, cols=len(panel_px), shared_yaxes=True,
						column_widths=[p / float(sum(panel_px)) for p in panel_px],
						horizontal_spacing=gap_px / float(plot_width), subplot_titles=subplot_titles)

	# ---- phylogeny panel
	label_x = max_depth * 1.03
	data_per_px = label_x / float(tree_core_px)
	fig.add_trace(go.Scatter(x=seg_x, y=seg_y, mode='lines', line=dict(color='#444444', width=1.2),
							 hoverinfo='skip', showlegend=False), row=1, col=1)
	guide_x, guide_y = [], []
	for i, name in enumerate(names):
		if leaf_x[name] < max_depth * 1.01:
			guide_x += [leaf_x[name], label_x * 0.995, None]
			guide_y += [i, i, None]
	fig.add_trace(go.Scatter(x=guide_x, y=guide_y, mode='lines', line=dict(color='#c8c8c8', width=0.7, dash='dot'),
							 hoverinfo='skip', showlegend=False), row=1, col=1)
	fig.add_trace(go.Scatter(x=[label_x] * n, y=list(range(n)), mode='text', text=names, textposition='middle right',
							 textfont=dict(size=11), hovertext=names, hoverinfo='text', showlegend=False), row=1, col=1)

	if tip_annotation:
		tip_values = tip_annotation.get('values', {})
		tip_colors = _assign_colors([str(v) for v in tip_values.values()], tip_annotation.get('colors'))
		first = True
		for cat in sorted(tip_colors, key=_natural_key):
			xs, ys, texts = [], [], []
			for i, name in enumerate(names):
				if name in tip_values and str(tip_values[name]) == cat:
					xs.append(leaf_x[name])
					ys.append(i)
					texts.append('%s<br>%s: %s' % (name, tip_annotation.get('title', ''), cat))
			if len(xs) == 0:
				continue
			fig.add_trace(go.Scatter(x=xs, y=ys, mode='markers', marker=dict(size=8, color=tip_colors[cat],
									 line=dict(color='#333333', width=0.5)), name=cat, hovertext=texts,
									 hoverinfo='text', legendgroup='tips',
									 legendgrouptitle_text=tip_annotation.get('title', '') if first else None),
						  row=1, col=1)
			first = False

	fig.update_xaxes(row=1, col=1, range=[-0.02 * max_depth, label_x + label_px * data_per_px],
					 tickvals=_nice_ticks(max_depth), showgrid=False, zeroline=False, title_text='Branch length')

	# ---- categorical heatmap panels
	col_index = 2
	for hi, h in enumerate(heatmaps):
		columns = h['columns']
		values = h['values']
		cat_labels = h.get('category_labels', {})
		ordered_cats, colors = _heatmap_categories(h, names)
		code = {c: i for i, c in enumerate(ordered_cats)}
		k = max(1, len(ordered_cats))

		z, text = [], []
		for leaf in names:
			zrow, trow = [], []
			for col in columns:
				v = values.get(leaf, {}).get(col)
				if v is None:
					zrow.append(None)
					trow.append('%s<br>%s: NA' % (leaf, col))
				else:
					zrow.append(code[str(v)])
					trow.append('%s<br>%s: %s' % (leaf, col, cat_labels.get(str(v), str(v))))
			z.append(zrow)
			text.append(trow)

		colorscale = []
		for i, c in enumerate(ordered_cats):
			colorscale.append([i / float(k), colors[c]])
			colorscale.append([(i + 1) / float(k), colors[c]])
		if len(colorscale) == 0:
			colorscale = [[0, '#ffffff'], [1, '#ffffff']]

		small = n * len(columns) <= 40000
		fig.add_trace(go.Heatmap(z=z, x=columns, y=list(range(n)), text=text, hovertemplate='%{text}<extra></extra>',
								 colorscale=colorscale, zmin=-0.5, zmax=k - 0.5, showscale=False,
								 xgap=1 if small else 0, ygap=1 if small else 0), row=1, col=col_index)
		if h.get('show_legend', True):
			for ci, c in enumerate(ordered_cats):
				fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
										 marker=dict(size=11, symbol='square', color=colors[c],
													 line=dict(color='#999999', width=0.5)),
										 name=cat_labels.get(c, c), legendgroup='heatmap_%d' % hi,
										 legendgrouptitle_text=h.get('title', '') if ci == 0 else None,
										 hoverinfo='skip'), row=1, col=col_index)
		fig.update_xaxes(row=1, col=col_index, tickangle=-90, showgrid=False, zeroline=False,
						 showticklabels=len(columns) <= 200, title_text=h.get('column_label', ''))
		col_index += 1

	# ---- gene architecture panel
	if gene_panel:
		genes = gene_panel.get('genes', {})
		max_len = 1.0
		for leaf in genes:
			for g in genes[leaf]:
				max_len = max(max_len, float(g['end']))
		head = max_len * 0.012
		half = 0.32
		color_polys = defaultdict(lambda: [[], []])
		hx, hy, htext = [], [], []
		for i, leaf in enumerate(names):
			for g in genes.get(leaf, []):
				s, e = float(min(g['start'], g['end'])), float(max(g['start'], g['end']))
				hl = min(head, 0.5 * (e - s))
				if g.get('forward', True):
					px = [s, e - hl, e, e - hl, s, s]
				else:
					px = [e, s + hl, s, s + hl, e, e]
				py = [i - half, i - half, i, i + half, i + half, i - half]
				cp = color_polys[g.get('color', '#dbdbdb')]
				cp[0] += px + [None]
				cp[1] += py + [None]
				hx.append((s + e) / 2.0)
				hy.append(i)
				htext.append('%s<br>%s<br>Ortholog group: %s<br>%d-%d (%s)' % (leaf, g.get('name', ''), g.get('group', ''),
							 int(s), int(e), '+' if g.get('forward', True) else '-'))
		for color, (px, py) in color_polys.items():
			fig.add_trace(go.Scatter(x=px, y=py, mode='lines', fill='toself', fillcolor=color,
									 line=dict(color='#333333', width=0.5), hoverinfo='skip', showlegend=False),
						  row=1, col=col_index)
		fig.add_trace(go.Scatter(x=hx, y=hy, mode='markers', marker=dict(size=10, opacity=0), hovertext=htext,
								 hoverinfo='text', showlegend=False), row=1, col=col_index)
		fig.update_xaxes(row=1, col=col_index, showgrid=False, zeroline=False, title_text='Position (bp)')

	fig.update_yaxes(range=[n - 0.5, -0.5], showticklabels=False, showgrid=False, zeroline=False)
	fig.update_layout(title=dict(text='%s<br><sup>%s</sup>' % (title, HELP_TEXT), x=0.0, xanchor='left'),
					  width=fig_width, height=fig_height, margin=margin, plot_bgcolor='white', paper_bgcolor='white',
					  legend=dict(x=1.0, xanchor='left', y=1.0, yanchor='top', groupclick='toggleitem'),
					  hoverlabel=dict(bgcolor='white'))

	base = os.path.splitext(os.path.basename(output_html))[0]
	fig.write_html(output_html, include_plotlyjs='cdn', full_html=True,
				   config={'displaylogo': False, 'responsive': False,
						   'toImageButtonOptions': {'format': 'svg', 'filename': base}})


# ---------------------------------------------------------------------------------------------------------------------
# Static (PDF) rendering with matplotlib.
# ---------------------------------------------------------------------------------------------------------------------

def _mpl_text(s):
	"""Escape dollar signs so matplotlib does not interpret text from data (e.g. sample names) as math."""
	return str(s).replace('$', r'\$')


def _hex_to_rgb(color):
	from matplotlib.colors import to_rgb
	return to_rgb(color)


def _legend_entries(heatmaps, tip_annotation, names):
	"""Collect legend sections as [title, [(label, color, kind)]] where kind is 'square' or 'circle'."""
	sections = []
	if tip_annotation:
		tip_values = tip_annotation.get('values', {})
		tip_colors = _assign_colors([str(v) for v in tip_values.values()], tip_annotation.get('colors'))
		present = set(str(tip_values[x]) for x in names if x in tip_values)
		entries = [(c, tip_colors[c], 'circle') for c in sorted(tip_colors, key=_natural_key) if c in present]
		if entries:
			sections.append([tip_annotation.get('title', ''), entries])
	for h in heatmaps:
		if not h.get('show_legend', True):
			continue
		ordered_cats, colors = _heatmap_categories(h, names)
		cat_labels = h.get('category_labels', {})
		entries = [(cat_labels.get(c, c), colors[c], 'square') for c in ordered_cats]
		if entries:
			sections.append([h.get('title', ''), entries])
	return sections


def _draw_legend(ax, sections, fontsize):
	"""Draw legend sections into an (invisible) axes."""
	from matplotlib.lines import Line2D
	from matplotlib.patches import Patch
	ax.axis('off')
	handles, labels, title_indices = [], [], []
	for si, (title, entries) in enumerate(sections):
		if si > 0:
			handles.append(Patch(facecolor='none', edgecolor='none'))
			labels.append('')
		title_indices.append(len(handles))
		handles.append(Patch(facecolor='none', edgecolor='none'))
		labels.append(_mpl_text(title or ''))
		for label, color, kind in entries:
			if kind == 'circle':
				handles.append(Line2D([0], [0], marker='o', color='none', markerfacecolor=color,
									  markeredgecolor='#333333', markeredgewidth=0.5, markersize=fontsize * 0.9))
			else:
				handles.append(Patch(facecolor=color, edgecolor='#999999', linewidth=0.5))
			labels.append(_mpl_text(label))
	if not handles:
		return
	rows_per_col = 45
	ncol = max(1, int(math.ceil(len(handles) / float(rows_per_col))))
	legend = ax.legend(handles, labels, loc='upper left', frameon=False, fontsize=fontsize, ncol=ncol,
					   handlelength=1.2, handletextpad=0.5, borderaxespad=0.0, columnspacing=1.0)
	texts = legend.get_texts()
	for ti in title_indices:
		texts[ti].set_fontweight('bold')


def create_phylo_heatmap_pdf(tree_file, output_pdf, title, heatmaps=None, gene_panel=None, tip_annotation=None,
							 midpoint_root=True, show_tip_labels=True):
	"""
	Create a static PDF phylogenetic heatmap (matplotlib) - takes the same figure specification as
	create_phylo_heatmap_html. The page grows with the number of leaves and columns (up to the maximum page size most PDF
	viewers support, after which rows/columns and fonts are shrunk to fit).
	"""
	from matplotlib.figure import Figure
	from matplotlib.collections import LineCollection, PolyCollection
	import numpy as np

	heatmaps = [h for h in (heatmaps or []) if len(h.get('columns', [])) > 0]
	names, seg_x, seg_y, leaf_x, max_depth = tree_layout(tree_file, midpoint_root=midpoint_root)
	n = len(names)
	if n == 0:
		raise RuntimeError('No leaves found in phylogeny: %s' % tree_file)
	if max_depth <= 0:
		max_depth = 1.0

	legend_sections = _legend_entries(heatmaps, tip_annotation, names)

	# ---- sizing (inches)
	top_in, bottom_in = 0.9, 0.6
	max_col_label = 0
	for h in heatmaps:
		if len(h['columns']) <= 200:
			max_col_label = max([max_col_label] + [len(str(c)) for c in h['columns']])
	bottom_in += min(2.5, 0.075 * max_col_label) + (0.3 if heatmaps or gene_panel else 0.0)

	row_in = 0.2
	if n * row_in + top_in + bottom_in > MAX_PDF_INCHES:
		row_in = (MAX_PDF_INCHES - top_in - bottom_in) / float(n)
	fontsize = max(1.5, min(8.0, row_in * 72 * 0.7))

	max_label_len = max(len(x) for x in names)
	label_in = min(4.5, 0.075 * max_label_len * (fontsize / 8.0) + 0.2) if show_tip_labels else 0.0
	tree_in = 3.0
	panel_in = [tree_in + label_in]
	cell_in = 0.18
	total_columns = sum(len(h['columns']) for h in heatmaps)
	fixed_in = sum(panel_in) + (8.0 if gene_panel else 0.0) + 0.5 * (len(heatmaps) + (1 if gene_panel else 0))
	legend_in = 0.0
	if legend_sections:
		longest = max(len(str(label)) for _, entries in legend_sections for label, _, _ in entries)
		longest = max([longest] + [len(t) for t, _ in legend_sections])
		legend_rows = sum(len(entries) + 2 for _, entries in legend_sections)
		legend_in = min(6.0, 0.08 * longest + 0.6) * max(1, int(math.ceil(legend_rows / 45.0)))
	if total_columns > 0 and fixed_in + legend_in + total_columns * cell_in > MAX_PDF_INCHES:
		cell_in = max(0.005, (MAX_PDF_INCHES - fixed_in - legend_in) / float(total_columns))
	for h in heatmaps:
		panel_in.append(max(1.2, len(h['columns']) * cell_in))
	if gene_panel:
		panel_in.append(8.0)
	gap_in = 0.5
	left_in, right_in = 0.3, 0.3 + legend_in

	fig_width = left_in + sum(panel_in) + gap_in * (len(panel_in) - 1) + right_in
	fig_height = top_in + n * row_in + bottom_in
	fig_width = min(MAX_PDF_INCHES, max(6.0, fig_width))
	fig_height = min(MAX_PDF_INCHES, max(4.0, fig_height))

	fig = Figure(figsize=(fig_width, fig_height))
	plot_h = n * row_in / fig_height
	# anchor the panels just below the title (any extra space from the minimum page size goes to the bottom)
	plot_b = max(0.0, fig_height - top_in - n * row_in) / fig_height

	axes = []
	x_cursor = left_in
	for pi, w in enumerate(panel_in):
		ax = fig.add_axes([x_cursor / fig_width, plot_b, w / fig_width, plot_h], sharey=axes[0] if axes else None)
		axes.append(ax)
		x_cursor += w + gap_in
	for ax in axes:
		ax.set_ylim(n - 0.5, -0.5)
		ax.tick_params(axis='y', left=False, labelleft=False)
		for side in ['top', 'right', 'left']:
			ax.spines[side].set_visible(False)

	# ---- phylogeny panel
	ax = axes[0]
	segments, current = [], []
	for x, y in zip(seg_x, seg_y):
		if x is None:
			if len(current) == 2:
				segments.append(current)
			current = []
		else:
			current.append((x, y))
	ax.add_collection(LineCollection(segments, colors='#444444', linewidths=0.8))
	label_x = max_depth * 1.03
	tree_fraction = tree_in / float(tree_in + label_in)
	ax.set_xlim(-0.02 * max_depth, label_x / tree_fraction)
	if show_tip_labels:
		guides = [[(leaf_x[name], i), (label_x * 0.995, i)] for i, name in enumerate(names)
				  if leaf_x[name] < max_depth * 1.01]
		ax.add_collection(LineCollection(guides, colors='#c8c8c8', linewidths=0.4, linestyles='dotted'))
		for i, name in enumerate(names):
			ax.text(label_x, i, _mpl_text(name), va='center', ha='left', fontsize=fontsize, clip_on=False)
	if tip_annotation:
		tip_values = tip_annotation.get('values', {})
		tip_colors = _assign_colors([str(v) for v in tip_values.values()], tip_annotation.get('colors'))
		xs, ys, cs = [], [], []
		for i, name in enumerate(names):
			if name in tip_values:
				xs.append(leaf_x[name])
				ys.append(i)
				cs.append(tip_colors[str(tip_values[name])])
		ax.scatter(xs, ys, c=cs, s=(fontsize * 1.1) ** 2, edgecolors='#333333', linewidths=0.4, zorder=3)
	from matplotlib.ticker import MaxNLocator
	ax.set_xticks([t for t in MaxNLocator(nbins=4).tick_values(0, max_depth) if 0 <= t <= max_depth * 1.0001])
	ax.spines['bottom'].set_bounds(0, max_depth)
	ax.tick_params(axis='x', labelsize=max(fontsize, 5))
	ax.set_xlabel('Branch length', fontsize=max(fontsize, 6))
	ax.set_title('Phylogeny', fontsize=max(fontsize + 1, 7))

	# ---- heatmap panels
	pi = 1
	for h in heatmaps:
		ax = axes[pi]
		columns = h['columns']
		ordered_cats, colors = _heatmap_categories(h, names)
		rgb_colors = {c: _hex_to_rgb(colors[c]) for c in ordered_cats}
		img = np.ones((n, len(columns), 3))
		for i, leaf in enumerate(names):
			leaf_values = h['values'].get(leaf, {})
			for j, col in enumerate(columns):
				v = leaf_values.get(col)
				if v is not None:
					img[i, j, :] = rgb_colors[str(v)]
		ax.imshow(img, aspect='auto', interpolation='nearest', extent=(-0.5, len(columns) - 0.5, n - 0.5, -0.5))
		if n * len(columns) <= 40000 and cell_in >= 0.05 and row_in >= 0.05:
			grid = [[(j - 0.5, -0.5), (j - 0.5, n - 0.5)] for j in range(1, len(columns))]
			grid += [[(-0.5, i - 0.5), (len(columns) - 0.5, i - 0.5)] for i in range(1, n)]
			ax.add_collection(LineCollection(grid, colors='white', linewidths=0.6))
		ax.set_ylim(n - 0.5, -0.5)
		if len(columns) <= 200:
			ax.set_xticks(range(len(columns)))
			ax.set_xticklabels([_mpl_text(c) for c in columns], rotation=90, fontsize=max(1.5, min(fontsize, cell_in * 72 * 0.75)))
		else:
			ax.set_xticks([])
		ax.set_xlabel(_mpl_text(h.get('column_label', '')), fontsize=max(fontsize, 6))
		ax.set_title(_mpl_text(h.get('title', '')), fontsize=max(fontsize + 1, 7))
		ax.spines['bottom'].set_visible(False)
		pi += 1

	# ---- gene architecture panel
	if gene_panel:
		ax = axes[pi]
		genes = gene_panel.get('genes', {})
		max_len = 1.0
		for leaf in genes:
			for g in genes[leaf]:
				max_len = max(max_len, float(g['end']))
		head = max_len * 0.012
		half = 0.32
		polys, face_colors = [], []
		for i, leaf in enumerate(names):
			for g in genes.get(leaf, []):
				s, e = float(min(g['start'], g['end'])), float(max(g['start'], g['end']))
				hl = min(head, 0.5 * (e - s))
				if g.get('forward', True):
					px = [s, e - hl, e, e - hl, s]
				else:
					px = [e, s + hl, s, s + hl, e]
				py = [i - half, i - half, i, i + half, i + half]
				polys.append(list(zip(px, py)))
				face_colors.append(g.get('color', '#dbdbdb'))
		ax.add_collection(PolyCollection(polys, facecolors=face_colors, edgecolors='#333333', linewidths=0.3))
		ax.set_xlim(-0.01 * max_len, max_len * 1.01)
		ax.set_ylim(n - 0.5, -0.5)
		ax.tick_params(axis='x', labelsize=max(fontsize, 5))
		ax.set_xlabel('Position (bp)', fontsize=max(fontsize, 6))
		ax.set_title(gene_panel.get('title', 'Gene architecture'), fontsize=max(fontsize + 1, 7))

	if legend_sections:
		lax = fig.add_axes([(fig_width - right_in + 0.3) / fig_width, plot_b, max(0.1, legend_in) / fig_width,
							plot_h])
		_draw_legend(lax, legend_sections, fontsize=max(6, min(9, fontsize + 1)))

	fig.suptitle(_mpl_text(title), x=left_in / fig_width, ha='left', y=1.0 - 0.25 / fig_height, va='top', fontsize=12)
	with MPL_LOCK:
		fig.savefig(output_pdf, format='pdf')


def create_legend_pdf(category_colors, output_pdf, title=''):
	"""
	Create a PDF with just a color legend.
	- category_colors: ordered list of [category, color] (or a dictionary of category -> color).
	"""
	from matplotlib.figure import Figure
	if isinstance(category_colors, dict):
		category_colors = list(category_colors.items())
	entries = [(str(c), col, 'square') for c, col in category_colors]
	rows = len(entries) + 2
	ncol = max(1, int(math.ceil(rows / 45.0)))
	longest = max([len(title)] + [len(e[0]) for e in entries]) if entries else len(title)
	width = min(MAX_PDF_INCHES, max(3.0, ncol * (0.09 * longest + 0.8)))
	height = min(MAX_PDF_INCHES, max(2.0, 0.22 * min(rows, 45) + 0.6))
	fig = Figure(figsize=(width, height))
	ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
	_draw_legend(ax, [[title, entries]], fontsize=9)
	with MPL_LOCK:
		fig.savefig(output_pdf, format='pdf')


def safely_create(func, *args, logObject=None, **kwargs):
	"""Run an interactive figure creation function, logging (but not raising) any issues."""
	try:
		func(*args, **kwargs)
		if logObject:
			logObject.info('Created interactive visualization: %s' % str(args[1] if len(args) > 1 else ''))
		return True
	except Exception:
		msg = 'Warning: issues creating interactive (HTML) visualization - skipping it.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		if logObject:
			logObject.warning(msg)
			logObject.warning(traceback.format_exc())
		return False


# ---------------------------------------------------------------------------------------------------------------------
# Tool-specific figure specifications, built from the track files produced by the lsaBGC programs, and their HTML/PDF
# renderings.
# ---------------------------------------------------------------------------------------------------------------------

def _read_tsv(path, header=True):
	rows = []
	with open(path) as of:
		for i, line in enumerate(of):
			line = line.rstrip('\n')
			if not line:
				continue
			ls = line.split('\t')
			if header and i == 0:
				continue
			rows.append(ls)
	return rows


def _gseef_spec(track_file):
	values = defaultdict(dict)
	colors = {}
	col_order = {}
	for gcf, gcf_index, label, annotation, color in _read_tsv(track_file):
		col_order[gcf] = int(gcf_index)
		if annotation != 'absent':
			values[label][gcf] = annotation
			colors[annotation] = color.strip('"')
	columns = sorted(col_order, key=lambda x: col_order[x])
	return dict(title='GSeeF: GCF presence across the species phylogeny',
				heatmaps=[{'title': 'GCF presence (BGC class)', 'columns': columns, 'values': values,
						   'colors': colors, 'column_label': 'GCF IDs'}])


def gseef_html(tree_file, track_file, output_html):
	"""GSeeF: GCF presence (colored by BGC annotation class) across samples."""
	create_phylo_heatmap_html(tree_file, output_html, **_gseef_spec(track_file))


def gseef_pdf(tree_file, track_file, output_pdf, legend_pdf=None):
	"""GSeeF static version (and, optionally, a separate PDF with just the annotation class legend)."""
	spec = _gseef_spec(track_file)
	create_phylo_heatmap_pdf(tree_file, output_pdf, **spec)
	if legend_pdf:
		colors = spec['heatmaps'][0]['colors']
		create_legend_pdf([[c, colors[c]] for c in sorted(colors, key=_natural_key)], legend_pdf,
						  title='BGC annotation class')


def _see_spec(gggenes_track_file, heatmap_track_file, title):
	genes = defaultdict(list)
	og_colors = {}
	for label, gene, start, end, forward, og, og_color in _read_tsv(gggenes_track_file):
		og_color = og_color.strip('"')
		if start == 'NA' or end == 'NA':
			continue
		if og != 'singleton':
			og_colors[og] = og_color
		genes[label].append({'start': int(start), 'end': int(end), 'forward': forward.upper() == 'TRUE',
							 'name': gene, 'group': og, 'color': og_color})
	values = defaultdict(dict)
	og_counts = {}
	for label, og, og_presence, og_count in _read_tsv(heatmap_track_file):
		og_counts[og] = int(og_count)
		if og_presence != 'Absent':
			values[label][og] = og_presence
	columns = sorted(og_counts, key=lambda x: (-og_counts[x], _natural_key(x)))
	return dict(title=title,
				heatmaps=[{'title': 'Ortholog group presence', 'columns': columns, 'values': values,
						   'colors': og_colors, 'show_legend': False, 'column_label': 'Ortholog groups'}],
				gene_panel={'title': 'Gene architecture', 'genes': genes})


def see_html(tree_file, gggenes_track_file, heatmap_track_file, output_html, title):
	"""lsaBGC-See: ortholog group presence and gene architecture for BGCs of a GCF."""
	create_phylo_heatmap_html(tree_file, output_html, **_see_spec(gggenes_track_file, heatmap_track_file, title))


def see_pdf(tree_file, gggenes_track_file, heatmap_track_file, output_pdf, title):
	create_phylo_heatmap_pdf(tree_file, output_pdf, **_see_spec(gggenes_track_file, heatmap_track_file, title))


def _comprehenseeive_spec(heatmap_track_file, detection_track_file, title):
	values = defaultdict(dict)
	og_samples = defaultdict(int)
	for label, og, og_copy in _read_tsv(heatmap_track_file):
		values[label][og] = og_copy
		if og_copy != '0':
			og_samples[og] += 1
	columns = sorted(set(og for leaf in values for og in values[leaf]), key=lambda x: (-og_samples[x], _natural_key(x)))
	detection = {}
	if detection_track_file and os.path.isfile(detection_track_file):
		for name, method in _read_tsv(detection_track_file):
			detection[name] = method
	return dict(title=title,
				heatmaps=[{'title': 'Ortholog group copy count', 'columns': columns, 'values': values,
						   'colors': {'0': '#ffffff', '1': '#949292', 'Multi': '#403f3f'},
						   'category_labels': {'0': 'Absent', '1': 'Single-copy', 'Multi': 'Multi-copy'},
						   'column_label': 'Ortholog groups'}],
				tip_annotation={'title': 'GCF detection', 'values': detection} if detection else None)


def comprehenseeive_html(tree_file, heatmap_track_file, detection_track_file, output_html, title):
	"""lsaBGC-ComprehenSeeIve: copy counts of GCF ortholog groups across all samples."""
	create_phylo_heatmap_html(tree_file, output_html,
							  **_comprehenseeive_spec(heatmap_track_file, detection_track_file, title))


def comprehenseeive_pdf(tree_file, heatmap_track_file, detection_track_file, output_pdf, title):
	create_phylo_heatmap_pdf(tree_file, output_pdf,
							 **_comprehenseeive_spec(heatmap_track_file, detection_track_file, title))


def _sociate_spec(track_file, title):
	values = defaultdict(dict)
	col_order = defaultdict(list)
	for allele, label, order, value in _read_tsv(track_file, header=False):
		values[label][allele] = value
		try:
			col_order[allele].append(float(order))
		except ValueError:
			col_order[allele].append(0.0)
	columns = sorted(col_order, key=lambda x: (sum(col_order[x]) / len(col_order[x]), _natural_key(x)))
	return dict(title=title,
				heatmaps=[{'title': 'Associated features', 'columns': columns, 'values': values,
						   'colors': {'focal': '#242424', '-': '#6f7070', 'positive': '#4b7dc9',
									  'negative': '#db6588'},
						   'category_labels': {'focal': 'Focal GCF', '-': 'No direction',
											   'positive': 'Positively associated',
											   'negative': 'Negatively associated'},
						   'column_label': 'GCFs / ortholog groups'}])


def sociate_html(tree_file, track_file, output_html, title):
	"""lsaBGC-Sociate: presence of features (dis-)associated with a focal GCF."""
	create_phylo_heatmap_html(tree_file, output_html, **_sociate_spec(track_file, title))


def sociate_pdf(tree_file, track_file, output_pdf, title):
	create_phylo_heatmap_pdf(tree_file, output_pdf, **_sociate_spec(track_file, title))


def read_color_map(color_map_file):
	"""Read a two-column (category, color) mapping file with a header, as made by util.mapColorsToCategories."""
	colors = {}
	if color_map_file and os.path.isfile(color_map_file):
		with open(color_map_file) as ocmf:
			header = []
			for i, line in enumerate(ocmf):
				ls = line.rstrip('\n').split('\t')
				if i == 0:
					header = ls
					continue
				row = dict(zip(header, ls))
				if 'category' in row and 'color' in row:
					colors[row['category']] = row['color'].strip('"')
	return colors


def _reconcile_spec(track_file, color_map_file, title):
	values = defaultdict(dict)
	columns = []
	for label, col_type, value in _read_tsv(track_file):
		values[label][col_type] = value
		if col_type not in columns:
			columns.append(col_type)
	return dict(title=title,
				heatmaps=[{'title': 'Population & GCF', 'columns': columns, 'values': values,
						   'colors': read_color_map(color_map_file)}])


def reconcile_html(tree_file, track_file, color_map_file, output_html, title):
	"""lsaBGC-Reconcile: population and GCF membership of genes across an ortholog group gene tree."""
	create_phylo_heatmap_html(tree_file, output_html, **_reconcile_spec(track_file, color_map_file, title))


def reconcile_pdf(tree_file, track_file, color_map_file, output_pdf, title):
	spec = _reconcile_spec(track_file, color_map_file, title)
	# the populations/GCFs color legends are written as separate PDFs by lsaBGC-Reconcile
	spec['heatmaps'][0]['show_legend'] = False
	create_phylo_heatmap_pdf(tree_file, output_pdf, **spec)


def legend_pdf_from_color_map(color_map_file, output_pdf, title=''):
	"""Create a legend PDF from a two-column (category, color) mapping file, preserving its order."""
	entries = []
	with open(color_map_file) as ocmf:
		header = []
		for i, line in enumerate(ocmf):
			ls = line.rstrip('\n').split('\t')
			if i == 0:
				header = ls
				continue
			row = dict(zip(header, ls))
			if 'category' in row and 'color' in row:
				entries.append([row['category'], row['color'].strip('"')])
	create_legend_pdf(entries, output_pdf, title=title)


def population_designation_pdf(tree_file, sample_to_population, output_pdf, title):
	"""popstrat: species phylogeny with leaves colored by population designation (for one identity cutoff)."""
	create_phylo_heatmap_pdf(tree_file, output_pdf, title,
							 tip_annotation={'title': 'Population', 'values': sample_to_population})


def population_overview_html(tree_file, output_html, title, population_tracks, gcf_presence=None):
	"""
	Overview of population designations (one column per designation, e.g. per core AAI cutoff) and, optionally, GCF
	presence/absence across the species phylogeny.
	- population_tracks: ordered list of [column_name, {sample: population}].
	- gcf_presence: optional ordered list of [gcf_id, set_of_samples].
	"""
	heatmaps = []
	if population_tracks:
		pop_values = defaultdict(dict)
		columns = []
		for column, sample_pops in population_tracks:
			columns.append(column)
			for sample, pop in sample_pops.items():
				pop_values[sample][column] = pop
		heatmaps.append({'title': 'Populations', 'columns': columns, 'values': pop_values, 'show_legend': False,
						 'column_label': 'Population designation'})
	if gcf_presence:
		gcf_values = defaultdict(dict)
		columns = []
		for gcf, samples in gcf_presence:
			columns.append(gcf)
			for sample in samples:
				gcf_values[sample][gcf] = 'Present'
		heatmaps.append({'title': 'GCF presence', 'columns': columns, 'values': gcf_values,
						 'colors': {'Present': '#08519c'}, 'column_label': 'GCF IDs'})
	create_phylo_heatmap_html(tree_file, output_html, title, heatmaps=heatmaps)
