import os
import traceback
import statistics
import random
import subprocess
import concurrent.futures
from ete3 import Tree
from operator import itemgetter
from collections import defaultdict
from lsaBGC.classes.Pan import Pan
from lsaBGC import util
from pomegranate import *
import math
import warnings
from Bio import SeqIO

warnings.filterwarnings('ignore')

SEED = 12345

class GCF(Pan):
	def __init__(self, bgc_genbanks_listing, gcf_id='GCF_X', logObject=None, lineage_name='Unnamed lineage'):
		super().__init__(bgc_genbanks_listing, lineage_name=lineage_name, logObject=logObject)
		self.gcf_id = gcf_id

		#######
		## Variables not set during initialization
		#######

		# General variables
		self.hg_to_color = None
		self.hg_order_scores = defaultdict(lambda: ['NA', 'NA'])
		self.specific_core_homologs =set([])
		self.scc_homologs = set([])
		self.core_homologs = set([])
		self.protocluster_core_homologs = set([])

		# Sequence and alignment directories
		self.nucl_seq_dir = None
		self.prot_seq_dir = None
		self.prot_alg_dir = None
		self.codo_alg_dir = None
		self.nucl_filt_seq_dir = None
		self.prot_filt_seq_dir = None
		self.prot_filt_alg_dir = None
		self.codo_filt_alg_dir = None

		# Set of samples with sequencing reads to avoid for reporting alleles and novel SNVs,
		# these samples do not exhibit enough support for harboring a full BGC for the GC
		self.avoid_samples = set([])

	def identifyKeyHomologGroups(self, all_primary=False):
		"""
		Function that is not used currently in lsaBGC-Pan - in the original lsaBGC - this function 
		computed which orthogroups were core to a BGC, which were single-copy core and which were 
		commonly found in BGC protocore regions. This is a little more tricky to integrate with zol.
		"""
		try:
			initial_samples_with_at_least_one_gcf_hg = set([])
			for hg in self.hg_genes:
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					if not gene_info['is_expansion_bgc'] or all_primary:
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						initial_samples_with_at_least_one_gcf_hg.add(sample_id)

			for hg in self.hg_genes:
				sample_counts = defaultdict(int)
				sample_with_hg_as_protocluster_core = 0
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					if not gene_info['is_expansion_bgc'] or all_primary:
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						sample_counts[sample_id] += 1
						if gene_info['core_overlap']:
							sample_with_hg_as_protocluster_core += 1


				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				samples_with_any_copy = set([s[0] for s in sample_counts.items() if s[1] > 0])

				# check that hg is single-copy-core or just core
				if len(samples_with_single_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.scc_homologs.add(hg)
				if len(samples_with_any_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.core_homologs.add(hg)

				if len(samples_with_any_copy) > 0 and float(sample_with_hg_as_protocluster_core)/len(samples_with_any_copy) >= 0.5:
					self.protocluster_core_homologs.add(hg)

			self.logObject.info("Conserved Core Set of Homolog Groups:\t" + '; '.join(self.core_homologs))
			self.logObject.info("Conserved Single-Copy Set of Homolog Groups:\t" + '; '.join(self.scc_homologs))
			self.logObject.info("Proto-Core / Rule Based Homolog Groups:\t" + '; '.join(self.protocluster_core_homologs))
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with identifying key homolog groups.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def aggregateProteins(self, gcf_prots_file, draft_mode=False):
		"""
		Function to aggregate protein sequences from BGC genbanks and output to a file in FASTA format.
		"""
		try:
			gpf_handle = open(gcf_prots_file, 'w')
			for bgc in self.bgc_genes:
				sample = self.bgc_sample[bgc]
				set_id = bgc
				if draft_mode:
					set_id = sample
				for lt in self.bgc_genes[bgc]:
					lt_prot_seq = self.comp_gene_info[lt]['prot_seq']
					gpf_handle.write('>' + set_id + '|' + lt + '\n' + lt_prot_seq + '\n')
			gpf_handle.close()
			assert(os.path.isfile(gcf_prots_file))
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with aggregating proteins from BGCs belonging to GCF to a FASTA file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def modifyPhylogenyForSamplesWithMultipleBGCs(self, input_phylogeny, result_phylogeny, prune_set=None):
		"""
		Definition:
		Function which takes in an input phylogeny and produces a replicate resulting phylogeny with samples/leafs which
		have multiple BGC instances for a GCF expanded.
		********************************************************************************************************************
		Parameters:
		- input_phylogeny: The input newick phylogeny file.
		- result_phylogeny: The resulting newick phylogeny file to create.
		********************************************************************************************************************
		"""
		try:
			number_of_added_leaves = 0
			t = Tree(input_phylogeny)
			if prune_set != None:
				t.prune(prune_set)

			for node in t.traverse('postorder'):
				if node.name in self.sample_bgcs and len(self.sample_bgcs[node.name]) > 1:
					og_node_name = node.name
					node.name = node.name + '_INNERNODE'
					for bgc_id in self.sample_bgcs[og_node_name]:
						# if bgc_id == node.name: continue
						node.add_child(name=bgc_id)
						child_node = t.search_nodes(name=bgc_id)[0]
						child_node.dist = 0
						if bgc_id != og_node_name: number_of_added_leaves += 1

			t.write(format=0, outfile=result_phylogeny)
			if self.logObject:
				self.logObject.info(
					"New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (
						number_of_added_leaves, result_phylogeny))
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def assignColorsToHGs(self, gene_to_hg, bgc_genes, outdir):
		"""
		Description:
		Simple function to associate each homolog group with a color for consistent coloring.
		********************************************************************************************************************
		Parameters:
		- gene_to_hg: gene to HG relationship.
		- bgc_genes:  set of genes per HG.
		********************************************************************************************************************
		"""

		hg_bgc_counts = defaultdict(int)
		for b in bgc_genes:
			for g in bgc_genes[b]:
				if g in gene_to_hg:
					hg_bgc_counts[gene_to_hg[g]] += 1

		hgs = set([])
		for c in hg_bgc_counts:
			if hg_bgc_counts[c] > 1:
				hgs.add(c)

		len_hgs = len(hgs)
		color_listing_file = outdir + 'colors_for_hgs.txt'
		util.generateColors(outdir, color_listing_file, len_hgs)

		# read in list of colors
		colors = []
		with open(color_listing_file) as ocf:
			colors = [x.strip() for x in ocf.readlines()]
		random.Random(SEED).shuffle(colors)

		hg_to_color = {}
		for i, c in enumerate(set(hgs)):
			hg_to_color[c] = colors[i]
		self.hg_to_color = hg_to_color

	def createItolBGCSeeTrack(self, result_track_file):
		"""
		Description:
		Function to create a track file for visualizing BGC gene architecture across a phylogeny in the interactive tree
		of life (iTol)
		********************************************************************************************************************
		Parameters:
		- self: GCF object
		- result_track_file: The path to the resulting iTol track file for BGC gene visualization.
		********************************************************************************************************************
		"""
		try:
			track_handle = open(result_track_file, 'w')

			if self.logObject:
				self.logObject.info("Writing iTol track file to: %s" % result_track_file)
				self.logObject.info("Track will have label: %s" % self.gcf_id)

			# write header for iTol track file
			track_handle.write('DATASET_DOMAINS\n')
			track_handle.write('SEPARATOR TAB\n')
			track_handle.write('DATASET_LABEL\t%s\n' % self.gcf_id)
			track_handle.write('COLOR\t#000000\n')
			track_handle.write('BORDER_WIDTH\t1\n')
			track_handle.write('BORDER_COLOR\t#000000\n')
			track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
			track_handle.write('DATA\n')

			# write the rest of the iTol track file for illustrating genes across BGC instances
			ref_hg_directions = {}
			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = [bgc, str(last_gene_end)]
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]
					shape = 'None'
					if ginfo['direction'] == '+':
						shape = 'TR'
					elif ginfo['direction'] == '-':
						shape = 'TL'
					gstart = ginfo['start']
					gend = ginfo['end']
					hg_color = "#dbdbdb"
					if hg in self.hg_to_color:
						hg_color = self.hg_to_color[hg]
					gene_string = '|'.join([str(x) for x in [shape, gstart, gend, hg_color, hg]])
					printlist.append(gene_string)
					if hg != 'singleton':
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					track_handle.write('\t'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# flip the genbank visual if necessary, first BGC processed is used as reference guide
					if flip_support > keep_support:
						flip_printlist = printlist[:2]
						for gene_string in printlist[2:]:
							gene_info = gene_string.split('|')
							new_shape = None
							if gene_info[0] == 'TR':
								new_shape = 'TL'
							elif gene_info[0] == 'TL':
								new_shape = 'TR'
							new_gstart = int(last_gene_end) - int(gene_info[2])
							new_gend = int(last_gene_end) - int(gene_info[1])
							new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
							flip_printlist.append(new_gene_info)
						track_handle.write('\t'.join(flip_printlist) + '\n')
					else:
						track_handle.write('\t'.join(printlist) + '\n')
			track_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had difficulties creating iTol track for visualization of BGC gene architecture.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def visualizeComprehenSeeIveGCFViaR(self, orthofinder_matrix_file, phylogeny_newick_file, heatmap_track_file,
										detection_track_file, result_pdf_file, plot_rscript):
		"""
		Definition:
		Function to create tracks for visualization of gene architecture of BGCs belonging to GCF and run Rscript bgSee.R
		to produce automatic PDFs of plots. In addition, bgSee.R also produces a heatmap to more easily identify homolog
		groups which are conserved across isolates found to feature GCF.
		********************************************************************************************************************
		Parameters:
		- self: GCF object
		- gggenes_track_file: Path to file with gggenes track information (will be created/written to by function, if it doesn't exist!)
		- heatmap_track_file: Path to file for heatmap visual component (will be created/written to by function, if it doesn't exist!)
		- phylogeny_file: Phylogeny to use for visualization.
		- result_pdf_file: Path to PDF file where plots from bgSee.R will be written to.
		********************************************************************************************************************
		"""
		try:
			if os.path.isfile(heatmap_track_file) or os.path.isfile(detection_track_file):
				os.system('rm -f %s %s' % (detection_track_file, heatmap_track_file))
			heatmap_track_handle = open(heatmap_track_file, 'w')
			detection_track_handle = open(detection_track_file, 'w')

			if self.logObject:
				self.logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
				self.logObject.info("Writing detection-method input file to: %s" % detection_track_file)

			# write header for track files
			heatmap_track_handle.write('label\tog\tog_copy\n')
			detection_track_handle.write('name\tdetection_method\n')

			gcf_relevant_hgs = set([])
			for bgc in self.bgc_hgs:
				for hg in self.bgc_hgs[bgc]:
					gcf_relevant_hgs.add(hg)

			sample_detection_methods = defaultdict(set)
			for bgc in self.bgc_gbk:
				gbk_path = self.bgc_gbk[bgc]
				detection_method = 'antiSMASH/GECCO'
				bgc_sample = self.bgc_sample[bgc]
				sample_detection_methods[bgc_sample].add(detection_method)

			for sample in sample_detection_methods:
				dm = 'Not Detected'
				if len(sample_detection_methods[sample]) == 1:
					dm = 'antiSMASH/GECCO'
				detection_track_handle.write(sample + '\t' + dm + '\n')

			detection_track_handle.close()

			sample_names = []
			sample_hg_counts = defaultdict(lambda: defaultdict(int))
			with open(orthofinder_matrix_file) as omf:
				for i, line in enumerate(omf):
					line = line.strip('\n')
					ls = line.split('\t')
					if i == 0:
						sample_names = ls[1:]
					else:
						hg = ls[0]
						if not hg in gcf_relevant_hgs: continue
						for j, lts in enumerate(ls[1:]):
							lts = lts.strip()
							sname = sample_names[j]
							if lts == '':
								sample_hg_counts[sname][hg] = 0
							else:
								num_lts = len(lts.split(','))
								sample_hg_counts[sname][hg] = num_lts

			tree_obj = Tree(phylogeny_newick_file)
			for node in tree_obj.traverse('postorder'):
				if not node.is_leaf(): continue
				sname = node.name
				for hg in gcf_relevant_hgs:
					copy_count = '0'
					if sample_hg_counts[sname][hg] == 1:
						copy_count = '1'
					elif sample_hg_counts[sname][hg] > 1:
						copy_count = 'Multi'
					heatmap_track_handle.write(sname + '\t' + hg + '\t' + copy_count + '\n')
			heatmap_track_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error("Had difficulties creating tracks for visualization of BGC gene architecture along phylogeny using R libraries.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())
		
		pr = open(plot_rscript, 'w')
		pr.write('library(ggplot2)\n')
		pr.write('library(ggtree)\n')
		pr.write('library(gggenes)\n')
		pr.write('library(ape)\n')
		pr.write('library(phytools)\n')
		pr.write('library(aplot)\n')
		pr.write('library(dplyr)\n\n')

		pr.write('phylo.tree_file <- "' + phylogeny_newick_file + '"\n') 
		pr.write('heatmap.data_file <- "' + heatmap_track_file + '"\n')
		pr.write('track_file <- "' + detection_track_file + '"\n')
		pr.write('pdf_file <- "' + result_pdf_file + '"\n\n')

		pr.write('phylo.tree <- read.tree(phylo.tree_file)\n')
		pr.write('phylo.tree <- midpoint.root(phylo.tree)\n')
		pr.write('heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n')
		pr.write('track.data <- read.table(track_file, header=T, sep="\\t")\n')
		pr.write('tree.labels <- phylo.tree$tip.label\n')
		pr.write('og_colors <- c("#FFFFFF", "#949292", "#403f3f")\n')
		pr.write('names(og_colors) <- c("0", "1", "Multi")\n\n')

		pr.write('pdf(pdf_file, height=30, width=30)\n')
		pr.write('gg_tr <- ggtree(phylo.tree)\n')
		pr.write('gg_tr <- gg_tr %<+% track.data + geom_tippoint(aes(color=detection_method), show.legend=T, size=3)\n')
		pr.write('gg_hm <- ggplot(heatmap.data, aes(x = og, y = label, fill=as.factor(og_copy))) +\n')
		pr.write('theme_classic() + scale_fill_manual("Copy Count", values=og_colors) +\n')
		pr.write('xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color="white") +\n')
		pr.write('theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))\n')
		pr.write('gg_hm %>% insert_left(gg_tr, width=0.3)\n')
		pr.write('dev.off()\n\n')

		pr.close()

		rscript_plot_cmd = ["Rscript", plot_rscript]
		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
		try:
			subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert(os.path.isfile(result_pdf_file))
			self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

		if self.logObject:
			self.logObject.info('Plotting completed (I think successfully)!')

	def visualizeGCFViaR(self, gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file, plot_rscript):
		"""
		Definition:
		Function to create tracks for visualization of gene architecture of BGCs belonging to GCF and run Rscript bgSee.R
		to produce automatic PDFs of plots. In addition, bgSee.R also produces a heatmap to more easily identify homolog
		groups which are conserved across isolates found to feature GCF.
		********************************************************************************************************************
		Parameters:
		- self: GCF object.
		- gggenes_track_file: Path to file with gggenes track information (will be created/written to by function, if it 
		                      doesn't exist!)
		- heatmap_track_file: Path to file for heatmap visual component (will be created/written to by function, if it 
		                      doesn't exist!)-
		- phylogeny_file: Phylogeny to use for visualization.
		- result_pdf_file: Path to PDF file where plots from bgSee.R will be written to.
		********************************************************************************************************************
		"""
		try:
			if os.path.isfile(gggenes_track_file) or os.path.isfile(heatmap_track_file):
				os.system('rm -f %s %s' % (gggenes_track_file, heatmap_track_file))
			gggenes_track_handle = open(gggenes_track_file, 'w')
			heatmap_track_handle = open(heatmap_track_file, 'w')
			if self.logObject:
				self.logObject.info("Writing gggenes input file to: %s" % gggenes_track_file)
				self.logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
			# write header for track files
			gggenes_track_handle.write('label\tgene\tstart\tend\tforward\tog\tog_color\n')
			heatmap_track_handle.write('label\tog\tog_presence\tog_count\n')

			ref_hg_directions = {}

			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			tree_obj = Tree(phylogeny_file)
			bgc_weights = defaultdict(int)
			all_bgcs_in_tree = set([])
			for leaf in tree_obj:
				all_bgcs_in_tree.add(str(leaf).strip('\n').lstrip('-'))
				bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

			bgc_hg_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
			hg_counts = defaultdict(int)
			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				if not bgc in all_bgcs_in_tree: continue
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = []
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]

					gstart = ginfo['start']
					gend = ginfo['end']
					forward = "FALSE"
					if ginfo['direction'] == '+': forward = "TRUE"

					hg_color = '"#dbdbdb"'
					if hg in self.hg_to_color:
						hg_color = '"' + self.hg_to_color[hg] + '"'

					gene_string = '\t'.join([str(x) for x in [bgc, lt, gstart, gend, forward, hg, hg_color]])
					printlist.append(gene_string)
					if hg != 'singleton':
						bgc_hg_presence[bgc][hg] = hg
						hg_counts[hg] += bgc_weights[bgc]
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					gggenes_track_handle.write('\n'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# flip the genbank visual if necessary, first BGC processed is used as reference guide
					if flip_support > keep_support:
						flip_printlist = []
						for gene_string in printlist:
							gene_info = gene_string.split('\t')
							new_forward = 'TRUE'
							if gene_info[4] == 'TRUE': new_forward = 'FALSE'
							new_gstart = int(last_gene_end) - int(gene_info[3])
							new_gend = int(last_gene_end) - int(gene_info[2])
							new_gene_string = '\t'.join([str(x) for x in
														 [gene_info[0], gene_info[1], new_gstart, new_gend, new_forward,
															gene_info[-2], gene_info[-1]]])
							flip_printlist.append(new_gene_string)
						gggenes_track_handle.write('\n'.join(flip_printlist) + '\n')
					else:
						gggenes_track_handle.write('\n'.join(printlist) + '\n')

			dummy_hg = None
			for bgc in bgc_hg_presence:
				for hg in hg_counts:
					dummy_hg = hg
					heatmap_track_handle.write('\t'.join([bgc, hg, bgc_hg_presence[bgc][hg], str(hg_counts[hg])]) + '\n')

			for i, bgc in enumerate(all_bgcs_in_tree):
				if not bgc in bgc_gene_counts.keys():
					gggenes_track_handle.write('\t'.join([bgc] + ['NA']*4 + ['Absent', '"#FFFFFF"']) + '\n')
					heatmap_track_handle.write('\t'.join([bgc, dummy_hg, 'Absent', '1']) + '\n')
				elif i == 0:
					gggenes_track_handle.write('\t'.join([bgc] + ['NA']*4 + ['Absent', '"#FFFFFF"']) + '\n')

			gggenes_track_handle.close()
			heatmap_track_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties creating tracks for visualization of BGC gene architecture along phylogeny using R libraries.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		
		pr = open(plot_rscript, 'w')

		pr.write('library(ggplot2)\n')
		pr.write('library(ggtree)\n')
		pr.write('library(gggenes)\n')
		pr.write('library(ape)\n')
		pr.write('library(phytools)\n')
		pr.write('library(aplot)\n')
		pr.write('library(dplyr)\n\n')

		pr.write('phylo.tree_file <- "' + phylogeny_file + '"\n')
		pr.write('genes.data_file <- "' + gggenes_track_file + '"\n')
		pr.write('heatmap.data_file <- "' + heatmap_track_file + '"\n')
		pr.write('pdf_file <- "' + result_pdf_file + '"\n\n')

		pr.write('phylo.tree <- read.tree(phylo.tree_file)\n')
		pr.write('phylo.tree <- midpoint.root(phylo.tree)\n')
		pr.write('genes.data <- read.table(genes.data_file, header=T, sep="\\t")\n')
		pr.write('heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n')
		pr.write('tree.labels <- phylo.tree$tip.label\n')
		pr.write('genes.sub.data <- distinct(genes.data[c("og", "og_color")])\n')
		pr.write('og_colors <- na.omit(c(genes.sub.data$og_color))\n')
		pr.write('names(og_colors) <- c(na.omit(c(genes.sub.data$og)))\n\n')

		pr.write('pdf(pdf_file, height=30, width=30)\n')
		pr.write('gg_tr <- ggtree(phylo.tree)#  + ggplot2::xlim(NA, 1)\n')
		pr.write('gg_gn <- ggplot(genes.data, aes(xmin = start, xmax = end, y = label, fill = og, forward = forward)) +\n')
		pr.write('geom_gene_arrow(show.legend=F) + theme_void() + scale_fill_manual(values=og_colors)\n')
		pr.write('gg_hm <- ggplot(heatmap.data, aes(x = reorder(og, -og_count), y = label, fill=og_presence)) +\n')
		pr.write('theme_classic() + scale_fill_manual(values=og_colors) +\n')
		pr.write('xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color="white", show.legend=F) +\n')
		pr.write('theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))\n')
		pr.write('gg_hm %>% insert_left(gg_tr, width=0.4) %>% insert_right(gg_gn, width=1.0)\n')
		pr.write('dev.off()\n\n')
		pr.close()
		
		rscript_plot_cmd = ["Rscript", plot_rscript]

		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
		try:
			subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert(os.path.isfile(result_pdf_file))
			self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

		if self.logObject:
			self.logObject.info('Plotting completed (I think successfully)!')

	def constructCodonAlignments(self, outdir, threads=1, only_scc=False, list_alignments=False, filter_outliers=False, use_ms5=True):
		"""
		Definition:
		Function to automate construction of codon alignments. This function first extracts protein and nucleotide sequnces
		from BGC Genbanks, then creates protein alignments for each homolog group using MAFFT, and finally converts those
		into codon alignments using PAL2NAL.
		********************************************************************************************************************
		Parameters:
		- self: GCF object.
		- outdir: Path to output/workspace directory. Intermediate files (like extracted nucleotide and protein
				  sequences, protein and codon alignments, will be writen to respective subdirectories underneath this
				  one).
		- threads: Number of threads/threads to use when fake-parallelizing jobs using multiprocessing.
		- only_scc: Whether to construct codon alignments only for homolog groups which are found to be core and in
					single copy for samples with the GCF. Note, if working with draft genomes and the BGC is fragmented
					this should be able to still identify SCC homolog groups across the BGC instances belonging to the
					GCF.
		********************************************************************************************************************
		"""

		nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences/'
		prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences/'
		prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments/'
		codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments/'
		if filter_outliers:
			nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences_MAD_Refined/'
			prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences_MAD_Refined/'
			prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments_MAD_Refined/'
			codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments_MAD_Refined/'

		if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
		if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
		if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
		if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

		parallel_jobs_4thread = max(math.floor(threads / 4), 1)
		multi_thread = 4
		if threads < 4:
			multi_thread = threads
			parallel_jobs_4thread = 1

		all_samples = set(self.bgc_sample.values())
		try:
			inputs = []
			for hg in self.hg_genes:
				# if len(self.hg_genes[hg]) < 2: continue
				sample_counts = defaultdict(int)
				gene_sequences = {}
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					bgc_id = gene_info['bgc_name']
					sample_id = self.bgc_sample[bgc_id]
					nucl_seq = gene_info['nucl_seq']
					prot_seq = gene_info['prot_seq']
					sample_counts[sample_id] += 1
					gid = sample_id + '|' + gene
					if only_scc:
						gid = sample_id
					gene_sequences[gid] = tuple([nucl_seq, prot_seq])
				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				# check that hg is single-copy-core
				if only_scc and len(samples_with_single_copy.symmetric_difference(all_samples)) > 0:
					continue
				elif only_scc and self.logObject:
					self.logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % hg)
				# check that hg is present in the original instances of GCF
				#if len([x for x in gene_sequences.keys() if len(x.split('|')[1].split('_')[0]) == 3]) == 0: continue
				if filter_outliers:
					gene_sequences = util.determineOutliersByGeneLength(gene_sequences, self.logObject)
				inputs.append([hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, multi_thread, use_ms5, self.logObject])

			with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_jobs_4thread) as executor:
				executor.map(create_codon_msas, inputs)
			
			if not filter_outliers:
				self.nucl_seq_dir = nucl_seq_dir
				self.prot_seq_dir = prot_seq_dir
				self.prot_alg_dir = prot_alg_dir
				self.codo_alg_dir = codo_alg_dir
			else:
				self.nucl_filt_seq_dir = nucl_seq_dir
				self.prot_filt_seq_dir = prot_seq_dir
				self.prot_filt_alg_dir = prot_alg_dir
				self.codo_filt_alg_dir = codo_alg_dir

			if list_alignments:
				codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.txt'
				if filter_outliers:
					codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.MAD_Refined.txt'
				codon_alg_listings_handle = open(codon_alg_listings_file, 'w')
				for f in os.listdir(codo_alg_dir):
					codon_alg_listings_handle.write(f.split('.msa.fna')[0] + '\t' + codo_alg_dir + f + '\n')
				codon_alg_listings_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def constructGCFPhylogeny(self, output_alignment, output_phylogeny, only_scc=False, ambiguious_position_cutoff=0.000001):
		"""
		Definition:
		Function to create phylogeny based on codon alignments of SCC homolog groups for GCF.
		********************************************************************************************************************
		Parameters:
		- self: GCF object.
		- output_alignment: Path to output file for concatenated SCC homolog group alignment.
		- output_phylogeny: Path to output file for approximate maximum-likelihood phylogeny produced by FastTree2 from
							concatenated SCC homolog group alignment.
		********************************************************************************************************************
		"""
		try:
			if only_scc:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					# concatenate gene alignments
					with open(hg_align_msa) as opm:
						for rec in SeqIO.parse(opm, 'fasta'):
							bgc_sccs['>' + rec.id] += str(rec.seq).upper()

				for b in bgc_sccs:
					fasta_data.append([b] + list(bgc_sccs[b]))

				for i, ls in enumerate(zip(*fasta_data)):
					if i == 0:
						fasta_data_tr.append(ls)
					else:
						n_count = len([x for x in ls if x == '-'])
						if (float(n_count) / len(ls)) < 0.1:
							fasta_data_tr.append(list(ls))

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()
			else:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					#print(f)
					# perform consensus calling
					sample_seqs = defaultdict(list)
					with open(hg_align_msa) as opm:
						for rec in SeqIO.parse(opm, 'fasta'):
							sample = rec.id.split('|')[0]
							sample_seqs[sample].append(list(str(rec.seq).upper()))

					for samp in sample_seqs:
						samp_seqs = sample_seqs[samp]
						consensus_seq = []
						for alleles in zip(*samp_seqs):
							valid_alleles = set([a for a in list(alleles) if a in set(['A', 'C', 'G', 'T'])])
							if len(valid_alleles) == 1:
								consensus_seq.append(list(valid_alleles)[0])
							else:
								consensus_seq.append('-')
						bgc_sccs['>' + samp] += "".join(consensus_seq)

				for b in bgc_sccs:
					fasta_data.append([b] + list(bgc_sccs[b]))

				for i, ls in enumerate(zip(*fasta_data)):
					if i == 0:
						fasta_data_tr.append(ls)
					else:
						n_count = len([x for x in ls if x == '-'])
						if (float(n_count) / len(ls)) < ambiguious_position_cutoff:
							fasta_data_tr.append(list(ls))

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error('Had issues with creating concatenated alignment of the SCC homolog groups.')
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		# use FastTree2 to construct phylogeny
		fasttree_cmd = ['fasttree', '-nt', output_alignment, '>', output_phylogeny]
		if self.logObject:
			self.logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
		try:
			subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(fasttree_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(fasttree_cmd))

def create_codon_msas(inputs):
	"""
	Definition:
	Helper function which is to be called from the constructCodonAlignments() function to parallelize construction
	of codon alignments for each homolog group of interest in the GCF.
	********************************************************************************************************************
	Parameters:
	- inputs: list of inputs passed in by GCF.constructCodonAlignments().
	********************************************************************************************************************
	"""
	hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, threads, use_ms5, logObject = inputs

	hg_nucl_fasta = nucl_seq_dir + '/' + hg + '.fna'
	hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
	hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
	hg_codo_msa = codo_alg_dir + '/' + hg + '.msa.fna'

	hg_nucl_handle = open(hg_nucl_fasta, 'w')
	hg_prot_handle = open(hg_prot_fasta, 'w')

	for s in sorted(gene_sequences):
		hg_nucl_handle.write('>' + s + '\n' + str(gene_sequences[s][0]) + '\n')
		hg_prot_handle.write('>' + s + '\n' + str(gene_sequences[s][1]) + '\n')
	hg_nucl_handle.close()
	hg_prot_handle.close()

	align_cmd = ['muscle', '-super5', hg_prot_fasta, '-output', hg_prot_msa, '-amino', '-threads', str(threads)]
	pal2nal_cmd = ['pal2nal.pl', hg_prot_msa, hg_nucl_fasta, '-output', 'fasta', '>', hg_codo_msa]

	if logObject:
		logObject.info('Running multiple sequence alignment with the following command: %s' % ' '.join(align_cmd))
	try:
		subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(align_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(align_cmd))

	if logObject:
		logObject.info('Running PAL2NAL with the following command: %s' % ' '.join(pal2nal_cmd))
	try:
		subprocess.call(' '.join(pal2nal_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(pal2nal_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(pal2nal_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(pal2nal_cmd))

	if logObject:
		logObject.info('Achieved codon alignment for homolog group %s' % hg)