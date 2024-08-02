import os
import sys
import logging
import traceback
import subprocess
from collections import defaultdict
from Bio import SeqIO
from lsaBGC.classes.BGC import BGC
from lsaBGC import util
import statistics
import random
import multiprocessing
import math
from operator import itemgetter
from scipy import stats
import decimal
import _pickle as cPickle

class Pan:
	def __init__(self, bgc_genbanks_listing, logObject=None, lineage_name='Unnamed lineage'):
		self.bgc_genbanks_listing = bgc_genbanks_listing
		self.lineage_name = lineage_name
		self.logObject = logObject

		#######
		## Variables not set during initialization
		#######

		# General variables
		self.pan_bgcs = {}
		self.comp_gene_info = {}
		self.bgc_info = {}
		self.bgc_gbk = {}
		self.bgc_genes = {}
		self.bgc_edgy = {}
		self.pan_genes = set([])
		self.bgc_sample = {}
		self.sample_bgcs = defaultdict(set)
		self.bgc_product = {}
		self.bgc_core_counts = {}
		self.bgc_population = None
		self.sample_population = None

		# homology related variables
		self.gene_to_hg = None
		self.hg_genes = None
		self.hg_median_copy_count = None
		self.hg_prop_multi_copy = None
		self.bgc_hgs = defaultdict(set)

		### other variables
		# Jaccard similarity distance between bgcs
		self.pairwise_relations = None
		# Absolute pearson correlation coefficient for syntenic similarity between bgcs
		self.pairwise_syntenic_relations = None
		# Containment relations between bgcs
		self.pairwise_containment_relations = None

		# variables containing location of files
		self.final_stats_file = None
		self.final_stats_expanded_singletons_file = None
		self.pair_relations_txt_file = None
		self.bgc_to_gcf_map_file = None

	def readInBGCGenbanks(self, comprehensive_parsing=True, prune_set=None, edge_dist_cutoff=5000):
		"""
		Function to parse file listing location of BGC Genbanks.

		:param comprehensive_parsing (optional): flag specifying whether to perform comprehensive extraction of information from Genbanks. default is True.
		"""
		sample_index = defaultdict(int)
		with open(self.bgc_genbanks_listing) as obsf:
			for i, line in enumerate(obsf):
				line = line.strip()
				try:
					assert (len(line.split('\t')) == 3)
				except Exception as e:
					msg = "More or less than three columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1)
					if self.logObject:
						self.logObject.error(msg)
						self.logObject.error(traceback.format_exc())
					sys.stderr.write(msg + '\n')
					sys.exit(1)

				sample, gbk, prediction_method, dist_to_edge = line.split('\t')

				if prune_set != None and not sample in prune_set: continue
				sample = util.cleanUpSampleName(sample)
				try:
					dist_to_edge = float(dist_to_edge)
					edgy_bgc = False
					if dist_to_edge <= edge_dist_cutoff:
						edgy_bgc = True
					
					if prune_set != None and not sample in prune_set: continue
					assert (util.is_genbank(gbk))
					bgc_id = sample
					if sample_index[sample] > 0:
						bgc_id = sample + '_' + str(sample_index[sample] + 1)
					sample_index[sample] += 1

					is_expansion_bgc = False

					# Parse genbank using

					BGC_Object = BGC(gbk, bgc_id, is_expansion_bgc=is_expansion_bgc, prediction_method=prediction_method)
					BGC_Object.parseGenbanks(comprehensive_parsing=comprehensive_parsing)

					self.pan_bgcs[bgc_id] = BGC_Object
					self.comp_gene_info.update(BGC_Object.gene_information)
					self.bgc_info[bgc_id] = BGC_Object.cluster_information
					self.bgc_genes[bgc_id] = set(BGC_Object.gene_information)
					self.pan_genes = self.pan_genes.union(BGC_Object.gene_information.keys())
					self.bgc_gbk[bgc_id] = gbk
					self.bgc_sample[bgc_id] = sample
					self.sample_bgcs[sample].add(bgc_id)
					self.bgc_product[bgc_id] = [x['product'] for x in BGC_Object.cluster_information]
					self.bgc_core_counts[bgc_id] = BGC_Object.cluster_information[0]['count_core_gene_groups']
					self.bgc_edgy[bgc_id] = edgy_bgc

					if self.logObject:
						self.logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
				except Exception as e:
					if self.logObject:
						self.logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis." % gbk)
						self.logObject.warning(traceback.format_exc())
					raise RuntimeWarning(traceback.format_exc())

	def inputHomologyInformation(self, gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy):
		"""
		Simple function to store OrthoFinder homology information (parsed inititally in lsaBGC.utils)

		:param gene_to_hg: dictionary mapping gene locus tags to homolog group identifiers
		:param hg_genes: dictionary of sets with keys corresponding to homolog groups and values corresponding to set
						 of gene locus tags belonging to that homolog group.
		:param hg_median_copy_count: dictionary for the median copy count of each homolog group
		:param hg_prop_multi_copy: dictionary for the proportion of samples with homolog group which have multipmle
								   genes assigned to homolog group (paralogs).
		"""
		try:
			self.gene_to_hg = dict(gene_to_hg)
			self.hg_genes = dict(hg_genes)
			self.hg_median_copy_count = dict(hg_median_copy_count)
			self.hg_prop_multi_copy = dict(hg_prop_multi_copy)

			for bgc_id in self.bgc_genes:
				for gene in self.bgc_genes[bgc_id]:
					if gene in self.gene_to_hg:
						self.bgc_hgs[bgc_id].add(self.gene_to_hg[gene])

			if self.logObject:
				self.logObject.info(
					"Successfully inputted homolog information from OrthoFinder into variables of Pan object!")
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had issues inputting homology information into Pan object instance.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def openStatsFile(self, outdir, run_parameter_tests=False):
		"""
		Simple function to initialize final report file with statistics on GCF clustering - main output from lsaBGC-Cluster.

		:param outdir: path to workspace directory
		:param run_parameter_tests: whether lsaBGC-Cluster is being run in "run_parameter_tests mode".
		"""
		try:
			final_stats_file = outdir + 'GCF_Details.txt'
			final_stats_expanded_singletons_file = outdir + 'GCF_Details_Expanded_Singletons.txt'
			sf_handle = open(final_stats_file, 'w')
			if run_parameter_tests:
				sf_handle.write('\t'.join(
					['MCL inflation parameter', 'Jaccard similarity cutoff', 'GCF id', 'number of BGCs',
					 'number of samples',
					 'samples with multiple BGCs in GCF', 'number of core OGs', 'mean number of OGs',
					 'stdev for number of OGs', 'number of core gene aggregates',
					 'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'annotations']) + '\n')
			else:
				sfes_handle = open(final_stats_expanded_singletons_file, 'w')
				sf_handle.write(
					'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
							   'number of core OGs', 'mean number of OGs', 'stdev for number of OGs',
							   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
				sfes_handle.write(
					'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
							   'number of core OGs', 'mean number of OGs', 'stdev for number of OGs',
							   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
				sfes_handle.close()

			sf_handle.close()
			self.final_stats_file = final_stats_file
			self.final_stats_expanded_singletons_file = final_stats_expanded_singletons_file
			if self.logObject:
				self.logObject.info(
					"Will be writing final stats report on GCF clustering to file %s" % final_stats_file)
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had issues initializing final stats report file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def calculateBGCPairwiseRelations(self, outdir, split_by_annotation=False):
		"""
		Function to calculate the Jaccard Similarity between pairs of BGCs based on homolog groups shared.

		:param outdir: path to workspace directory.
		:param split_by_annotation: partition BGCs foremost by general annotation category of BGC by AntiSMASH (similar
		                            to what BiG-Scape performs.
		"""
		try:
			pair_relations_txt_file = outdir + 'bgc_pair_relationships.tmp.txt'
			prf_handle = open(pair_relations_txt_file, 'w')

			pairwise_relations = defaultdict(lambda: defaultdict(float))
			pairwise_containment_relations = defaultdict(lambda: defaultdict(float))
			pairwise_syntenic_relations = defaultdict(lambda: defaultdict(lambda: 0.0))

			bgc_hg_ranked_orders = defaultdict(dict)
			for bgc in self.bgc_hgs:
				sc_hg_starts = []
				for hg in self.bgc_hgs[bgc]:
					gene_counts = 0
					hg_start = 0
					for g in self.hg_genes[hg]:
						if self.comp_gene_info[g]['bgc_name'] == bgc:
							hg_start = self.comp_gene_info[g]['start']
							gene_counts += 1
					if gene_counts == 1:
						sc_hg_starts.append([hg, hg_start])

				for hg_ord, hg_its in enumerate(sorted(sc_hg_starts, key=itemgetter(1))):
					bgc_hg_ranked_orders[bgc][hg_its[0]] = (hg_ord + 1)

			for i, bgc1 in enumerate(self.bgc_hgs):
				bgc1_hgs = self.bgc_hgs[bgc1]
				bgc1_hgs_ranked_orders = bgc_hg_ranked_orders[bgc1]
				for j, bgc2 in enumerate(self.bgc_hgs):
					if i < j:
						bgc2_hgs = self.bgc_hgs[bgc2]
						bgc2_hgs_ranked_orders = bgc_hg_ranked_orders[bgc2]

						overlap_metric = float(len(bgc1_hgs.intersection(bgc2_hgs))) / float(len(bgc1_hgs.union(bgc2_hgs)))
						overlap_metric_scaled = 100.00 * overlap_metric
						pairwise_relations[bgc1][bgc2] = overlap_metric_scaled
						pairwise_relations[bgc2][bgc1] = overlap_metric_scaled

						bgc1_containment_in_bgc2 = (float(len(bgc1_hgs.intersection(bgc2_hgs))) / float(len(bgc1_hgs)))*100.0
						bgc2_containment_in_bgc1 = (float(len(bgc1_hgs.intersection(bgc2_hgs))) / float(len(bgc2_hgs)))*100.0
						pairwise_containment_relations[bgc1][bgc2] = bgc1_containment_in_bgc2
						pairwise_containment_relations[bgc2][bgc1] = bgc2_containment_in_bgc1

						sc_hgs_intersect = set(bgc1_hgs_ranked_orders.keys()).intersection(set(bgc2_hgs_ranked_orders.keys()))
						if len(sc_hgs_intersect) >= 3:
							bgc1_orders = []
							bgc2_orders = []
							for hg in sc_hgs_intersect:
								bgc1_orders.append(bgc1_hgs_ranked_orders[hg])
								bgc2_orders.append(bgc2_hgs_ranked_orders[hg])
							r, pval = stats.pearsonr(bgc1_orders, bgc2_orders)
							abs_correlation_coefficient = abs(r)
							pairwise_syntenic_relations[bgc1][bgc2] = abs_correlation_coefficient
							pairwise_syntenic_relations[bgc2][bgc1] = abs_correlation_coefficient

						products_union_not_empty = len(set(self.bgc_product[bgc1]).union(set(self.bgc_product[bgc2]))) > 0
						products_intersection_matches = float(len(set(self.bgc_product[bgc1]).intersection(set(self.bgc_product[bgc2]))))/float(len(set(self.bgc_product[bgc1]).union(set(self.bgc_product[bgc2])))) == 1.0
						if (not split_by_annotation) or (split_by_annotation and products_union_not_empty and products_intersection_matches):
							prf_handle.write('%s\t%s\t%f\n' % (bgc1, bgc2, overlap_metric_scaled))
			prf_handle.close()

			if self.logObject:
				self.logObject.info("Calculated pairwise relations and wrote to: %s" % pair_relations_txt_file)
			self.pairwise_relations = pairwise_relations
			self.pairwise_containment_relations = pairwise_containment_relations
			self.pairwise_syntenic_relations = pairwise_syntenic_relations
			self.pair_relations_txt_file = pair_relations_txt_file
			self.bgc_to_gcf_map_file = outdir + 'BGC_to_GCF_Mapping.txt'
		except Exception as e:
			if self.logObject:
				self.logObject.error("Problem with creating relations file between pairs of BGCs.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runMCLAndReportGCFs(self, mip, jcp, sccp, ccp, outdir, run_parameter_tests=False, threads=1):
		"""
		Function to run MCL and report the GCFs (gene-cluster families) of homologous BGCs identified.

		:param mip: MCL inflation parameter.
		:param jcp: Jaccard similarity threshold for homology between two BGCs to be considered.
		:param outdir: path to workspace directory.
		:param run_parameter_tests: True
		:param threads: number of threads/threads to use for MCL.
		"""
		pair_relations_filt_txt_file = outdir + 'bgc_pair_relationships.txt' 
		try:
			prftf_handle = open(pair_relations_filt_txt_file, 'w')
			with open(self.pair_relations_txt_file) as oprtf:
				for line in oprtf:
					line = line.strip('\n')
					b1, b2, jaccard_sim = line.split('\t')
					if self.pairwise_containment_relations[b1][b2] >= sccp:
						if float(jaccard_sim) >= jcp:
							prftf_handle.write(line + '\n')
						elif self.bgc_edgy[b1] == True and self.pairwise_containment_relations[b1][b2] >= ccp:
							prftf_handle.write(line + '\n')
						elif self.bgc_edgy[b2] == True and self.pairwise_containment_relations[b2][b1] >= ccp:
							prftf_handle.write(line + '\n')
			prftf_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.warning("Problem with parsing paired sample Jaccard similarity relationship file.")
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning(traceback.format_exc())
			return
		
		pair_relations_mci_file = outdir + 'bgc_pair_relationships.mci'
		pair_relations_tab_file = outdir + 'bgc_pair_relationships.tab'
		relations_mcl_file = outdir + 'mcl.' + str(mip).replace('.', '_') + '.out'
		mcxdump_out_file = outdir + 'final_mcl.' + str(mip).replace('.', '_') + '.out'

		mcxload_cmd = ['mcxload', '-abc', pair_relations_filt_txt_file, '--stream-mirror', '-write-tab',
					   pair_relations_tab_file,
					   '-o', pair_relations_mci_file]
		mcl_cmd = ['mcl', pair_relations_mci_file, '-I', str(mip), '-o', relations_mcl_file, '-te', str(threads)]
		mcxdump_cmd = ['mcxdump', '-icl', relations_mcl_file, '-tabr', pair_relations_tab_file, '-o',
					   mcxdump_out_file]

		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(mcxload_cmd))
		try:
			subprocess.call(' '.join(mcxload_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(mcxload_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.warning('Had an issue running: %s' % ' '.join(mcxload_cmd))
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning('Had an issue running: %s' % ' '.join(mcxload_cmd))
			return
		
		if self.logObject:
			self.logObject.info('Converted format of pair relationship file via mxcload.')

		if self.logObject:
			self.logObject.info('Running MCL and MCXDUMP with inflation parameter set to %f' % mip)
			self.logObject.info('Running the following command: %s' % ' '.join(mcl_cmd))
		try:
			subprocess.call(' '.join(mcl_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(relations_mcl_file) and os.path.getsize(relations_mcl_file) > 50)
			self.logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.warning('Had an issue running: %s' % ' '.join(mcl_cmd))
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning('Had an issue running: %s' % ' '.join(mcl_cmd))
			return
		
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(mcxdump_cmd))
		try:
			subprocess.call(' '.join(mcxdump_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(mcxdump_out_file) and os.path.getsize(mcxdump_out_file) > 50)
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.warning('Had an issue running: %s' % ' '.join(mcxdump_cmd))
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning('Had an issue running: %s' % ' '.join(mcxdump_cmd))
			return
		
		if self.logObject:
			self.logObject.info('Successfully ran MCL and MCXDUMP with inflatimulti_same_sampleon parameter set to %f!' % mip)

		try:
			sf_handle = open(self.final_stats_file, 'a+')
			sfes_handle = open(self.final_stats_expanded_singletons_file, 'a+')

			clustered_bgcs = set([])
			gcf_identifier = 1
			with open(mcxdump_out_file) as omo:
				for gcf in omo:
					gcf = gcf.strip()
					gcf_mems = gcf.split()
					if len(gcf_mems) < 2: continue
					diffs = set([])
					samp_counts = defaultdict(int)
					samp_ogs = defaultdict(set)
					products = defaultdict(float)
					core_gene_cluster_counts = 0
					for a, bgc1 in enumerate(gcf_mems):
						samp1 = self.bgc_sample[bgc1]
						samp_counts[samp1] += 1
						for prod in self.bgc_product[bgc1]: products[prod] += 1.0 / len(self.bgc_product[bgc1])
						if core_gene_cluster_counts < self.bgc_core_counts[bgc1]: core_gene_cluster_counts = self.bgc_core_counts[bgc1]
						samp_ogs[samp1] = samp_ogs[samp1].union(self.bgc_hgs[bgc1])
						clustered_bgcs.add(bgc1)
						for b, bgc2 in enumerate(gcf_mems):
							if a < b:
								diffs.add(self.pairwise_relations[bgc1][bgc2])
					multi_same_sample = 0
					num_ogs = []
					core_hgs = set([])
					for si, s in enumerate(samp_counts):
						if samp_counts[s] > 1:
							multi_same_sample += 1
						if si == 0:
							core_hgs = samp_ogs[s]
						else:
							core_hgs = core_hgs.intersection(samp_ogs[s])
						num_ogs.append(len(samp_ogs[s]))
					stdev = "NA"
					mean = "NA"
					try:
						mean = statistics.mean(num_ogs)
						stdev = statistics.stdev(num_ogs)
					except:
						pass
					gcf_stats = ['GCF_' + str(gcf_identifier), len(gcf_mems), len(samp_counts.keys()), multi_same_sample,
								 len(core_hgs),
								 mean, stdev, min(diffs),
								 max(diffs), core_gene_cluster_counts,
								 '; '.join([x[0] + ':' + str(x[1]) for x in products.items()])]
					if run_parameter_tests:
						gcf_stats = [mip, jcp] + gcf_stats
					else:
						sfes_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
					sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
					gcf_identifier += 1

			singleton_bgcs = []
			for bgc in self.bgc_hgs:
				if not bgc in clustered_bgcs:
					singleton_bgcs.append(bgc)
					if not run_parameter_tests:
						products = {}
						for prod in self.bgc_product[bgc]: products[prod] = 1.0 / len(self.bgc_product[bgc])
						mean_og_count = len(self.bgc_hgs[bgc])
						gcf_stats = ['GCF_' + str(gcf_identifier), 1, 1, 0, "NA", mean_og_count, 0, 0, 0, 0,  '; '.join([x[0] + ':' + str(x[1]) for x in products.items()]) ]
						sfes_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
						gcf_identifier += 1

			singleton_stats = ['singletons', len(singleton_bgcs)] + (['NA'] * 7)
			if run_parameter_tests:
				singleton_stats = [mip, jcp] + singleton_stats
			sf_handle.write('\t'.join([str(x) for x in singleton_stats]) + '\n')

			sf_handle.close()
			sfes_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.warning(
					"Problem appending information on GCFs for current parameter combination to GCF statistics report file.")
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning(traceback.format_exc())
			return
		
		if not run_parameter_tests:
			try:
				if self.logObject:
					self.logObject.info(
					"Writing list of BGCs for each GCF, which will be used as input for downstream programs in the suite!")
				gcf_listing_dir = outdir + 'GCF_Listings/'
				if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
				gcf_identifier = 1
				with open(mcxdump_out_file) as omo:
					for gcf in omo:
						gcf = gcf.strip()
						gcf_mems = gcf.split()
						if len(gcf_mems) < 2: continue
						outf_list = open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w')
						for bgc in gcf_mems:
							outf_list.write('%s\n' % (self.bgc_gbk[bgc]))
						outf_list.close()
						gcf_identifier += 1
				for bgc in singleton_bgcs:
					outf_list = open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w')
					outf_list.write('%s\n' % (self.bgc_gbk[bgc]))
					gcf_identifier += 1

				if self.logObject:
					self.logObject.info("Successfully wrote lists of BGCs for each GCF.")
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing GCF lists.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError(traceback.format_exc())
		else:
			try:
				if self.logObject:
					self.logObject.info(
					"Writing list of BGCs for each GCF for each clustering parameter combination into single plot.")
				btgmf_handle = open(self.bgc_to_gcf_map_file, 'a+')

				with open(mcxdump_out_file) as omo:
					for j, gcf in enumerate(omo):
						gcf = gcf.strip()
						gcf_mems = gcf.split()
						if len(gcf_mems) < 2: continue
						for bgc in gcf_mems:
							sname = self.bgc_sample[bgc]
							btgmf_handle.write('%f\t%f\tGCF_%d\t%s\t%s\n' % (mip, jcp, j + 1, sname, bgc))
				for sbgc in singleton_bgcs:
					sname = self.bgc_sample[sbgc]
					btgmf_handle.write('%f\t%f\tGCF_singletons\t%s\t%s\n' % (mip, jcp, sname, sbgc))
				btgmf_handle.close()
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing BGC to GCF relationship.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError(traceback.format_exc())

	def plotResultsFromUsingDifferentParameters(self, outdir):
		"""
		Function to create an 8x11 in PDF with figures aimed to give the user an idea on how different values of the
		MCL inflation parameter and the Jaccard similarity cutoff might impact clustering to help them choose the best
		parameter values.

		:param outdir: path to workspace directory.
		"""
		try:
			plot_input_dir = outdir + 'plotting_input/'
			if not os.path.isdir(plot_input_dir): os.system('mkdir %s' % plot_input_dir)

			singleton_counts = defaultdict(int)
			clustered_counts = defaultdict(int)
			cluster_counts = defaultdict(int)
			cluster_mixedannot_counts = defaultdict(int)
			cluster_mixedannot_wohypo_counts = defaultdict(int)
			all_annotation_classes = set([])
			with open(self.final_stats_file) as ogdf:
				for i, line in enumerate(ogdf):
					line = line.strip()
					ls = line.split("\t")
					if i == 0: continue
					mcl = round(float(ls[0]), 2)
					jcp = round(float(ls[1]), 2)
					param_id = '%f_%f' % (mcl, jcp)
					gcf_id = ls[2]
					if gcf_id == 'singletons':
						singleton_counts[param_id] = int(ls[3])
					else:
						clustered_counts[param_id] += int(ls[3])
						cluster_counts[param_id] += 1
						annotations = set(
							[x.split(':')[0].replace('-', '.') for x in ls[-1].split('; ') if
							 x.split(':')[0] != 'NA'])
						if len(ls[-1].split('; ')) > 1: cluster_mixedannot_counts[param_id] += 1
						if len(annotations) > 1: cluster_mixedannot_wohypo_counts[param_id] += 1
						all_annotation_classes = all_annotation_classes.union(annotations)

			plot_overview_file = plot_input_dir + 'plotting_input_1.txt'
			plot_input_1_file = plot_input_dir + 'plotting_input_2.txt'
			plot_input_2_file = plot_input_dir + 'plotting_input_3.txt'
			plot_annot_file = plot_input_dir + 'plotting_input_4.txt'
			plot_input_3_file = plot_input_dir + 'plotting_input_5.txt'

			plot_input_1_handle = open(plot_input_1_file, 'w')
			plot_input_2_handle = open(plot_input_2_file, 'w')
			plot_annot_handle = open(plot_annot_file, 'w')
			plot_overview_handle = open(plot_overview_file, 'w')
			plot_input_3_handle = open(plot_input_3_file, 'w')

			plot_annot_handle.write('\n'.join(['Unknown'] + sorted(list(all_annotation_classes))))
			plot_annot_handle.close()

			coord_tuples = defaultdict(set)
			plot_input_1_handle.write('\t'.join(
				['GCF', 'Parameters', 'JaccardSim', 'Inflation', 'Samples', 'CoreExists', 'CoreSize', 'CoreGeneClusters', 'AvgGeneCount', 'StdDevGeneCount', 'Unknown'] + sorted(list(all_annotation_classes))) + '\n')
			plot_input_2_handle.write('\t'.join(['GCF', 'Parameters', 'Annotation', 'Count', 'Total.BGCs']) + '\n')
			plot_input_3_handle.write('\t'.join(['GCF', 'Parameters', 'Total.BGCs', 'Category', 'Sample_Count']) + '\n')

			plot_overview_handle.write('\t'.join(
				['Parameters', 'JaccardSim', 'Inflation', 'Clustered', 'Singletons', 'SingletonToClustered',
				 'MixedAnnotations', 'MixedAnnotationsWoHypothetical', 'NumberClusters',
				 'MixedAnnotationProportion',
				 'MixedAnnotationWoHypotheticalProportion']) + '\n')

			jcp_gcfs = defaultdict(lambda: defaultdict(set))
			with open(self.bgc_to_gcf_map_file) as obgmf:
				for line in obgmf:
					line = line.strip()
					mip, jcp, gcf_id, sname, gbkpath = line.split('\t')
					jcp = float(jcp)
					mip = float(mip)
					param_mod_id = 'Inf: %.1f JacSim: %.1f' % (mip, jcp)
					gcf_name = gcf_id + ' - ' + param_mod_id
					jcp_gcfs[jcp][gcf_name].add(gbkpath)
						
			with open(self.final_stats_file) as ogdf:
				for i, line in enumerate(ogdf):
					line = line.strip()
					ls = line.split("\t")
					if i == 0: continue
					mcl = round(float(ls[0]), 2)
					jcp = round(float(ls[1]), 2)
					param_id = '%f_%f' % (mcl, jcp)
					gcf_id = ls[2]
					if gcf_id == 'singletons': continue

					param_mod_id = 'MCL Inflation: %.1f; Jaccard Similarity Cutoff: %.1f' % (mcl, jcp)
					# get stats for plotting

					samples = int(ls[4])
					samples_with_multi_bgcs = int(ls[5])
					size_of_core = int(ls[6])
					number_of_core_genes = ls[-2]
					avg_gene_count = ls[-6]
					stdev_gene_count = ls[-5]
					annotations = ls[-1].split('; ')

					unique_coord = False
					random_offset_1 = random.randint(-50, 51) / 100.0
					random_offset_2 = random.randint(-50, 51) / 100.0
					coord_tuple = tuple([samples + random_offset_1, samples_with_multi_bgcs + random_offset_2])
					while not unique_coord:
						if not coord_tuple in coord_tuples[param_mod_id]:
							unique_coord = True
						else:
							random_offset_1 = random.randint(-50, 51) / 100.0
							random_offset_2 = random.randint(-50, 51) / 100.0
							coord_tuple = tuple(
								[samples + random_offset_1, samples_with_multi_bgcs + random_offset_2])
					coord_tuples[param_mod_id].add(coord_tuple)
					core_exists = "Core DNE"
					if size_of_core > 0:
						core_exists = "Core Exists"

					annot_count = defaultdict(float)
					for an in annotations:
						ans = an.split(':')
						annot_count[ans[0].replace('-', '.')] = float(ans[1])

					annotation_class_abd = [annot_count['NA']]
					for ac in sorted(all_annotation_classes):
						annotation_class_abd.append(annot_count[ac])

					plot_input_1_handle.write('\t'.join([str(x) for x in
														 [gcf_id + ' - ' + param_mod_id, param_mod_id, jcp, mcl,
														  samples, core_exists,
														  size_of_core, number_of_core_genes, avg_gene_count,
														  stdev_gene_count] + annotation_class_abd]) + '\n')
	
					tot_count = sum(annot_count.values())
					plot_input_3_handle.write('\t'.join([str(x) for x in [gcf_id + ' - ' + param_mod_id, param_mod_id, tot_count, 'Single-Copy', (samples-samples_with_multi_bgcs)]]) + '\n')
					plot_input_3_handle.write('\t'.join([str(x) for x in [gcf_id + ' - ' + param_mod_id, param_mod_id, tot_count, 'Multi-Copy', samples_with_multi_bgcs]]) + '\n')
	
					if int(gcf_id.split('_')[1]) <= 50:
						for an in annot_count:
							an_count = annot_count[an]
							if an == 'NA':
								an = 'Unknown'
							plot_input_2_handle.write('\t'.join([str(x) for x in
																 [gcf_id + ' - ' + param_mod_id, param_mod_id, an,
																  an_count, tot_count]]) + '\n')

			for pid in clustered_counts:
				mcl, jcp = [float(x) for x in pid.split('_')]
				param_mod_id = 'MCL Inflation: %.1f; Jaccard Similarity Cutoff: %.1f' % (mcl, jcp)

				plot_overview_handle.write(
					'\t'.join([str(x) for x in [param_mod_id, jcp, mcl, clustered_counts[pid],
												singleton_counts[pid],
												singleton_counts[pid] / clustered_counts[pid],
												cluster_mixedannot_counts[pid],
												cluster_mixedannot_wohypo_counts[pid],
												cluster_counts[pid],
												cluster_mixedannot_counts[pid] / float(
													cluster_counts[pid]),
												cluster_mixedannot_wohypo_counts[pid] / float(
													cluster_counts[pid])]]) + '\n')

				for sgcf_id in range(1, singleton_counts[pid] + 1):
					sgcf_name = 'SGCF_' + str(sgcf_id)
					plot_input_1_handle.write('\t'.join(
						[sgcf_name + ' - ' + param_mod_id, param_mod_id, str(jcp), str(mcl), '1',
						 'Not Applicable'] + ['NA'] * (5 + len(all_annotation_classes))) + '\n')

			plot_input_1_handle.close()
			plot_input_2_handle.close()
			plot_input_3_handle.close()
			plot_overview_handle.close()

			plot_pdf_file = outdir + 'Plots_Depicting_Parameter_Influence_on_GCF_Clustering.pdf'
			
			self.createParamInfluencePlot(plot_input_dir, plot_input_1_file, plot_annot_file, plot_input_2_file, plot_overview_file, plot_input_3_file, plot_pdf_file)

			if self.logObject:
				self.logObject.info('Plotting completed!')

		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Problem creating plot(s) for assessing best parameter choices for GCF clustering.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def convertGenbanksIntoFastas(self, fasta_dir, fasta_listing_file, sample_retention_set=None):
		"""
		Function to convert Genbanks for BGC instances into FASTA format and listing files. Note, there will
		be one FASTA per sample not per BGC Genbank (in case sample's have more than one associated Genbank).

		:param listing_file: tab-delimited file with two columns: (1) sample name (2) path to BGC Genbank
		:param fasta_listing_file: tab-delimited file with two columns: (1) sample name (2) path to BGC FASTA
		"""
		try:
			sample_index = defaultdict(int)

			all_samples = set([])
			with open(self.bgc_genbanks_listing) as obsf:
				for i, line in enumerate(obsf):
					line = line.strip()
					try:
						assert (len(line.split('\t')) >= 2)
					except Exception as e:
						if self.logObject:
							self.logObject.error(
							"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
									i + 1))
							self.logObject.error(traceback.format_exc())
						raise RuntimeError(
							"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
									i + 1))
					if len(line.split('\t')) == 2:
						sample, gbk = line.split('\t')
					else:
						sample, gbk = line.split('\t')[:-1]
					sample = util.cleanUpSampleName(sample)
					if sample_retention_set != None and sample not in sample_retention_set: continue
					try:
						assert (util.is_genbank(gbk))
						bgc_id = sample
						if sample_index[sample] > 0:
							bgc_id = sample + '_' + str(sample_index[sample] + 1)
						sample_index[sample] += 1

						sample_fasta = fasta_dir + sample + '.fasta'
						sample_handle = open(sample_fasta, 'a+')
						with open(gbk) as ogbk:
							for ri, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
								sample_handle.write('>' + bgc_id + '|' + str(ri) + '\n' + str(rec.seq) + '\n')
						sample_handle.close()
						all_samples.add(sample)
						if self.logObject:
							self.logObject.info("Added Genbank %s into sample %s's GCF relevant FASTA." % (gbk, sample))
					except Exception as e:
						if self.logObject:
							self.logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
							self.logObject.warning(traceback.format_exc())
						raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

			outf = open(fasta_listing_file, 'w')
			for s in all_samples:
				outf.write(s + '\t' + fasta_dir + s + '.fasta\n')
			outf.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues converting BGC Genbank listing into FASTA listing.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def createParamInfluencePlot(self, workspace, plot_input_1_file, plot_annot_file, plot_input_2_file, plot_overview_file, plot_input_3_file, plot_pdf_file):
		rscript_file = workspace + 'params_impact_rscript.R' 
		rfh = open(rscript_file, 'w')

		rfh.write('library(ggplot2)\n')
		rfh.write('library(RColorBrewer)\n')
		rfh.write('library(cowplot)\n\n')
		
		rfh.write('# read in plotting inputs\n')
		rfh.write('data <- read.table("' + plot_input_1_file + '", header=T, sep="\t")\n')
		rfh.write('annot_classes <- scan("' + plot_annot_file + '", what="", sep="\n")\n')
		rfh.write('largest_gcfs_data <- read.table("' + plot_input_2_file + '", header=T, sep="\t")\n')
		rfh.write('overview_data <- read.table("' + plot_overview_file + '", header=T, sep="\t")\n')
		rfh.write('cc_gcfs_data <- read.table("' + plot_input_3_file + '", header=T, sep="\t")\n')
		rfh.write('pdf_file <- "' + plot_pdf_file + '"\n\n')

		rfh.write('# define color scheme for plots\n')
		rfh.write('nb.cols <- length(annot_classes)-1\n')
		rfh.write('annot_colors <- colorRampPalette(brewer.pal(nb.cols, "Set2"))(nb.cols)\n')
		rfh.write('annot_colors <- c("#808080", annot_colors)\n')
		rfh.write('names(annot_colors) <- annot_classes\n')
		rfh.write('core_colors <- c("#737574", "#222422", "#ed665f")\n')
		rfh.write('names(core_colors) <- c("Not Applicable", "Core Exists", "Core DNE")\n')

		rfh.write('# open pdf for the report\n')
		rfh.write('pdf(pdf_file, height=8, width=11)\n')

		rfh.write('\n')
		rfh.write('# Title page + overview scattermap\n')
		rfh.write('repname <- ggplot()+theme_void() + ggtitle("Parameter Impact on\nGCF Delineations") + theme(plot.title = element_text(hjust = 0.5, size=70, face="bold"))\n')
		rfh.write('scatter_plot <- ggplot(overview_data, aes(x=NumberClusters, y=SingletonToClustered)) +\n')
		rfh.write('  geom_jitter(aes(color=as.factor(Inflation), shape=as.factor(JaccardSim)), alpha=0.7, size=2) + theme_classic() +\n')
		rfh.write('  ggtitle("Tradeoff between Cluster Granularity vs. Singletons Incurred") + xlab("Number of GCFs") +\n')
		rfh.write('  ylab("Ratio of BGCs which are\nsingletons to those in GCFs") +\n')
		rfh.write('  guides(shape=guide_legend(title="Jaccard Similarity Cutoff"), color=guide_legend(title="MCL Inflation"))+\n')
		rfh.write('  theme(legend.position = "bottom")\n')
		rfh.write('print(plot_grid(repname, scatter_plot, rel_heights=c(1,2), ncol=1))\n')
		rfh.write('\n')
		rfh.write('# Overview Heatmaps\n')
		rfh.write('gheat1 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=SingletonToClustered)) +\n')
		rfh.write('  geom_tile() + theme_classic() +  scale_fill_gradient(low = "lightblue", high = "darkblue") + xlab("MCL Inflation") +\n')
		rfh.write('  ylab("Jaccard Similarity Cutoff") + ggtitle("Ratio of BGCs which are\nsingletons to those in GCFs")+ theme(legend.title = element_blank())\n')
		rfh.write('gheat2 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=NumberClusters)) +\n')
		rfh.write('  geom_tile() + theme_classic() +   scale_fill_gradient(low = "lightgreen", high = "darkgreen") + xlab("MCL Inflation") +\n')
		rfh.write('  ylab("Jaccard Similarity Cutoff") + ggtitle("Number of GCFs")+ theme(legend.title = element_blank())\n')
		rfh.write('gheat3 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=MixedAnnotationProportion)) +\n')
		rfh.write('  geom_tile() + theme_classic() +   scale_fill_gradient(low = "yellow", high = "gold") + xlab("MCL Inflation") +\n')
		rfh.write('  ylab("Jaccard Similarity Cutoff") + ggtitle("Proportion of GCFs which feature\nmultiple BGC annotation categories\n(including unannotated)")+ theme(legend.title = element_blank())\n')
		rfh.write('gheat4 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=MixedAnnotationWoHypotheticalProportion)) +\n')
		rfh.write('  geom_tile() + theme_classic() +   scale_fill_gradient(low = "pink", high = "darkred") + xlab("MCL Inflation") +\n')
		rfh.write('  ylab("Jaccard Similarity Cutoff") + ggtitle("Proportion of GCFs which feature\nmultiple BGC annotation categories\n(not including unannotated)")+ theme(legend.title = element_blank())\n')
		rfh.write('\n')
		rfh.write('gheatmaps_top <- plot_grid(gheat1, gheat2)\n')
		rfh.write('gheatmaps_bottom <- plot_grid(gheat3, gheat4)\n')
		rfh.write('print(plot_grid(gheatmaps_top, gheatmaps_bottom, rel_heights=c(1,1), ncol=1))\n')
		rfh.write('\n')
		rfh.write('# Boxplot views\n')
		rfh.write('boxpl1 <- ggplot(data, aes(x=as.factor(Inflation), y=CoreSize, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Number of Core Homolog Groups per GCF") + scale_fill_brewer(palette="Set2")# + scale_y_log10()\n')
		rfh.write('boxpl2 <- ggplot(data, aes(x=as.factor(Inflation), y=Samples, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Number of Samples with GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()\n')
		rfh.write('boxpl3 <- ggplot(data, aes(x=as.factor(Inflation), y=AvgGeneCount, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Avg. Gene Count per BGC in GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()\n')
		rfh.write('boxpl4 <- ggplot(data, aes(x=as.factor(Inflation), y=StdDevGeneCount, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Std. Deviation of Gene Count for BGCs in GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()\n')
		rfh.write('\n')
		rfh.write('top_row <- plot_grid(boxpl1, boxpl2)\n')
		rfh.write('bot_row <- plot_grid(boxpl3, boxpl4)\n')
		rfh.write('\n')
		rfh.write('print(plot_grid(boxpl1, boxpl2, boxpl3, boxpl4, rel_heights=c(1,1,1,1), ncol=1))\n')
		rfh.write('\n')
		rfh.write('# determine number of parameter combinations tested and axes boundaries for scatter plots\n')
		rfh.write('params <- unique(data$Parameters)\n')
		rfh.write('npages <- length(params)\n')
		rfh.write('\n')
		rfh.write('for(i in 1:npages) {\n')
		rfh.write('  curr_param = params[i]\n')
		rfh.write('\n')
		rfh.write('  data_filt <- data[data$Parameters == curr_param,]\n')
		rfh.write('  data_filt_scatter <- data_filt[!is.na(data_filt$Samples_Offset),]\n')
		rfh.write('  largest_gcfs_data_filt <- largest_gcfs_data[largest_gcfs_data$Parameters == curr_param,]\n')
		rfh.write('  cc_gcfs_data_filt <- cc_gcfs_data[cc_gcfs_data$Parameters == curr_param,]\n')
		rfh.write('\n')
		rfh.write('  paramname <- ggplot()+theme_void() + ggtitle(curr_param) + theme(plot.title = element_text(hjust = 0.5, size=32, face="bold"))\n')
		rfh.write('\n')
		rfh.write('  ghist_1 <- ggplot(data_filt, aes(x=Samples, fill=CoreExists)) + geom_histogram(color="black", size=0.3) + theme_classic() +\n')
		rfh.write('    theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values=core_colors) + xlab("Sample Count") + guides(fill=guide_legend(title=""))\n')
		rfh.write('\n')
		rfh.write('  ghist_2 <- ggplot(data_filt, aes(x=StdDevGeneCount)) + xlab("Std. Deviation in Gene Count") + geom_histogram(color="black", fill="black") + theme_classic() + scale_y_log10()\n')
		rfh.write('\n')
		rfh.write('  gbar_1 <- ggplot(cc_gcfs_data_filt, aes(x=reorder(GCF, Total.BGCs), y=Sample_Count, fill=Category)) +\n')
		rfh.write('    geom_bar(stat="identity", color="black", show.legend=F) + theme_classic() + scale_fill_manual(values=c("black", "grey")) +\n')
		rfh.write('    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("All GCFs (grey=single-copy, black=multi-copy)") + ylab("Sample Count")\n')
		rfh.write('\n')
		rfh.write('  gbar_2 <- ggplot(largest_gcfs_data_filt, aes(x=reorder(GCF, Total.BGCs), y=Count, fill=Annotation)) +\n')
		rfh.write('    geom_bar(stat="identity", color="black", show.legend=F) + theme_classic() + scale_fill_manual(values=annot_colors) +\n')
		rfh.write('    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Top 50 GCFs with Most BGCs (color = BGC product type)") + ylab("BGC Count")\n')
		rfh.write('\n')
		rfh.write('  mid_row <- plot_grid(ghist_1, ghist_2)\n')
		rfh.write('  print(plot_grid(paramname, mid_row, gbar_1, gbar_2, rel_heights=c(1,2,2,2), ncol=1))\n')
		rfh.write('}\n')
		rfh.write('\n')

		rfh.write('dev.off()\n')

		rfh.close()

		rscript_plot_cmd = ["Rscript", rscript_file]
		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
		try:
			subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=sys.stderr,
							stderr=sys.stderr,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
	

	def readInPopulationsSpecification(self, pop_specs_file, prune_set=None):
		"""
		Read in population specifications per sample and assign population to each BGC.

		:param pop_specs_file: path to file which has tab separated sample (1st column) and corresponding population (2nd column)
		"""
		try:
			self.bgc_population = defaultdict(lambda: "NA")
			self.sample_population = defaultdict(lambda: "NA")
			with open(pop_specs_file) as opsf:
				for i, line in enumerate(opsf):
					if i == 0 and line.startswith('name'): continue
					line = line.strip()
					sample, population = line.split('\t')
					if prune_set != None and not sample in prune_set: continue
					self.sample_population[sample] = population
					for bgc in self.sample_bgcs[sample]:
						self.bgc_population[bgc] = population
			self.logObject.info("Successfully parsed population specifications file. There are %d populations." % len(population))
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issue in parsing population specifications per sample and associating with each BGC.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

