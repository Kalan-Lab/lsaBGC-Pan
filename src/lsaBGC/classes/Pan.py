import os
import contextlib
import sys
import logging
import traceback
import subprocess
from collections import defaultdict
from Bio import SeqIO
from lsaBGC.classes.BGC import BGC
from lsaBGC import util
from lsaBGC import paramplots
import statistics
import random
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
		Definition:
		Function to parse file listing location of BGC Genbanks.
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF, or BGC object.
		- comprehensive_parsing (optional): flag specifying whether to perform comprehensive extraction of information from 
		                                    Genbanks. default is True.
		- prune_set (optional): set of samples to prune out and disregard BGCs from. default is None. 
		- edge_dist_cutoff (optional): the distance from scaffold/contig edge for a BGC to be regarded as incomplete. default 
		                               is 5000
		********************************************************************************************************************
		"""
		sample_index = defaultdict(int)
		with open(self.bgc_genbanks_listing) as obsf:
			for i, line in enumerate(obsf):
				line = line.strip()
				try:
					assert (len(line.split('\t')) == 4)
				except Exception as e:
					msg = "More or less than four columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1)
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
					raise RuntimeWarning("Unable to validate %s as Genbank. Skipping its incorporation into analysis." % gbk) from e

	def inputHomologyInformation(self, gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy):
		"""
		Definition:
		Simple function to store OrthoFinder homology information (parsed initially in lsaBGC.utils)
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF, or BGC object.
		- gene_to_hg: dictionary mapping gene locus tags to homolog group identifiers
		- hg_genes: dictionary of sets with keys corresponding to homolog groups and values corresponding to set
					of gene locus tags belonging to that homolog group.
		- hg_median_copy_count: dictionary for the median copy count of each homolog group
		- hg_prop_multi_copy: dictionary for the proportion of samples with homolog group which have multiple
							  genes assigned to homolog group (paralogs).
		********************************************************************************************************************
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
			raise RuntimeError("Had issues inputting homology information into Pan object instance.") from e

	def openStatsFile(self, outdir, run_parameter_tests=False):
		"""
		Definition:
		Simple function to initialize final report file with statistics on GCF clustering - main output from lsaBGC-Cluster.
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF or BGC object.
		- outdir: path to workspace directory
		- run_parameter_tests: whether lsaBGC-Cluster is being run in "run_parameter_tests mode".
		********************************************************************************************************************
		"""
		try:
			# Parameter testing results are written to their own table (13 columns) and final clustering results to
			# GCF_Details.txt / GCF_Details_Expanded_Singletons.txt (11 columns). Column order matches the order in 
			# which values are written in runMCLAndReportGCFs.
			final_stats_expanded_singletons_file = None
			if run_parameter_tests:
				final_stats_file = outdir + 'GCF_Details.Parameter_Tests.txt'
			else:
				final_stats_file = outdir + 'GCF_Details.txt'
				final_stats_expanded_singletons_file = outdir + 'GCF_Details_Expanded_Singletons.txt'
			with open(final_stats_file, 'w') as sf_handle:
				if run_parameter_tests:
					sf_handle.write('\t'.join(
						['MCL inflation parameter', 'Jaccard similarity cutoff', 'GCF id', 'number of BGCs',
						 'number of samples',
						 'samples with multiple BGCs in GCF', 'number of core OGs', 'mean number of OGs',
						 'stdev for number of OGs', 'min pairwise Jaccard similarity',
						 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
				else:
					with open(final_stats_expanded_singletons_file, 'w') as sfes_handle:
						sf_handle.write(
							'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
									   'number of core OGs', 'mean number of OGs', 'stdev for number of OGs',
									   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
						sfes_handle.write(
							'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
									   'number of core OGs', 'mean number of OGs', 'stdev for number of OGs',
									   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')

			self.final_stats_file = final_stats_file
			self.final_stats_expanded_singletons_file = final_stats_expanded_singletons_file
			if self.logObject:
				self.logObject.info(
					"Will be writing final stats report on GCF clustering to file %s" % final_stats_file)
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had issues initializing final stats report file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError("Had issues initializing final stats report file.") from e

	def calculateBGCPairwiseRelations(self, outdir, split_by_annotation=False):
		"""
		Function to calculate the Jaccard Similarity between pairs of BGCs based on homolog groups shared.

		:param outdir: path to workspace directory.
		:param split_by_annotation: partition BGCs foremost by general annotation category of BGC by AntiSMASH (similar
		                            to what BiG-Scape performs.
		"""
		try:
			pair_relations_txt_file = outdir + 'bgc_pair_relationships.tmp.txt'
			with open(pair_relations_txt_file, 'w') as prf_handle:

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
			raise RuntimeError("Problem with creating relations file between pairs of BGCs.") from e

	def runMCLAndReportGCFs(self, mip, jcp, sccp, ccp, outdir, run_parameter_tests=False, threads=1):
		"""
		Description:
		Function to run MCL and report the GCFs (gene-cluster families) of homologous BGCs identified.
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF or BGC object.
		- mip: MCL inflation parameter.
		- jcp: Jaccard similarity threshold for homology between two BGCs to be considered.
		- outdir: path to workspace directory.
		- run_parameter_tests: True
		- threads: number of threads/threads to use for MCL.
		********************************************************************************************************************
		"""
		pair_relations_filt_txt_file = outdir + 'bgc_pair_relationships.txt' 
		try:
			with open(pair_relations_filt_txt_file, 'w') as prftf_handle:
				pairs_kept = 0
				pairs_failing_synteny = 0
				with open(self.pair_relations_txt_file) as oprtf:
					for line in oprtf:
						line = line.strip('\n')
						b1, b2, jaccard_sim = line.split('\t')
						passes_overlap = (float(jaccard_sim) >= jcp or
										  (self.bgc_edgy[b1] == True and self.pairwise_containment_relations[b1][b2] >= ccp) or
										  (self.bgc_edgy[b2] == True and self.pairwise_containment_relations[b2][b1] >= ccp))
						if not passes_overlap: continue
						# require syntenic correspondence: absolute correlation of gene order for single-copy ortholog
						# groups shared between the pair (computed only when >= 3 such ortholog groups are shared,
						# otherwise it is 0.0).
						if self.pairwise_syntenic_relations[b1][b2] >= sccp:
							prftf_handle.write(line + '\n')
							pairs_kept += 1
						else:
							pairs_failing_synteny += 1
			if self.logObject:
				self.logObject.info('For MCL inflation %s and Jaccard cutoff %s: %d BGC pairs retained for MCL, %d pairs meeting overlap criteria were removed for not meeting the syntenic correlation cutoff of %s.' % (str(mip), str(jcp), pairs_kept, pairs_failing_synteny, str(sccp)))
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
			self.logObject.info('Converted format of pair relationship file via mcxload.')

		if self.logObject:
			self.logObject.info('Running MCL and MCXDUMP with inflation parameter set to %f' % mip)
			self.logObject.info('Running the following command: %s' % ' '.join(mcl_cmd))
		try:
			subprocess.call(' '.join(mcl_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(relations_mcl_file)) # and os.path.getsize(relations_mcl_file) > 50)
			if self.logObject:
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
			assert (os.path.isfile(mcxdump_out_file)) # and os.path.getsize(mcxdump_out_file) > 50)
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(mcxdump_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.warning('Had an issue running: %s' % ' '.join(mcxdump_cmd))
				self.logObject.warning(traceback.format_exc())
			#raise RuntimeWarning('Had an issue running: %s' % ' '.join(mcxdump_cmd))
			return
		
		if self.logObject:
			self.logObject.info('Successfully ran MCL and MCXDUMP with inflation parameter set to %f!' % mip)

		try:
			sfes_context = contextlib.nullcontext()
			if not run_parameter_tests:
				sfes_context = open(self.final_stats_expanded_singletons_file, 'a+')
			with open(self.final_stats_file, 'a+') as sf_handle, sfes_context as sfes_handle:

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
						except Exception as e:
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

				singleton_stats = ['singletons', len(singleton_bgcs)] + (['NA'] * 9)
				if run_parameter_tests:
					singleton_stats = [mip, jcp] + singleton_stats
				sf_handle.write('\t'.join([str(x) for x in singleton_stats]) + '\n')

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
				os.makedirs(gcf_listing_dir, exist_ok=True)
				gcf_identifier = 1
				with open(mcxdump_out_file) as omo:
					for gcf in omo:
						gcf = gcf.strip()
						gcf_mems = gcf.split()
						if len(gcf_mems) < 2: continue
						with open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w') as outf_list:
							for bgc in gcf_mems:
								outf_list.write('%s\n' % (self.bgc_gbk[bgc]))
						gcf_identifier += 1
				for bgc in singleton_bgcs:
					with open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w') as outf_list:
						outf_list.write('%s\n' % (self.bgc_gbk[bgc]))
					gcf_identifier += 1

				if self.logObject:
					self.logObject.info("Successfully wrote lists of BGCs for each GCF.")
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing GCF lists.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError("Problem with writing GCF lists.") from e
		else:
			try:
				if self.logObject:
					self.logObject.info(
					"Writing list of BGCs for each GCF for each clustering parameter combination into single plot.")
				with open(self.bgc_to_gcf_map_file, 'a+') as btgmf_handle:

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
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing BGC to GCF relationship.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError("Problem with writing BGC to GCF relationship.") from e

	def plotResultsFromUsingDifferentParameters(self, outdir):
		"""
		Description:
		Function to create an 8x11 in PDF with figures aimed to give the user an idea on how different values of the
		MCL inflation parameter and the Jaccard similarity cutoff might impact clustering to help them choose the best
		parameter values.
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF or BGC object.
		- outdir: path to workspace directory.
		********************************************************************************************************************
		"""
		try:
			plot_input_dir = outdir + 'plotting_input/'
			os.makedirs(plot_input_dir, exist_ok=True)

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

			with open(plot_input_1_file, 'w') as plot_input_1_handle, \
					open(plot_input_2_file, 'w') as plot_input_2_handle, \
					open(plot_annot_file, 'w') as plot_annot_handle, \
					open(plot_overview_file, 'w') as plot_overview_handle, \
					open(plot_input_3_file, 'w') as plot_input_3_handle:

				plot_annot_handle.write('\n'.join(['Unknown'] + sorted(list(all_annotation_classes))))

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


			plot_pdf_file = outdir + 'Plots_Depicting_Parameter_Influence_on_GCF_Clustering.pdf'
			
			self.createParamInfluencePlot(plot_input_dir, plot_input_1_file, plot_annot_file, plot_input_2_file, plot_overview_file, plot_input_3_file, plot_pdf_file)

			if self.logObject:
				self.logObject.info('Plotting completed!')

		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Problem creating plot(s) for assessing best parameter choices for GCF clustering.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError("Problem creating plot(s) for assessing best parameter choices for GCF clustering.") from e

	def createParamInfluencePlot(self, workspace, plot_input_1_file, plot_annot_file, plot_input_2_file, plot_overview_file, plot_input_3_file, plot_pdf_file):
		"""
		Description:
		Function to create the parameter impact assessment PDF report for lsaBGC-Cluster.
		********************************************************************************************************************
		Parameters:
		- self: A Pan, GCF or BGC object.
		- workspace: The output/workspace with intermediate files (retained for backwards compatibility).
		- plot_input_1_file: The input table  
		- plot_annot_file: The input table used to create the heatmap visuals on the second page of the PDF report.
		- plot_input_2_file: The input table used to assess homogeneity in product/BGC-types per GCF clustering.
		- plot_overview_file: The input table with overview information on clustering using a variety of parameters.
		- plot_input_3_file: The input table with information on copy-count of GCFs per sample.
		- plot_pdf_file: The resulting path to the report PDF file.
		********************************************************************************************************************
		"""
		try:
			paramplots.create_param_influence_report(plot_input_1_file, plot_annot_file, plot_input_2_file,
													 plot_overview_file, plot_input_3_file, plot_pdf_file)
			if self.logObject:
				self.logObject.info('Created parameter impact report: %s' % plot_pdf_file)
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue creating the parameter impact report: %s' % plot_pdf_file)
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue creating the parameter impact report: %s' % plot_pdf_file) from e
