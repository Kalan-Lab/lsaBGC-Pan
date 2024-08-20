import os
import sys
import json
from Bio import SeqIO
import logging
import subprocess
import statistics
from operator import itemgetter
from collections import defaultdict
import traceback
import concurrent.futures
from ete3 import Tree
import concurrent.futures
import numpy as np
import gzip
import warnings
warnings.simplefilter('ignore')
import pkg_resources  # part of setuptools
import pandas as pd 

version = pkg_resources.require("lsaBGC-Pan")[0].version

valid_alleles = set(['A', 'C', 'G', 'T'])

def reformatOrthologInfo(ortholog_matrix_file, zol_results_dir, logObject):
	ortholog_listing_file = zol_results_dir + 'LocusTag_to_Ortholog_Relations.txt'
	try:
		outf = open(ortholog_listing_file, 'w')
		with open(ortholog_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: continue
				og = ls[0]
				for lts in ls[1:]:
					for lt in lts.split(', '):
						if lt.strip() != '':
							outf.write(lt + '\t' + og + '\n')
		outf.close()
		return (ortholog_listing_file)
	except:
		msg = 'Issue reformatting ortholog group matrix to table format for use as input for zol.'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def determineNonRepBGCs(sample, gcf, sample_gcf_bgcs, gcf_bgcs, bgc_pairwise_relations, edgy_bgcs, logObject):
	try:
		complete_bgcs = []
		for sample_bgc in sample_gcf_bgcs:
			if not sample_bgc in edgy_bgcs:
				complete_bgcs.append(sample_bgc)
		
		nonrep_bgcs = set([])
		if len(complete_bgcs) == 1:
			for sample_bgc in sample_gcf_bgcs:
				if sample_bgc != complete_bgcs[0]:
					nonrep_bgcs.add(sample_bgc)

		elif len(complete_bgcs) >= 2:
			sample_bgc_scores = []
			for sample_bgc in set(sample_gcf_bgcs).difference(edgy_bgcs):
				sample_bgc_relation_to_gcf_bgcs = 0.0
				for bgc in gcf_bgcs:
					sample_bgc_relation_to_gcf_bgcs += bgc_pairwise_relations[sample_bgc][bgc]
				sample_bgc_scores.append([sample_bgc, sample_bgc_relation_to_gcf_bgcs])

			rep_bgc = None
			for i, sbi in enumerate(sorted(sample_bgc_scores, key=itemgetter(1), reverse=True)): 
				if i == 0:
					rep_bgc = sbi[0]
				
			for sample_bgc in sample_gcf_bgcs:
				if sample_bgc != rep_bgc:
					nonrep_bgcs.add(sample_bgc)

		elif len(complete_bgcs) == 0:
			sample_bgc_scores = []
			for sample_bgc in sample_gcf_bgcs:
				sample_bgc_relation_to_gcf_bgcs = 0.0
				for bgc in gcf_bgcs:
					sample_bgc_relation_to_gcf_bgcs += bgc_pairwise_relations[sample_bgc][bgc]
				sample_bgc_scores.append([sample_bgc, sample_bgc_relation_to_gcf_bgcs])

			for i, sbi in enumerate(sorted(sample_bgc_scores, key=itemgetter(1), reverse=True)): 
				if i > 0:
					nonrep_bgcs.add(sbi[0])

		assert(len(nonrep_bgcs) > 0)
		return(nonrep_bgcs)
	except:
		msg = 'Issue selecting representative BGC for sample %s for GCF %s' % (sample, gcf)
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)


def runMIBiGMapper(detailed_BGC_listing_with_Pop_and_GCF_map_file, ortholog_matrix_file, mibig_dir, multi_thread, parallel_jobs_4thread, logObject):
	try:
		gcfs = set([])
		with open(detailed_BGC_listing_with_Pop_and_GCF_map_file) as odbl:
			for i, line in enumerate(odbl):
				if i == 0: continue
				line = line.strip()
				sample, population, method, genome_path, bgc_id, bgc_path, bgc_type, gcf_id, scaffold, start, end, bgc_length, dist_to_edge = line.split('\t')
				gcfs.add(gcf_id)

		lsabgc_map_cmds = []
		for gcf in gcfs:
			resdir = mibig_dir + gcf + '/'
			lsabgc_map_cmd = ['lsaBGC-MIBiGMapper', '-g', gcf, '-l', detailed_BGC_listing_with_Pop_and_GCF_map_file, 
					 		 '-m', ortholog_matrix_file, '-o', resdir, '-c', str(multi_thread), logObject]
			lsabgc_map_cmds.append(lsabgc_map_cmd)
		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_jobs_4thread) as executor:
				executor.map(multiProcess, lsabgc_map_cmds)
		except:
			msg = 'Issues with parallel running lsaBGC-MIBiGMapper commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

	except:
		msg = 'Issues running lsaBGC-MIBiGMapper'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def runSeeAndComprehenSeeIve(detailed_BGC_listing_with_Pop_and_GCF_map_file, species_tree, ortholog_matrix_file, see_dir, compsee_dir, threads, logObject):
	try:
		gcfs = set([])
		with open(detailed_BGC_listing_with_Pop_and_GCF_map_file) as odbl:
			for i, line in enumerate(odbl):
				if i == 0: continue
				line = line.strip()
				sample, population, method, genome_path, bgc_id, bgc_path, bgc_type, gcf_id, scaffold, start, end, bgc_length, dist_to_edge = line.split('\t')
				gcfs.add(gcf_id)

		lsabgc_see_and_csi_cmds = []
		for gcf in gcfs:
			gcf_see_resdir = see_dir + gcf + '/'
			gcf_csi_resdir = compsee_dir + gcf + '/'
			lsabgc_see_cmd = ['lsaBGC-See', '-g', gcf, '-l', detailed_BGC_listing_with_Pop_and_GCF_map_file, 
					          '-s', species_tree, '-m', ortholog_matrix_file, '-o', gcf_see_resdir, logObject]
			lsabgc_csi_cmd = ['lsaBGC-ComprehenSeeIve', '-g', gcf, '-l', detailed_BGC_listing_with_Pop_and_GCF_map_file, 
					          '-s', species_tree, '-m', ortholog_matrix_file, '-o', gcf_csi_resdir, logObject]
			lsabgc_see_and_csi_cmds.append(lsabgc_see_cmd)
			lsabgc_see_and_csi_cmds.append(lsabgc_csi_cmd)

		os.environ["OMP_NUM_THREADS"] = "1"
		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
				executor.map(multiProcess, lsabgc_see_and_csi_cmds)
		except:
			msg = 'Issues with parallel running of lsaBGC-See and lsaBGC-ComprehenSeeIve commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)
		os.environ["OMP_NUM_THREADS"] = "4"
	except:
		msg = 'Issues running lsaBGC-See and lsaBGC-ComprehenSeeIve'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def computeConservationOfOGWithinGCFContext(inputs):
	bgc_paths, og_listing_file, output_file, logObject = inputs
	
	try:
		lt_to_og = {}
		with open(og_listing_file) as oolf:
			for line in oolf:
				line = line.strip()
				ls = line.split('\t')
				lt_to_og[ls[0]] = ls[1]
		
		og_bgcs = defaultdict(set)
		for bgc in bgc_paths:
			with open(bgc) as obgbk:
				for rec in SeqIO.parse(obgbk, 'genbank'):
					for feat in rec.features:
						if feat.type != 'CDS': continue
						lt = feat.qualifiers.get('locus_tag')[0]
						og = lt_to_og[lt]
						og_bgcs[og].add(bgc)
		
		outf = open(output_file, 'w')
		for og in og_bgcs:
			og_bgc_prop = len(og_bgcs[og])/float(len(bgc_paths))
			outf.write(og + '\t' + str(og_bgc_prop) + '\n')
		outf.close()	
	except:
		msg = 'Issues computing conservation of orthogroups for complete instances of GCF %s' % output_file.split('/')[-1].split('.txt')[0]
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def runZol(detailed_BGC_listing_with_Pop_and_GCF_map_file, ortholog_listing_file, pairwise_relations, zol_comp_results_dir, zol_full_results_dir, cgc_results_dir, zol_parameters, zol_high_quality_preset, zol_edge_distance, zol_keep_multi_copy, threads, multi_thread, parallel_jobs_4thread, logObject):
	try:
		gcf_sample_bgcs = defaultdict(lambda: defaultdict(list))
		edgy_bgcs = set([])
		gcf_bgcs = defaultdict(list)
		with open(detailed_BGC_listing_with_Pop_and_GCF_map_file) as odbl:
			for i, line in enumerate(odbl):
				if i == 0: continue
				line = line.strip()
				sample, population, method, genome_path, bgc_id, bgc_path, bgc_type, gcf_id, scaffold, start, end, bgc_length, dist_to_edge = line.split('\t')
				gcf_sample_bgcs[gcf_id][sample].append(bgc_path)
				if float(dist_to_edge) <= zol_edge_distance:
					edgy_bgcs.add(bgc_path)
					gcf_bgcs[gcf_id].append(bgc_path)

		zol_cmds = []
		cgc_cmds = []
		complete_conservation_inputs = []
	
		zp = zol_parameters
		if zol_high_quality_preset:
			zp = '-qa -b -s'

		for gcf in gcf_sample_bgcs:
			gcf_bgcs_to_input = []
			for sample in gcf_sample_bgcs[gcf]:
				
				if len(gcf_sample_bgcs[gcf][sample]) > 1:
					if not zol_keep_multi_copy:
						samp_gcf_nonrep_bgcs = determineNonRepBGCs(sample, gcf, gcf_sample_bgcs[gcf][sample], gcf_bgcs[gcf], pairwise_relations, edgy_bgcs, logObject)
						samp_gcf_rep_bgcs = list(set(gcf_sample_bgcs[gcf][sample]).difference(samp_gcf_nonrep_bgcs))
						gcf_bgcs_to_input = gcf_bgcs_to_input + samp_gcf_rep_bgcs
					else:
						gcf_bgcs_to_input = gcf_bgcs_to_input + gcf_sample_bgcs[gcf][sample]
				else:
					gcf_bgcs_to_input.append(gcf_sample_bgcs[gcf][sample][0])

			gcf_full_bgcs_to_input = list(set(gcf_bgcs_to_input).difference(edgy_bgcs))
			if len(gcf_bgcs_to_input) > 0:
				zol_comp_cmd = ['zol', '-c', str(multi_thread), '-i', ' '.join(gcf_bgcs_to_input), zp, '-po', ortholog_listing_file, 
					            '-o', zol_comp_results_dir + gcf + '/', logObject]
				cgc_comp_cmd = ['cgc', '-i', zol_comp_results_dir + gcf + '/', '-t', 'conservation', 'conservation', 
					            'entropy', '-c', 'grey', 'white_to_black', 'light_to_dark_green', '-rh', '1.7', 
								'1', '1', '-sl', '-p', '-sc', '-ld', '0.04', '-lts', '3.0', '-b', '0.5', '-w', '8', 
								'-o', cgc_results_dir + gcf + '/', logObject]
				zol_cmds.append(zol_comp_cmd)
				cgc_cmds.append(cgc_comp_cmd)
			else:
				os.mkdir(zol_comp_results_dir + gcf)

			if len(gcf_full_bgcs_to_input) > 0:
				complete_conservation_input = [gcf_full_bgcs_to_input, ortholog_listing_file, zol_full_results_dir + gcf + '.txt', logObject]
				complete_conservation_inputs.append(complete_conservation_input)
			else:
				os.mkdir(zol_full_results_dir + gcf)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
				executor.map(computeConservationOfOGWithinGCFContext, complete_conservation_inputs)			
		except:
			msg = 'Issues with parallel computing of orthogroup conservations across complete instances.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_jobs_4thread) as executor:
				executor.map(multiProcess, zol_cmds)
		except:
			msg = 'Issues with parallel running of zol commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
				executor.map(multiProcess, cgc_cmds)
		except:
			msg = 'Issues with parallel running of cgc commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

	except:
		msg = 'Issues running zol or cgc'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)


def runCmdViaSubprocess(cmd, logObject=None, check_files=[], check_directories=[]):
	if logObject != None:
		logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		if logObject != None:
			logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		if logObject != None:
			logObject.error('Had an issue running: %s' % ' '.join(cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


def mapColorsToCategories(categories_set, colors_file, colors_mapping_file):
	try:
		colors = []
		with open(colors_file) as ocf:
			for line in ocf:
				line = line.strip()
				colors.append(line)
		
		assert(len(colors) >= len(categories_set))

		cmf_handle = open(colors_mapping_file, 'w')
		cmf_handle.write('category\tcolor\n')
		for i, cat in enumerate(sorted(categories_set)):
			cmf_handle.write(cat + '\t"' + colors[i] + '"\n')
		cmf_handle.close()
	
	except:
		msg = 'Issue mapping colors from file to categories, potentially because too few colors provided.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def pairwiseDistancesFromTree(tree_file, logObject):
	pairwise_distances = defaultdict(lambda: defaultdict(lambda: None))
	try:
		try:
			t = Tree(tree_file)
		except:
			return(pairwise_distances)
		
		leaves = set([])
		for node in t.traverse('postorder'):
			if node.is_leaf():
				leaves.add(node.name)

		for i, l1 in enumerate(sorted(leaves)):
			for j, l2 in enumerate(sorted(leaves)):
				if i == j:
					pairwise_distances[l1][l2] = 0.0
				elif i < j:
					dist = t.get_distance(l1, l2)
					pairwise_distances[l1][l2] = dist
					pairwise_distances[l2][l1] = dist
		
		return(pairwise_distances)
	except:	
		msg = 'Issue with calculating pairwise distances between leaves in the tree file: %s' % tree_file
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def generateColors(workspace, outfile, color_count, palette='Spectral', palette_color_count=12):
	try:
		assert(os.path.isdir(workspace))
		workspace = os.path.abspath(workspace) + '/'
	
		color_brew_script = workspace + 'brewColors.R'
		cbs_handle = open(color_brew_script, 'w')
		cbs_handle.write('library(RColorBrewer)\n')
		cbs_handle.write('mycolors <- colorRampPalette(brewer.pal(' + str(palette_color_count) + ', "' + palette + '"))(' + str(color_count) + ')\n')
		cbs_handle.write('write(mycolors, file="' + outfile + '")\n')		
		cbs_handle.close()

		rscript_cmd = ['Rscript', color_brew_script]
		runCmdViaSubprocess(rscript_cmd, check_files=[outfile])

	except:
		msg = 'Issues generating colors!'
		sys.stderr.write(msg + '\n')
		sys.exit(1)

def parseCDSCoord(str_gbk_loc):
	try:
		start = None
		end = None
		direction = None
		all_coords = []
		is_multi_part = False
		if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
			start = min([int(x.strip('>').strip('<')) for x in
						 str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
			end = max([int(x.strip('>').strip('<')) for x in
					   str(str_gbk_loc)[1:].split(']')[0].split(':')])
			direction = str(str_gbk_loc).split('(')[1].split(')')[0]
			all_coords.append([start, end, direction])
		elif 'order' in str(str_gbk_loc):
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		else:
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		return(all_coords, start, end, direction, is_multi_part)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())

def cleanUpSampleName(original_name):
	return original_name.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ',
																										'_').replace(
		':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(',
																													'').replace(
		')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')

def parseGECCOGBKForFunction(bgc_gbk, logObject):
	try:
		rec = SeqIO.read(bgc_gbk, 'genbank')
		product = 'unknown'
		try:
			product = rec.annotations['structured_comment']['GECCO-Data']['biosyn_class']
		except:
			try:
				product = rec.annotations['structured_comment']['GECCO-Data']['cluster_type']
			except:
				pass
		if product == 'Unknown':
			product = 'unknown'
		return(product)
	except:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		raise RuntimeError()

def parseAntiSMASHGBKForFunction(bgc_gbk, logObject, compress_multi=True):
	product = 'unknown'
	try:
		products = set([])
		with open(bgc_gbk) as obg:
			for rec in SeqIO.parse(obg, 'genbank'):
				for feat in rec.features:
					if feat.type == 'protocluster':
						try:
							products.add(feat.qualifiers.get('product')[0])
						except:
							pass
		if len(products) == 1:
			product = list(products)[0]
		elif len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products:
				product = 'NRPS'
		elif len(products) >= 2:
			if compress_multi:
				product = 'multi-type'
			else:
				product = ' | '.join(sorted(products))
	except:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		raise RuntimeError()
	return(product)

def parseOrthoFinderMatrix(orthofinder_matrix_file, relevant_gene_lts, all_primary=False):
	"""
	Function to parse and return information from OrthoFinderV2 de novo homolog group identification.

	:param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file
	:param relevant_gene_lts: set of all the relevant gene locus tag identifiers found in BGC Genbanks

	:return gene_to_hg: dictionary mapping gene locus tags to homolog group
	:return hg_genes: dictionary with set of gene locus tags for each homolog group
	:return hg_median_gene_counts: median copy count for each homolog group
	:return hg_multicopy_proportion: proportion of samples with homolog group which have multiple (paralogous) genes in the homolog group.
	"""
	gene_to_hg = {}
	hg_genes = defaultdict(set)
	hg_multicopy_proportion = defaultdict(lambda: 'NA')
	hg_median_gene_counts = defaultdict(lambda: 'NA')
	with open(orthofinder_matrix_file) as ofm:
		for i, line in enumerate(ofm):
			if i == 0: continue
			line = line.strip('\n')
			ls = line.split('\t')
			hg = ls[0]
			flag_in_bgc = False
			for sgs in ls[1:]:
				for g in sgs.split(', '):
					if g in relevant_gene_lts:
						flag_in_bgc = True
						gene_to_hg[g] = hg
						hg_genes[hg].add(g)
			if flag_in_bgc:
				gene_counts = []
				for sgs in ls[1:]:
					# critical for calculating homolog group stats, like median gene counts, multicopy proportion
					# use only genes from the original set of genomes used to conduct full orthofinder analysis.
					gene_counts.append(len([x for x in sgs.split(', ') if (len(x.split('_')[0]) == 3 or all_primary)]))

				hg_multicopy_proportion[hg] = float(sum([1 for x in gene_counts if x > 1])) / sum(
					[1 for x in gene_counts if x > 0])
				hg_median_gene_counts[hg] = statistics.median(gene_counts)

	return ([gene_to_hg, hg_genes, hg_median_gene_counts, hg_multicopy_proportion])

def run_cmd(cmd, logObject, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	"""
	Simple function to run a single command through subprocess with logging.
	"""
	logObject.info('Running the following command: %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr, executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))

def multiProcessNoLogWithTimeout(input):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list).
	"""
	input_cmd = input[:-1]
	TIMEOUT = input[-1] 
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash', timeout=TIMEOUT)
	except subprocess.TimeoutExpired:
		sys.stderr.write('Command timed out: ' + ' '.join(input_cmd))
		sys.stderr.write(traceback.format_exc())
		return

def multiProcessNoLog(input_cmd):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list).
	"""
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
	except Exception as e:
		sys.stderr.write(traceback.format_exc())

def multiProcess(input):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	"""
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except Exception as e:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))
		logObject.warning(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())

def addLocusTagsToGBKs(inputs):
	sample = inputs[0]
	locus_tag_prefix = inputs[1]
	resdir = inputs[2]
	genome_gbk = inputs[3]
	region_gbks = inputs[4:-2]
	keep_ids_flag = inputs[-2]
	logObject = inputs[-1]
	try:
		samp_resdir = resdir + sample + '/'
		setupReadyDirectory([samp_resdir])
		location_to_cds = {}
		cds_iter = 0
		updated_genome_gbk = samp_resdir + genome_gbk.split('/')[-1]	
		ugg_handle = open(updated_genome_gbk, 'w')
		with open(genome_gbk) as ogg:
			for rec in SeqIO.parse(ogg, 'genbank'):
				updated_rec = rec
				updated_features = []
				scaff = rec.id
				for feat in rec.features:
					if feat.type != 'CDS': continue
					cds_lt = None
					# first try setting locus tag to just protein_id
					if keep_ids_flag:
						try:
							cds_lt = feat.qualifiers.get('locus_tag')[0]
						except:
							try:
								cds_lt = feat.qualifiers.get('protein_id')[0]
							except:
								msg = 'Issue finding either locus_tag or protein_id associated with CDS feature at %s' % str(feat.location)
								sys.stderr.write(msg + '\n')
								sys.exit(1)

					all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feat.location))
					loctup = tuple([scaff, start])

					if cds_lt == None:
						cds_lt = locus_tag_prefix
						if (cds_iter+1) < 10:
							cds_lt += '_00000'+str(cds_iter+1)
						elif (cds_iter+1) < 100:
							cds_lt += '_0000'+str(cds_iter+1)
						elif (cds_iter+1) < 1000:
							cds_lt += '_000'+str(cds_iter+1)
						elif (cds_iter+1) < 10000:
							cds_lt += '_00'+str(cds_iter+1)
						elif (cds_iter+1) < 100000:
							cds_lt += '_0'+str(cds_iter+1)
						else:
							cds_lt += '_' + str(cds_iter+1)
						cds_iter += 1
					
					feat.qualifiers['locus_tag'] = cds_lt
					updated_features.append(feat)
					location_to_cds[loctup] = cds_lt

				updated_rec.features = updated_features
				SeqIO.write(updated_rec, ugg_handle, 'genbank')
		ugg_handle.close()

		for bgc_gbk in region_gbks:
			bgc_starts = []
			with open(bgc_gbk) as ogbf:
				for line in ogbf:
					line = line.strip()
					if line.startswith('Orig. start  ::'):
						bgc_starts.append(int(line.split()[-1].replace('>', '').replace('<', '')))			
			assert(len(bgc_starts) == 1)
			bgc_start = bgc_starts[0]

			updated_bgc_gbk = samp_resdir + bgc_gbk.split('/')[-1]	
			ubg_handle = open(updated_bgc_gbk, 'w')
			with open(bgc_gbk) as ogbf:
				for rec in SeqIO.parse(ogbf, 'genbank'):
					scaff = rec.id
					updated_rec = rec
					updated_features = []
					for feat in rec.features:
						if feat.type != 'CDS': 
							updated_features.append(feat)
							continue
						cds_lt = None
						
						all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feat.location))
											
						loctup = tuple([scaff, start+bgc_start])
						cds_lt = location_to_cds[loctup]
						feat.qualifiers['locus_tag'] = cds_lt
						updated_features.append(feat)

					updated_rec.features = updated_features
					SeqIO.write(updated_rec, ubg_handle, 'genbank')
			ubg_handle.close()

	except Exception as e:
		msg = "Problem processing one of the GenBank files for sample %s to add CDS locus tags." % sample
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.warning(msg)
		logObject.warning(traceback.format_exc())

def checkCDSHaveLocusTags(inputs):
	sample, gbk_type, gbk, outf, logObject = inputs
	try:
		cds_count = 0
		all_cds_have_lt = True
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feat in rec.features:
					if feat.type != 'CDS': continue
					locus_tag = None
					try:
						locus_tag = feat.qualifiers.get('locus_tag')[0]
					except:
						pass
					cds_count += 1
					if locus_tag == None:
						all_cds_have_lt = False

		outh = open(outf, 'w')
		outh.write(gbk_type + '\t' + str(sample) + '\t' + gbk + '\t' + str(cds_count) + '\t' + str(all_cds_have_lt) + '\n')
		outh.close()
	except Exception as e:
		logObject.warning("Problem processing GenBank file %s for sample %s" % (gbk, sample))
		logObject.warning(traceback.format_exc())

def findAntiSMASHBGCInFullGenbank(inputs):
	try:
		full_gbk, bgc_gbk, outf = inputs

		bgc_starts = []
		with open(bgc_gbk) as ogbf:
			for line in ogbf:
				line = line.strip()
				if line.startswith('Orig. start'):
					bgc_starts.append(int(line.split()[-1].replace('>', '').replace('<', '')))			

		bgc_scaffold = None
		bgc_length = 0
		with open(bgc_gbk) as obg:
			for rec in SeqIO.parse(obg, 'genbank'):
				bgc_scaffold = rec.id
				bgc_length = len(str(rec.seq))

		full_scaff_length = None
		with open(full_gbk) as ofg:
			for rec in SeqIO.parse(ofg, 'genbank'):
				if rec.id == bgc_scaffold:
					full_scaff_length = len(str(rec.seq))

		assert(full_scaff_length != None)				
		outf_handle = open(outf, 'w')
		for start_coord in bgc_starts:
			end_coord = start_coord + bgc_length - 1
			outf_handle.write(bgc_scaffold + '\t' + str(start_coord) + '\t' + str(end_coord) + '\t' + str(min([start_coord, full_scaff_length-end_coord])) + '\n')
		outf_handle.close()

	except Exception as e:
		raise RuntimeWarning(traceback.format_exc())

def processAntiSMAHSBGCtoGenomeMappingResults(process_data, logObject):
	try:
		loc_lists = defaultdict(list)
		for	pd in process_data:
			sample = pd[0]
			bgc_file = pd[2]
			result_file = pd[3]
			with open(result_file) as orf:
				for line in orf:
					line = line.strip()
					scaff, start_coord, end_coord, dist_to_edge = line.split('\t')
					bgc_length = int(end_coord) - int(start_coord) + 1
					loc_list = [sample, scaff, int(start_coord), int(end_coord), bgc_length, int(dist_to_edge), 'antismash']
					loc_lists[bgc_file].append(loc_list)

		bgc_locations = []
		for bgc in loc_lists:
			if len(loc_lists[bgc]) != 1: 
				sys.stderr.write('Warning: no location or more than 2 potential locations for BGC region %s. Skipping incorporation.' % bgc)
				continue
			else:
				bgc_locations.append([bgc] + loc_lists[bgc][0])
		return bgc_locations		
	
	except Exception as e:
		msg = 'Issue processing location mapping information.'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def extractProteinsFromGenBank(inputs):
	input_genbank, output_proteome, logObject = inputs
	try:
		op_handle = open(output_proteome, 'w')
		with open(input_genbank) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						prot_lt = feature.qualifiers.get('locus_tag')[0]
						prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
						op_handle.write('>' + prot_lt + '\n' + prot_seq + '\n')
		op_handle.close()
	except Exception as e:
		logObject.error("Issues with parsing out protein sequences from GenBank file %s." % input_genbank) 
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def determineAsofName(asof_index):
	asof_index_str = str(asof_index)
	asof_name = None
	if len(asof_index_str) == 1:
		asof_name = '000000' + asof_index_str
	elif len(asof_index_str) == 2:
		asof_name = '00000' + asof_index_str
	elif len(asof_index_str) == 3:
		asof_name = '0000' + asof_index_str
	elif len(asof_index_str) == 4:
		asof_name = '000' + asof_index_str
	elif len(asof_index_str) == 5:
		asof_name = '00' + asof_index_str
	elif len(asof_index_str) == 6:
		asof_name = '0' + asof_index_str
	else:
		asoof_name = asof_index_str
	assert(asof_name != None)
	return(asof_name)

def runOrthoFinder2Full(prot_directory, orthofinder_outdir, run_msa, logObject, threads=1):
	result_file = orthofinder_outdir + 'Final_Orthogroups.tsv'
	try:
		orthofinder_cmd = ['orthofinder', '-f', prot_directory, '-t', str(threads)]

		if run_msa:
			orthofinder_cmd += ['-M', 'msa']

		logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
		subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran OrthoFinder!')
		tmp_orthofinder_dir = os.path.abspath(
			[prot_directory + 'OrthoFinder/' + f for f in os.listdir(prot_directory + 'OrthoFinder/') if
			 f.startswith('Results')][0]) + '/'

		os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		n0_file = orthofinder_outdir + 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
		putative_xenologs_dir = orthofinder_outdir + 'Putative_Xenologs/'
		phylo_misplaced_genes_dir = orthofinder_outdir + 'Phylogenetically_Misplaced_Genes/'
		
		result_handle = open(result_file, 'w')

		gene_to_hog = {}
		gene_to_og = {}
		with open(n0_file) as n0f:
			for i, line in enumerate(n0f):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes = ls[3:]
					result_handle.write('Orthogroup\t' + '\t'.join(genomes) + '\n')
				else:
					hog = ls[0].split('N0.')[1]
					og = ls[1]
					for j, gs in enumerate(ls[3:]):
						genome = genomes[j]
						for gene in gs.split(', '):
							gene = gene.strip()
							gene_to_hog[gene] = hog
							gene_to_og[gene] = og

		genome_misplaced_genes = defaultdict(set)
		fg_to_genome = {}
		for f in os.listdir(phylo_misplaced_genes_dir):
			genome = f.split('.txt')[0]
			with open(phylo_misplaced_genes_dir + f) as opmgdf:
				for i, line in enumerate(opmgdf):
					line = line.strip()
					genome_misplaced_genes[genome].add(line)
					fg_to_genome[line] = genome
					
		close_hogs = defaultdict(lambda: defaultdict(int))
		for f in os.listdir(putative_xenologs_dir):
			genome = f.split('.tsv')[0]
			with open(putative_xenologs_dir + f) as opxf:
				for i, line in enumerate(opxf):
					if i == 0: continue
					line = line.strip()
					og, focal_genes, other_genes = line.split('\t')
					for fg in focal_genes.split(', '):
						if not fg in genome_misplaced_genes[genome]: continue
						for otg in other_genes.split(', '):
							if otg in gene_to_hog:
								close_hogs[fg][gene_to_hog[otg]] += 1

		hog_missing_to_add = defaultdict(lambda: defaultdict(set))
		for fg in close_hogs:
			max_value = max(close_hogs[fg].values())
			if max_value == 0: continue
			top_hits = 0 
			top_hog = None
			for hog in close_hogs[fg]:
				if close_hogs[fg][hog] == max_value:
					top_hog = hog
					top_hits += 1
			if top_hits == 1:
				hog_missing_to_add[top_hog][genome].add(fg)

		genomes = []
		genome_genes_accounted = defaultdict(set)
		result_handle = open(result_file, 'w')
		with open(n0_file) as n0f:
			for i, line in enumerate(n0f):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes = ls[3:]
					result_handle.write('Orthogroup\t' + '\t'.join(genomes) + '\n')
				else:
					hog = ls[0].split('N0.')[1]
					og = ls[1]
					printlist = []
					for j, gs in enumerate(ls[3:]):
						genome = genomes[j]
						gss = set(gs.split(', '))
						gss_with_missing = gss.union(hog_missing_to_add[hog][genome]) 
						printlist.append(', '.join(sorted(gss_with_missing)))
						for gene in gss_with_missing:
							genome_genes_accounted[genome].add(gene)
					result_handle.write(hog + '\t' + '\t'.join(printlist) + '\n')
		
		genomes_sf = []
		with open(singletons_file) as osf:
			for i, line in enumerate(osf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes_sf = ls[1:]
				else:
					for j, gs in enumerate(ls[1:]):
						genome = genomes_sf[j]
						for gene in gs.split(', '):
							genome_genes_accounted[genome].add(gene)
					result_handle.write(line + '\n')

		genomes_mf = []
		asof_index = 0
		with open(main_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes_mf = ls[1:]
				else:
					og = ls[0]
					printlist = [og]
					value_count = 0
					for j, gs in enumerate(ls[1:]):
						genome = genomes_mf[j]
						for gene in gs.split(', '):
							if not gene in genome_genes_accounted[genome]:
								asof_name = 'ASOF' + determineAsofName(asof_index)	
								printlist = [asof_name]
								for gen in genomes_mf:
									if gen == genome:
										printlist.append(gene)
									else:
										printlist.append('')
								result_handle.write('\t'.join(printlist) + '\n')
								asof_index += 1
		result_handle.close()

		assert (genomes == genomes_mf and genomes == genomes_sf)
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return result_file


def runOrthoFinder2FullFungal(prot_directory, orthofinder_outdir, run_msa, logObject, threads=1):
	result_file = orthofinder_outdir + 'Final_Orthogroups.tsv'
	try:
		orthofinder_cmd = ['orthofinder', '-f', prot_directory, '-t', str(threads)]

		if run_msa:
			orthofinder_cmd += ['-M', 'msa']

		logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
		subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran OrthoFinder!')
		tmp_orthofinder_dir = os.path.abspath(
			[prot_directory + 'OrthoFinder/' + f for f in os.listdir(prot_directory + 'OrthoFinder/') if
			 f.startswith('Results')][0]) + '/'

		os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		n0_file = orthofinder_outdir + 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
		
		result_handle = open(result_file, 'w')

		hog_ids = set([])
		hog_prots = set([])
		genomes = []
		with open(n0_file) as n0f:
			for i, line in enumerate(n0f):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes = ls[3:]
					result_handle.write('Orthogroup\t' + '\t'.join(genomes) + '\n')
				else:
					hog = ls[0].split('N0.')[1]
					hog_ids.add(int(hog.split('HOG')[1]))
					result_handle.write('\t'.join([hog] + ls[3:]) + '\n')
					for lts in ls[3:]:
						for lt in lts.split(', '):
							if lt.strip() != '':
								hog_prots.add(lt)

		hog_id_iter = max(hog_ids) + 1
		for prot_file in os.listdir(prot_directory):
			if not os.path.isfile(prot_directory + prot_file) or not prot_file.endswith('.faa'): continue
			genome_id = '.faa'.join(prot_file.split('.faa')[:-1])
			with open(prot_directory + prot_file) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					if not rec.id in hog_prots:
						hog_id_iter_str = 'HOG'
						if (hog_id_iter) < 10:
							hog_id_iter_str += '00000'+str(hog_id_iter)
						elif (hog_id_iter) < 100:
							hog_id_iter_str += '0000'+str(hog_id_iter)
						elif (hog_id_iter) < 1000:
							hog_id_iter_str += '000'+str(hog_id_iter)
						elif (hog_id_iter) < 10000:
							hog_id_iter_str += str(hog_id_iter)
						elif (hog_id_iter) < 100000:
							hog_id_iter_str += str(hog_id_iter)
						else:
							hog_id_iter_str += str(hog_id_iter)
						printlist = [hog_id_iter_str]
						hog_id_iter += 1
						for genome_iter in genomes:
							if genome_id == genome_iter:
								printlist.append(rec.id)
							else:
								printlist.append('')
						result_handle.write('\t'.join(printlist) + '\n')	
		result_handle.close()
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return result_file


def runOrthoFinder2(prot_directory, orthofinder_outdir, logObject, threads=1):
	result_file = orthofinder_outdir + 'Final_Orthogroups.tsv'
	try:
		orthofinder_cmd = ['orthofinder', '-f', prot_directory, '-t', str(threads), '-og']
		logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
		subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran OrthoFinder!')
		tmp_orthofinder_dir = os.path.abspath(
			[prot_directory + 'OrthoFinder/' + f for f in os.listdir(prot_directory + 'OrthoFinder/') if
			 f.startswith('Results')][0]) + '/'
		os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		result_handle = open(result_file, 'w')
		with open(main_file) as omf:
			for line in omf:
				result_handle.write(line)
		with open(singletons_file) as osf:
			for i, line in enumerate(osf):
				if i == 0: continue
				result_handle.write(line)
		result_handle.close()
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return result_file

def runPanaroo(detailed_BGC_listing_file, panaroo_input_dir, results_directory, logObject, threads=1, panaroo_options='--clean-mode moderate --remove-invalid-genes'):
	result_file = results_directory + 'Final_Orthogroups.tsv'
	try:
		sample_genomes = {}
		with open(detailed_BGC_listing_file) as odlf:
			for i, line in enumerate(odlf):
				if i == 0: continue
				line = line.strip()
				sample, _, genome_gbk = line.split('\t')[:3]
				sample_genomes[sample] = genome_gbk

		reformat_cmds = []
		panaroo_inputs = []
		for sample in sample_genomes:
			inf = sample_genomes[sample]
			outf = panaroo_input_dir + '.'.join(inf.split('/')[-1].split('.')[:-1]) + '.gff'
			reformat_cmd = ['genbankToProkkaGFF.py', '-i', inf, '-o', outf, logObject]
			reformat_cmds.append(reformat_cmd)
			panaroo_inputs.append(outf)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
				executor.map(multiProcess, reformat_cmds)
		except:
			msg = 'Issues creating Prokka-like GFF files for Panaroo.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.exit(1)

		panaroo_cmd = ['panaroo', '-t', str(threads), panaroo_options, '-i', ' '.join(panaroo_inputs), '-o', results_directory]

		try:
			logObject.info('Running the following command: %s' % ' '.join(panaroo_cmd))
			subprocess.call(' '.join(panaroo_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran Panaroo!')
		except:
			msg = 'Issue running Panaroo with the following command: %s' % ' '.join(panaroo_cmd)
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.exit(1)

		main_ortho_file = results_directory + 'gene_presence_absence.csv'
		
		result_handle = open(result_file, 'w')
		clustered_lts = set([])
		cds_iter = 1
		genomes = []
		with open(main_ortho_file) as omof:
			for i, line in enumerate(omof):
				line = line.strip('\n')
				ls = line.split(',')
				if i == 0: 
					result_handle.write('\t'.join(['OrthoGroup'] + ls[3:]) + '\n')
					genomes = ls[3:]
				else:
					og_id = 'OG' + determineAsofName(cds_iter)
					for lts in ls[3:]:
						for lt in lts.split(';'):
							lt = lt.strip()
							if lt != '': clustered_lts.add(lt)
					result_handle.write('\t'.join([og_id] + [', '.join(x.split(';')) for x in ls[3:]]) + '\n')
					cds_iter += 1
		
		for sample in sample_genomes:
			gbk = sample_genomes[sample]
			with open(gbk) as ogbk:
				for rec in SeqIO.parse(ogbk, 'genbank'):
					for feat in rec.features:
						if not feat.type == 'CDS': continue
						lt = feat.qualifiers.get('locus_tag')[0]
						if not lt in clustered_lts:
							og_id = 'OG' + determineAsofName(cds_iter)
							cds_iter += 1
							printlist = [og_id]
							for genome in genomes:
								if genome == sample:
									printlist.append(lt)
								else:
									printlist.append('')
							result_handle.write('\t'.join(printlist) + '\n')
		result_handle.close()

		assert (os.path.isfile(result_file))
	except Exception as e:
		msg = "Problem with determining orthogroups using Panaroo."
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	return result_file

def selectFinalResultsAndCleanUp(outdir, fin_outdir, logObject):

	""" 
	TODO: Refactor for lsaBGC-PAN
	"""
	try:
		delete_set = set(['BLASTing_of_Ortholog_Groups', 'OrthoFinder2_Results', 'KOfam_Annotations',
						  'Prodigal_Gene_Calling_Additional', 'Predicted_Proteomes_Initial',
						  'Prodigal_Gene_Calling', 'Genomic_Genbanks_Initial'])
		for fd in os.listdir(outdir):
			if os.path.isfile(fd): continue
			subdir = outdir + fd + '/'
			if fd in delete_set:
				os.system('rm -rf %s' % subdir)
		if os.path.isdir(outdir + 'BGCs_Retagged_and_Updated'):
			os.system('rm -rf %s' % outdir + 'BGCs_Retagged')
		if os.path.isdir(outdir + 'lsaBGC_AutoExpansion_Results/'):
			os.system('ln -s %s %s' % (
			outdir + 'lsaBGC_AutoExpansion_Results/Updated_GCF_Listings/', fin_outdir + 'Expanded_GCF_Listings'))
			os.system('ln -s %s %s' % (
			outdir + 'lsaBGC_AutoExpansion_Results/Orthogroups.expanded.tsv', fin_outdir + 'Expanded_Orthogroups.tsv'))
			os.system('ln -s %s %s' % (outdir + 'lsaBGC_AutoExpansion_Results/Sample_Annotation_Files.txt',
									   fin_outdir + 'Expanded_Sample_Annotation_Files.txt'))
			if os.path.isfile(outdir + 'Intermediate_Files/GToTree_output.tre'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/GToTree_output.tre', fin_outdir))
			if os.path.isfile(outdir + 'Intermediate_Files/GToTree_Expected_Similarities.txt'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/GToTree_Expected_Similarities.txt', fin_outdir))
			if os.path.isfile(outdir + 'Intermediate_Files/Samples_in_GToTree_Tree.txt'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/Samples_in_GToTree_Tree.txt', fin_outdir))
		else:
			os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/*', fin_outdir))
			if os.path.isdir(outdir + 'BiG_SCAPE_Results_Reformatted/'):
				os.system('ln -s %s %s' % (
				outdir + 'BiG_SCAPE_Results_Reformatted/GCF_Listings/', fin_outdir + 'GCF_Listings'))
				os.system('ln -s %s %s' % (
				outdir + 'BiG_SCAPE_Results_Reformatted/GCF_Details.txt', fin_outdir + 'GCF_Details.txt'))
			elif os.path.isdir(outdir + 'lsaBGC_Cluster_Results/'):
				os.system(
					'ln -s %s %s' % (outdir + 'lsaBGC_Cluster_Results/GCF_Listings/', fin_outdir + 'GCF_Listings'))
				os.system('ln -s %s %s' % (
				outdir + 'lsaBGC_Cluster_Results/GCF_Details_Expanded_Singletons.txt', fin_outdir + 'GCF_Details.txt'))
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def setupReadyDirectory(directories):
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				os.system('rm -rf %s' % d)
			os.system('mkdir %s' % d)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def p_adjust_bh(p):
	"""
	Benjamini-Hochberg p-value correction for multiple hypothesis testing.
	"""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

def is_newick(newick):
	"""
	Function to validate if Newick phylogeny file is correctly formatted.
	"""
	try:
		t = Tree(newick)
		return True
	except:
		return False


def is_fastq(fastq):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		with open(fastq) as of:
			SeqIO.parse(of, 'fastq')
		return True
	except:
		return False


def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				SeqIO.parse(ogf, 'fasta')
		else:
			with open(fasta) as of:
				SeqIO.parse(of, 'fasta')
		return True
	except:
		return False


def is_genbank(gbk):
	"""
	Function to check in Genbank file is correctly formatted.
	"""
	try:
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				SeqIO.parse(ogf, 'genbank')
		else:
			with open(gbk) as of:
				SeqIO.parse(of, 'genbank')
		return True
	except:
		return False

def parseVersionFromSetupPy():
	"""
	Parses version from setup.py program.
	"""
	return(str(version))

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""
	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger


def closeLoggerObject(logObject):
	"""
	Function which closes/terminates loggerObject.
	:param logObject: logging logger object to close
	"""
	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)


def logParameters(parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to std.stderr
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')


def logParametersToFile(parameter_file, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()


def logParametersToObject(logObject, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))


def determineSeqSimProteinAlignment(protein_alignment_file, use_only_core=True):
	protein_sequences = {}
	with open(protein_alignment_file) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			protein_sequences[rec.id] = str(rec.seq).upper()

	pair_seq_matching = defaultdict(lambda: defaultdict(lambda: 0.0))
	for i, g1 in enumerate(sorted(protein_sequences)):
		s1 = g1.split('|')[0]
		g1s = protein_sequences[g1]
		for j, g2 in enumerate(sorted(protein_sequences)):
			if i >= j: continue
			s2 = g2.split('|')[0]
			if s1 == s2: continue
			g2s = protein_sequences[g2]
			tot_comp_pos = 0
			match_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a != '-' or g2a != '-':
					if not use_only_core or (use_only_core and g1a != '-' and g2a != '-'):
						tot_comp_pos += 1
						if g1a == g2a:
							match_pos += 1
			general_matching_percentage = 0.0
			if tot_comp_pos > 0:
				general_matching_percentage = float(match_pos) / float(tot_comp_pos)
			if pair_seq_matching[s1][s2] < general_matching_percentage and pair_seq_matching[s2][
				s1] < general_matching_percentage:
				pair_seq_matching[s1][s2] = general_matching_percentage
				pair_seq_matching[s2][s1] = general_matching_percentage

	return pair_seq_matching

def castToNumeric(x):
	"""
	Description:
	This function attempts to cast a variable into a float. A special exception is whether "< 3 segregating sites!" is
	the value of the variable, which will simply be retained as a string.
	********************************************************************************************************************
	Parameters:
	- x: Input variable.
	********************************************************************************************************************
	Returns:
	- A float casting of the variable's value if numeric or "nan" if not.
	********************************************************************************************************************
	"""
	try:
		if x == '< 3 segregating sites!':
			return(x)
		else:
			x = float(x)
			return (x)
	except:
		return float('nan')

def loadTableInPandaDataFrame(input_file, numeric_columns):
	"""
	Description:
	This function formats reads a TSV file and stores it as a pandas dataframe.
	********************************************************************************************************************
	Parameters:
	- input_file: The input TSV file, with first row corresponding to the header.
	- numeric_columns: Set of column names which should have numeric data.
	********************************************************************************************************************
	Returns:
	- panda_df: A pandas DataFrame object reprsentation of the input TSV file.
	********************************************************************************************************************
	"""
	import pandas as pd
	panda_df = None
	try:
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			key = ls[0]
			cast_vals = ls[1:]
			if key in numeric_columns:
				cast_vals = []
				for val in ls[1:]:
					cast_vals.append(castToNumeric(val))
			panda_dict[key] = cast_vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return panda_df

def createFinalSpreadsheets(detailed_BGC_listing_with_Pop_and_GCF_map_file, zol_results_dir, zol_high_qual_flag,
							     mibig_dir, recon_result_file, recon_og_result_file, recon_pop_color_file,
								 sociate_result_file, final_spreadsheet_xlsx, scratch_dir, logObject):
	try:
		# generate Excel spreadsheet
		writer = pd.ExcelWriter(final_spreadsheet_xlsx, engine='xlsxwriter')
		workbook = writer.book

		er_sheet = workbook.add_worksheet('Explanation of Results')
		er_sheet.write(0, 0, 'An explanation on the interpretation of the different sheets in this spreadsheet can be found online at:')
		er_sheet.write(1, 0, 'https://github.com/Kalan-Lab/lsaBGC-Pan/wiki/8.-explanation-of-final-spreadsheet-and-visual-reports')

		# specify different types of cell formatting 

		warn_format = workbook.add_format({'bg_color': '#bf241f', 'bold': True, 'font_color': '#FFFFFF'})
		na_format = workbook.add_format({'font_color': '#a6a6a6', 'bg_color': '#FFFFFF', 'italic': True})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})

		# create detailed BGC info sheet
		
		bo_numeric_columns = set(['start', 'end', 'bgc_length', 'dist_to_edge'])
		bo_data = loadTableInPandaDataFrame(detailed_BGC_listing_with_Pop_and_GCF_map_file, bo_numeric_columns)
		bo_data.to_excel(writer, sheet_name='BGC Overview', index=False, na_rep="NA")
		bo_sheet =  writer.sheets['BGC Overview']
		bo_sheet.conditional_format('A1:BA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# create zol results spreadsheet

		zol_combined_tsv_file = scratch_dir + 'zol_results.tsv'

		zol_sheet_header = ['GCF ID', 'Ortholog Group (OG) ID', 'OG is Single Copy?', 'Proportion of Total Gene Cluster Instances with OG', 
					        'Proprtion of Complete Gene Cluster Instances with OG', 'OG Median Length (bp)', 'OG Consensus Order', 
							'OG Consensus Direction', 'Tajima\'s D', 'Proportion of Filtered Codon Alignment is Segregating Sites', 
							'Entropy', 'Upstream Region Entropy', 'Median Beta-RD-gc', 'Max Beta-RD-gc', 
							'Proportion of sites which are highly ambiguous in codon alignment', 
							'Proportion of sites which are highly ambiguous in trimmed codon alignment', 'Median GC', 'Median GC Skew']

		if zol_high_qual_flag:
			zol_sheet_header += ['GARD Partitions Based on Recombination Breakpoints',
			           'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
				       'Average delta(Beta, Alpha) by FUBAR across sites',
				       'Proportion of Sites Under Selection which are Positive'] 
					
		zol_sheet_header += ['KO Annotation (E-value)', 'PGAP Annotation (E-value)', 'PaperBLAST Annotation (E-value)', 'CARD Annotation (E-value)',
							'IS Finder (E-value)', 'MIBiG Annotation (E-value)', 'VOG Annotation (E-value)', 'VFDB Annotation (E-value)', 
							'Pfam Domains', 'CDS Locus Tags', 'OG Consensus Sequence']
		
		# ^ basically added two columns (GCF id and complete instances conservation) and took away one (custom db annotation)

		zol_full_dir = zol_results_dir + 'Comprehensive/'
		gcf_comp_cons_dir = zol_results_dir + 'Complete_Instances/'	

		comp_cons = defaultdict(lambda: defaultdict(lambda: 'NA'))
		for f in os.listdir(gcf_comp_cons_dir):
			gcf_comp_cons_file = gcf_comp_cons_dir + f
			gcf = f.split('.txt')[0]
			if not os.path.isfile(gcf_comp_cons_file): continue			
			with open(gcf_comp_cons_file) as ogrf:
				for i, line in enumerate(ogrf):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					og = ls[0]
					cons = ls[1]
					comp_cons[gcf][og] = cons

		zctf_handle = open(zol_combined_tsv_file, 'w')
		zctf_handle.write('\t'.join(zol_sheet_header) + '\n')
		num_rows = 1
		for gcf in os.listdir(zol_full_dir):
			gcf_result_file = zol_full_dir + gcf + '/Final_Results/Consolidated_Report.tsv'
			if not os.path.isfile(gcf_result_file): continue
			with open(gcf_result_file) as ogrf:
				for i, line in enumerate(ogrf):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					row = [gcf, ls[0], ls[1], ls[2], comp_cons[gcf][ls[0]]] + ls[3:16] + ls[17:]
					if zol_high_qual_flag:
						row = [gcf, ls[0], ls[1], ls[2], comp_cons[gcf][ls[0]]] + ls[3:20] + ls[21:]
					zctf_handle.write('\t'.join(row) + '\n')
					num_rows += 1
		zctf_handle.close()

		zr_numeric_columns = set(['Proportion of Total Gene Cluster Instances with OG', 'Proprtion of Complete Gene Cluster Instances with OG', 
							      'OG Median Length (bp)', 'OG Consensus Order', 'Tajima\'s D', 'GARD Partitions Based on Recombination Breakpoints',
			           			  'Number of Sites Identified as Under Positive or Negative Selection by FUBAR', 'Average delta(Beta, Alpha) by FUBAR across sites',
				     		      'Proportion of Sites Under Selection which are Positive', 'Proportion of Filtered Codon Alignment is Segregating Sites',
								  'Entropy', 'Upstream Region Entropy', 'Median Beta-RD-gc', 'Max Beta-RD-gc', 'Proportion of sites which are highly ambiguous in codon alignment', 
								  'Proportion of sites which are highly ambiguous in trimmed codon alignment', 'Median GC', 'Median GC Skew'])
		
		zr_data = loadTableInPandaDataFrame(zol_combined_tsv_file, zr_numeric_columns)
		zr_data.to_excel(writer, sheet_name='zol Results', index=False, na_rep="NA")
		zr_sheet =  writer.sheets['zol Results']

		zr_sheet.conditional_format('C2:C' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"False"', 'format': warn_format})
		zr_sheet.conditional_format('A2:CA' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"NA"', 'format': na_format})
		zr_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# prop gene-clusters with hg
		zr_sheet.conditional_format('D2:D' + str(num_rows), {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
		zr_sheet.conditional_format('E2:E' + str(num_rows), {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		# gene-lengths
		zr_sheet.conditional_format('F2:F' + str(num_rows), {'type': '2_color_scale', 'min_color': "#a3dee3", 'max_color': "#1ebcc9", "min_value": 100, "max_value": 2500, 'min_type': 'num', 'max_type': 'num'})
		
		# taj-d
		zr_sheet.conditional_format('I2:I' + str(num_rows),
										{'type': '3_color_scale', 'min_color': "#f7a09c", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										'max_color': "#87cefa", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0})

		# prop seg sites
		zr_sheet.conditional_format('J2:J' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#eab3f2", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#a37ba8", "min_value": 0.0, "max_value": 1.0})

		# entropy
		zr_sheet.conditional_format('K2:K' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

		# upstream region entropy
		zr_sheet.conditional_format('L2:L' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

		# median beta-rd gc
		zr_sheet.conditional_format('M2:M' + str(num_rows),
										{'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})
		# max beta-rd gc
		zr_sheet.conditional_format('N2:N' + str(num_rows),
										{'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})

		# ambiguity full ca
		zr_sheet.conditional_format('O2:O' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

		# ambiguity trim ca
		zr_sheet.conditional_format('P2:P' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

		# GC
		zr_sheet.conditional_format('Q2:Q' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#abffb7", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#43bf55", "min_value": 0.0, "max_value": 1.0})

		# GC Skew
		zr_sheet.conditional_format('R2:R' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#c7afb4", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#965663", "min_value": -2.0, "max_value": 2.0})


		# create MIBiG mapping spreadsheet
		mibig_json_tar_url = 'https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz'
		mibig_json_tar_dir = scratch_dir + 'mibig_json_3.1/'
		mibig_json_tar_file = scratch_dir + 'mibig_json_3.1.tar.gz'

		if not os.path.isdir(mibig_json_tar_dir):
			if os.path.isfile(mibig_json_tar_file):
				os.remove(mibig_json_tar_file)
			wget_cmd = ['axel', mibig_json_tar_url, '-o', mibig_json_tar_file]

			try:
				subprocess.call(' '.join(wget_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
				assert (os.path.isfile(mibig_json_tar_file))
			except:
				sys.stderr.write('Had an issue running: %s\n' % ' '.join(wget_cmd))
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)

			uncompress_cmd = ['tar', '-zxvf', mibig_json_tar_file, '-C', scratch_dir]
			try:
				os.mkdir(mibig_json_tar_dir)
				subprocess.call(' '.join(uncompress_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
				assert (os.path.isdir(mibig_json_tar_dir))
			except:
				sys.stderr.write('Had an issue running: %s\n' % ' '.join(uncompress_cmd))
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)

		mibig_bgc_compounds = defaultdict(list)
		if os.path.isdir(mibig_json_tar_dir):
			for f in os.listdir(mibig_json_tar_dir):
				if not f.endswith('.json'): continue
				bgc = f.split('.')[0]
				json_bgc_file = mibig_json_tar_dir + f
				with open(json_bgc_file) as json_data:
					data = json.load(json_data)
					for comp in data['cluster']['compounds']:
						for key in comp:
							if key == 'compound':
								mibig_bgc_compounds[bgc].append(comp[key])


		mb_combined_tsv_file = scratch_dir + 'mibig_map_results.tsv'
		mb_header = ['GCF ID', 'MIBiG BGC ID', 'GCF OG ID', 'MIBiG Protein Matching', 'MIBiG Compound(s)']
		mb_outf = open(mb_combined_tsv_file, 'w')
		mb_outf.write('\t'.join(mb_header) + '\n')
		for gcf in os.listdir(mibig_dir):
			map_file = mibig_dir + gcf + '/GCF_to_MIBiG_Relations.txt'
			if not os.path.isfile(map_file): continue
			with open(map_file) as omf:
				for i, line in enumerate(omf):
					if i == 0: continue
					line = line.strip('\n')
					ls = line.split('\t')
					mb_outf.write('\t'.join(ls + [' | '.join(mibig_bgc_compounds[ls[1]])]) + '\n')
		mb_outf.close()

		mb_data = loadTableInPandaDataFrame(mb_combined_tsv_file, set([]))
		mb_data.to_excel(writer, sheet_name='lsaBGC-MIBiGMapper Results', index=False, na_rep="NA")
		mb_sheet =  writer.sheets['lsaBGC-MIBiGMapper Results']
		mb_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# create reconcile spreadsheet

		rs_numeric_columns = set([])
		rs_non_numeric_columns = set(['orthogroup', 'found in non-BGC context', 'GCFs'])

		with open(recon_result_file) as orrf:
			for i, line in enumerate(orrf):
				if i == 0:
					line = line.strip('\n')
					colnames = set(line.split('\t'))
					rs_numeric_columns = colnames.difference(rs_non_numeric_columns)

		rs_data = loadTableInPandaDataFrame(recon_result_file, rs_numeric_columns)
		rs_data.to_excel(writer, sheet_name='lsaBGC-Reconcile Results', index=False, na_rep="NA")
		rs_sheet =  writer.sheets['lsaBGC-Reconcile Results']
		rs_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# create OG by sample/population spreadsheet
		# TODO: add colors for populations

		og_numeric_columns = set([])
		og_non_numeric_columns = set(['sample'])

		num_rows = 1
		with open(recon_result_file) as orrf:
			for i, line in enumerate(orrf):
				if i == 0:
					line = line.strip('\n')
					colnames = set(line.split('\t'))
					og_numeric_columns = colnames.difference(og_non_numeric_columns)
					num_rows += 1
		
		og_data = loadTableInPandaDataFrame(recon_og_result_file, og_numeric_columns)
		og_data.to_excel(writer, sheet_name='BGC OG by Sample Matrix', index=False, na_rep="NA")
		og_sheet =  writer.sheets['BGC OG by Sample Matrix']
		og_sheet.conditional_format('A1:CA2', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})
		og_sheet.conditional_format('B3:CA' + str(num_rows),
										{'type': '2_color_scale', 'min_color': "#f8f5fc", 'min_type': 'num', 'max_type': 'num',
										'max_color': "#5a739c", "min_value": 0.0, "max_value": 2.0})		

		# create lsaBGC-sociate spreadsheet
		so_numeric_columns = set(['allele frequency', 'frequency in samples with focal GCF', 'frequency in samples without focal GCF','pvalue', 
					  'phylgoenetically corrected pvalue', 'beta', 'beta-std-err', 'variant_h2'])
		so_data = loadTableInPandaDataFrame(sociate_result_file, so_numeric_columns)
		so_data.to_excel(writer, sheet_name='lsaBGC-Sociate Results', index=False, na_rep="NA")
		so_sheet =  writer.sheets['lsaBGC-Sociate Results']
		so_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# close workbook 
		workbook.close()
	except:
		msg = 'Difficulties creating final multi-sheet spreadsheet XLSX file!'
		logObject.error(msg)
		sys.stderr.write(msg +  '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
