import os
import sys
import math
import json
from Bio import SeqIO
import logging
import subprocess
import shutil
import statistics
from operator import itemgetter
from collections import defaultdict, OrderedDict
from Bio.SeqRecord import SeqRecord
import traceback
import concurrent.futures
from ete3 import Tree
from tqdm import tqdm
import numpy as np
import gzip
import warnings
warnings.simplefilter('ignore')
import importlib.metadata
import pandas as pd 

try:
	package_name = "lsaBGC-Pan"
	package_version = str(importlib.metadata.version(package_name))
except importlib.metadata.PackageNotFoundError:
	package_version = "NA"
	
valid_alleles = set(['A', 'C', 'G', 'T'])

# characters in paths which would break commands run through the shell (whitespace and shell special characters)
SHELL_UNSAFE_CHARACTERS = set([' ', '\t', '\n', ';', '&', '|', '<', '>', '(', ')', '$', '`', "'", '"', '\\', '*', '?', '!', '[', ']', '{', '}'])

def reformatOrthologInfo(ortholog_matrix_file, zol_results_dir, logObject):
	"""
	Description:
	Function to reformat sample by ortholog group matrix to a format that can be fed into zol directly as precomputed.
	********************************************************************************************************************
	Parameters:
	- ortholog_matrix_file: Path to sample by ortholog matrix file.
	- zol_results_dir: The workspace where to write the resulting re-formatted ortholog info file.
	- logObject: A python logging object
	********************************************************************************************************************
	Returns:
	- ortholog_listing_file: Path to the reformatted ortholog information file in table format.
	********************************************************************************************************************
	"""
	ortholog_listing_file = zol_results_dir + 'LocusTag_to_Ortholog_Relations.txt'
	try:
		with open(ortholog_listing_file, 'w') as outf, open(ortholog_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: continue
				og = ls[0]
				for lts in ls[1:]:
					for lt in lts.split(', '):
						if lt.strip() != '':
							outf.write(lt + '\t' + og + '\n')
		return (ortholog_listing_file)
	except Exception as e:
		msg = 'Issue reformatting ortholog group matrix to table format for use as input for zol.'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def determineNonRepBGCs(sample, gcf, sample_gcf_bgcs, gcf_bgcs, bgc_pairwise_relations, edgy_bgcs, logObject):
	"""
	Description:
	Determine which BGC instances to retain for zol analysis for sample with multiple copies of a single GCF.
	********************************************************************************************************************
	Parameters:
	- sample: sample name
	- gcf: GCF identifier
	- sample_gcf_bgcs: The set of all BGC instances for the sample-GCF pairing.
	- gcf_bgcs: The set of all BGC instances for the GCF.
	- bgc_pairwise_relations: A dictionary with sequence differences between pairs of BGCs.
	- edgy_bgcs: The set of BGCs which are near scaffold/contig edges.
	- logObject: A python logging object
	********************************************************************************************************************
	Returns:
	- nonrep_bgcs: The set of non-representative BGCs for the sample-GCF pairing determined.
	********************************************************************************************************************
	"""
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
	except Exception as e:
		msg = 'Issue selecting representative BGC for sample %s for GCF %s' % (sample, gcf)
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)


def runMIBiGMapper(detailed_BGC_listing_with_Pop_and_GCF_map_file, ortholog_matrix_file, mibig_dir, multi_thread, parallel_jobs_4thread, logObject):
	"""
	Description:
	Void wrapper function to run lsaBGC-MIBiGMapper for all GCFs in parallel. 
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_with_pop_and_GCF_map_file: The detailed BGC listing file produced in lsaBGC-Pan with much of 
	                                                  the information for sample mapping, GCF mapping, and context of 
													  each BGC.
	- ortholog_matrix_file: The sample by ortholog matrix.
	- mibig_dir: The results directory where to write lsaBGC-MIBiGMapper results to. 
	- multi_threads: The number of threads per job.
	- parallel_jobs_4thread: The number of parallel jobs to run whereby each is assumed to be using 4 threads.
	- logObject: A python logging object
	********************************************************************************************************************
	"""
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
				submitted_jobs = [executor.submit(multiProcess, arg) for arg in lsabgc_map_cmds]
				for future in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Running lsaBGC-MIBiGMapper"):
					future.result()
		except Exception as e:
			msg = 'Issues with parallel running lsaBGC-MIBiGMapper commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

	except Exception as e:
		msg = 'Issues running lsaBGC-MIBiGMapper'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def runSeeAndComprehenSeeIve(detailed_BGC_listing_with_Pop_and_GCF_map_file, species_tree, ortholog_matrix_file, see_dir, compsee_dir, threads, logObject):
	"""
	Description:
	Void wrapper function to run lsaBGC-See and lsaBGC-ComprehenSeeIve for all GCFs in parallel. 
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_with_pop_and_GCF_map_file: The detailed BGC listing file produced in lsaBGC-Pan with much of 
	                                                  the information for sample mapping, GCF mapping, and context of 
													  each BGC.
	- species_tree: A species tree in Newick format.
	- ortholog_matrix_file: The sample by ortholog matrix.
	- see_dir: The results directory where to write lsaBGC-See results to. 
	- compsee_dir: The results directory where to write lsaBGC-ComprehenSeeIve results to.
	- threads: The number of parallel jobs to run whereby each is assumed to be using 1 thread.
	- logObject: A python logging object
	********************************************************************************************************************
	"""
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
				submitted_jobs = [executor.submit(multiProcess, arg) for arg in lsabgc_see_and_csi_cmds]
				for future in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Running lsaBGC-See and lsaBGC-ComprehenSeeIve"):
					future.result()
		except Exception as e:
			msg = 'Issues with parallel running of lsaBGC-See and lsaBGC-ComprehenSeeIve commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)
		os.environ["OMP_NUM_THREADS"] = "4"
	except Exception as e:
		msg = 'Issues running lsaBGC-See and lsaBGC-ComprehenSeeIve'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def computeConservationOfOGWithinGCFContext(inputs):
	"""
	Description:
	Void function which computes the conservation of an orthogroup across BGC instances within a GCF. 
	********************************************************************************************************************
	Parameters:
	- inputs: A list for four items:
		- bgc_paths: A list of BGC paths belonging to the focal GCF.
		- og_listing_file: The locus tag to orthogroup mapping table.
		- output_file: The output file which maps orthogroups to their conservation.
		- logObject: A python logging object.
	********************************************************************************************************************
	"""
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
		
		with open(output_file, 'w') as outf:
			for og in og_bgcs:
				og_bgc_prop = len(og_bgcs[og])/float(len(bgc_paths))
				outf.write(og + '\t' + str(og_bgc_prop) + '\n')
	except Exception as e:
		msg = 'Issues computing conservation of orthogroups for complete instances of GCF %s' % output_file.split('/')[-1].split('.txt')[0]
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def runZol(detailed_BGC_listing_with_Pop_and_GCF_map_file, ortholog_listing_file, pairwise_relations, zol_comp_results_dir, zol_full_results_dir, cgc_results_dir, zol_parameters, zol_high_quality_preset, zol_edge_distance, zol_keep_multi_copy, threads, multi_thread, parallel_jobs_4thread, logObject):
	"""
	Description:
	Void wrapper function to run zol for all GCFs in parallel. 
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_with_pop_and_GCF_map_file: The detailed BGC listing file produced in lsaBGC-Pan with much of 
	                                                  the information for sample mapping, GCF mapping, and context of 
													  each BGC.
	- ortholog_listing_file: The locus tag to orthogroup mapping table.
	- pairwise_relations: A dictionary with pairwise relations between pairs of BGCs
	- zol_comp_results_dir: The results directory where to write zol comprehensive results to. 
	- zol_full_results_dir: The results directory where to write orthogroup conservations to based on full/complete BGC
	                        instances belonging to a GCF.
	- cgc_results_dir: The results directory where to write cgc results to.
	- zol_parameters: The parameters with which to run zol analysis. 
	- zol_high_quality_preset: Flag for whether to run a more full/high-quality zol analysis.
	- zol_edge_distance: The distance from the edge of a contig/scaffold for a BGC instance to be considered partial.
	- zol_keep_multi_copy: Flag for whether to retain multi-copy instances of a GCF per sample.
	- threads: The number of parallel jobs to run assuming each job uses 1 thread.
	- multi_thread: The number of threads to use per multi-threaded job - typically hardcoded to 4.
	- parallel_jobs_4thread: The number of parallel jobs to run assuming each job uses 4 threads.
	- logObject: A python logging object.
	********************************************************************************************************************
	"""
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
				gcf_bgcs[gcf_id].append(bgc_path)
				if float(dist_to_edge) <= zol_edge_distance:
					edgy_bgcs.add(bgc_path)

		zol_cmds = []
		cgc_cmds = []
		complete_conservation_inputs = []
	
		zp = zol_parameters
		if zol_high_quality_preset:
			# zol >= 1.7.0 replaced -qa/--quality-align with -am/--aligner-method (MUSCLE5 -align is the accurate option)
			zp = '-am muscle -s -b'

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
				submitted_jobs = [executor.submit(computeConservationOfOGWithinGCFContext, arg) for arg in complete_conservation_inputs]
				for future in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Computing conservation of orthogroups"):
					future.result()
		except Exception as e:
			msg = 'Issues with parallel computing of orthogroup conservations across complete instances.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_jobs_4thread) as executor:
				submitted_jobs = [executor.submit(multiProcess, arg) for arg in zol_cmds]
				for _ in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Running zol"):
					pass
		except Exception as e:
			msg = 'Issues with parallel running of zol commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

		try:
			with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
				submitted_jobs = [executor.submit(multiProcess, arg) for arg in cgc_cmds]
				for future in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Running cgc"):
					future.result()
		except Exception as e:
			msg = 'Issues with parallel running of cgc commands.'
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			sys.stderr.write(traceback.format_exc())
			logObject.error(traceback.format_exc())
			sys.exit(1)

	except Exception as e:
		msg = 'Issues running zol or cgc'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)


def prioritizeEnvironmentBinPath():
	"""
	Description:
	Places the bin/ directory of the python interpreter running lsaBGC-Pan (i.e. the conda environment) at the front
	of PATH, so that dependencies started via "#!/usr/bin/env python3" shebangs (e.g. OrthoFinder, which bofasa calls)
	use the environment's python and its packages rather than another python (e.g. Homebrew's) found earlier on PATH.
	"""
	env_bin = os.path.dirname(os.path.abspath(sys.executable))
	path_dirs = [d for d in os.environ.get('PATH', '').split(os.pathsep) if d]
	if len(path_dirs) == 0 or os.path.abspath(path_dirs[0]) != env_bin:
		path_dirs = [env_bin] + [d for d in path_dirs if os.path.abspath(d) != env_bin]
		os.environ['PATH'] = os.pathsep.join(path_dirs)


def runCmdViaSubprocess(cmd, logObject=None, check_files=[], check_directories=[], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL):
	"""
	Description:
	Void function to run a command using python's subprocess. 
	********************************************************************************************************************
	Parameters:
	- cmd: A list corresponding to a command.
	- logObject [optional]: A python logging object.
	- check_files [optional]: A list of files to simply check are made after running the command.
	- check_directories [optional]: A list of directories to simply check are made after running the command.
	********************************************************************************************************************
	"""
	if logObject != None:
		logObject.info('Running %s' % ' '.join(cmd))
	opened_handles = []
	log_files = [x for x in [stdout, stderr] if isinstance(x, str)]
	try:
		if isinstance(stdout, str):
			stdout = open(stdout, 'w')
			opened_handles.append(stdout)
		if isinstance(stderr, str):
			stderr = open(stderr, 'w')
			opened_handles.append(stderr)

		return_code = subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr,
						executable='/bin/bash')
		if return_code != 0:
			# if expected outputs are specified, rely on them to determine success (some tools can
			# exit non-zero after producing their results), otherwise treat the non-zero exit as a failure.
			if len(check_files) == 0 and len(check_directories) == 0:
				raise RuntimeError('Command exited with non-zero return code %d' % return_code)
			msg = 'Warning: command exited with non-zero return code %d: %s' % (return_code, ' '.join(cmd))
			sys.stderr.write(msg + '\n')
			if logObject != None:
				logObject.warning(msg)
		for cf in check_files:
			if not os.path.isfile(cf):
				raise RuntimeError('Expected output file was not created: %s' % cf)
		for cd in check_directories:
			if not os.path.isdir(cd):
				raise RuntimeError('Expected output directory was not created: %s' % cd)
		if logObject != None:
			logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except Exception as e:
		msg = 'Had an issue running: %s' % ' '.join(cmd)
		if log_files:
			msg += '\nCheck the tool\'s log files for details: %s' % ', '.join(log_files)
		if logObject != None:
			logObject.error(msg)
			logObject.error(traceback.format_exc())
		raise RuntimeError(msg) from e
	finally:
		for handle in opened_handles:
			handle.close()


def mapColorsToCategories(categories_set, colors_file, colors_mapping_file):
	"""
	Description:
	Void function to map colors to categories for various plotting purposes.
	********************************************************************************************************************
	Parameters:
	- categories_set: Set of categories to have colors assigned to.
	- colors_file: A file where each line corresponds to a hex code for a color.
	- colors_mapping_file: The result file where the first column has a category and the second the mapped color.
	********************************************************************************************************************
	"""
	try:
		colors = []
		with open(colors_file) as ocf:
			for line in ocf:
				line = line.strip()
				colors.append(line)
		
		assert(len(colors) >= len(categories_set))

		with open(colors_mapping_file, 'w') as cmf_handle:
			cmf_handle.write('category\tcolor\n')
			for i, cat in enumerate(sorted(categories_set)):
				cmf_handle.write(cat + '\t"' + colors[i] + '"\n')
	
	except Exception as e:
		msg = 'Issue mapping colors from file to categories, potentially because too few colors provided.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def pairwiseDistancesFromTree(tree_file, logObject):
	"""
	Description:
	Function to compute the branch distance between leaves in a newick tree.
	********************************************************************************************************************
	Parameters:
	- tree_file: The input tree file in Newick format.
	- logObject: a python logging object
	********************************************************************************************************************
	Returns:
	-  pairwise_distances: A dictionary showing branch distances between pairs of leaves.
	********************************************************************************************************************
	"""
	pairwise_distances = defaultdict(lambda: defaultdict(lambda: None))
	try:
		try:
			t = Tree(tree_file)
		except Exception as e:
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
	except Exception as e:	
		msg = 'Issue with calculating pairwise distances between leaves in the tree file: %s' % tree_file
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

# ColorBrewer palettes (https://colorbrewer2.org) used for coloring in lsaBGC-Pan (largest version of each).
BREWER_PALETTES = {
	'Spectral': ['#9E0142', '#D53E4F', '#F46D43', '#FDAE61', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5',
				 '#3288BD', '#5E4FA2'],
	'Set1': ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'],
	'Set2': ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3'],
	'Pastel1': ['#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4', '#FED9A6', '#FFFFCC', '#E5D8BD', '#FDDAEC', '#F2F2F2']}
# qualitative palettes share colors across their n-color versions, the diverging Spectral palette only has a fixed
# version here (the 11-color one), which is what was always requested in lsaBGC-Pan.
QUALITATIVE_PALETTES = set(['Set1', 'Set2', 'Pastel1'])


def brewerColorRamp(palette, color_count, palette_color_count=None):
	"""
	Description:
	Function to get colors interpolated along a ColorBrewer palette - equivalent to the R code:
	colorRampPalette(brewer.pal(palette_color_count, palette))(color_count)
	********************************************************************************************************************
	Parameters:
	- palette: The ColorBrewer palette name (Spectral, Set1, Set2 or Pastel1).
	- color_count: The count of desired colors.
	- palette_color_count [Optional]: The number of palette colors to interpolate between [Default is color_count].
	********************************************************************************************************************
	Returns:
	- A list of hex color codes.
	********************************************************************************************************************
	"""
	color_count = int(color_count)
	if color_count <= 0:
		return []
	base = BREWER_PALETTES[palette]
	if palette_color_count is None:
		palette_color_count = color_count
	if palette in QUALITATIVE_PALETTES:
		# brewer.pal returns at least 3 and at most the maximum number of colors for a palette
		base = base[:max(3, min(int(palette_color_count), len(base)))]
	# mirror R's arithmetic exactly: colorRamp() linearly interpolates channels on a 0-1 scale at positions
	# seq(0, 1, length.out=n) using approx(), the result is multiplied by 255 and rgb(..., maxColorValue=255) truncates.
	knots_y = [[int(c[i:i + 2], 16) / 255.0 for c in base] for i in (1, 3, 5)]
	m = len(base)
	knots_x = [i * (1.0 / (m - 1)) for i in range(m)] if m > 1 else [0.0]

	def _approx(v, y):
		if m == 1:
			return y[0]
		i, j = 0, m - 1
		if v <= knots_x[0]:
			return y[0]
		if v >= knots_x[-1]:
			return y[-1]
		while i < j - 1:
			ij = (i + j) // 2
			if v < knots_x[ij]:
				j = ij
			else:
				i = ij
		if v == knots_x[j]:
			return y[j]
		if v == knots_x[i]:
			return y[i]
		return y[i] + (y[j] - y[i]) * ((v - knots_x[i]) / (knots_x[j] - knots_x[i]))

	if color_count == 1:
		positions = [0.0]
	else:
		by = 1.0 / (color_count - 1)
		positions = [k * by for k in range(color_count)]
	colors = []
	for v in positions:
		channels = [min(max(_approx(v, knots_y[ch]), 0.0), 1.0) * 255 for ch in range(3)]
		colors.append('#%02X%02X%02X' % tuple(int(x) for x in channels))
	return colors


def generateColors(workspace, outfile, color_count, palette='Spectral', palette_color_count=12):
	"""
	Description:
	Void function to generate an output file with colors in hex code for various plotting purposes.
	********************************************************************************************************************
	Parameters:
	- workspace: The workspace where to write temporary files (retained for backwards compatibility; unused).
	- outfile: The output file where hex codes of colors will be written to.
	- color_count: The count of desired colors.
	- palette [Optional]: The ColorBrewer color palette to use.
	- palette_color_count [Optional]: The color palette's default color count.
	********************************************************************************************************************
	"""
	try:
		colors = brewerColorRamp(palette, color_count, palette_color_count=palette_color_count)
		with open(outfile, 'w') as of:
			if colors:
				of.write('\n'.join(colors) + '\n')
	except Exception as e:
		msg = 'Issues generating colors!'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def parseCDSCoord(location_string, record_length=None):
	"""
	Description:
	Function to process a location string from a GenBank file (as produced by str(feature.location) in Biopython).
	Multi-part locations (e.g. spliced genes) are supported, including genes spanning the origin of circular records,
	e.g. join{[900:1000](+), [0:60](+)}. 
	********************************************************************************************************************
	Parameters:
	- location_string: The location string to process.
	- record_length (optional): Length of the record the feature is on. Used for features spanning the origin, whose
	                            end coordinate is then reported as continuing past the end of the record (e.g. 1060 for
	                            the example above), such that start < end and end - start + 1 is the feature's length.
	********************************************************************************************************************
	Returns:
	- A list of items: [all_coords, start, end, direction, is_multi_part], where all_coords is a list of [start, end, 
	  direction] for each part ordered along the genome (for features spanning the origin, starting with the part 
	  before the origin), start/end are 1-based inclusive coordinates.
	********************************************************************************************************************
	"""
	try:
		all_starts = []
		all_ends = []
		all_coords = []
		direction = None
		if "order" in str(location_string):
			all_directions = [] # type: ignore
			for exon_coord in location_string[6:-1].split(", "):
				start = (min([int(x.strip(">").strip("<")) for x in exon_coord[1:].split("]")[0].split(":")])+ 1)
				end = max([int(x.strip(">").strip("<")) for x in exon_coord[1:].split("]")[0].split(":")])
				direction = exon_coord.split("(")[1].split(")")[0]
				all_starts.append(start)
				all_ends.append(end)
				all_directions.append(direction)
				all_coords.append([start, end, direction])
			start = min(all_starts)
			end = max(all_ends)
			assert len(set(all_directions)) == 1
			direction = all_directions[0]
		elif "join" in location_string:
			all_directions = []
			for exon_coord in location_string[5:-1].split(", "):
				start = (min([int(x.strip(">").strip("<")) for x in exon_coord[1:].split("]")[0].split(":")])+ 1)
				end = max([int(x.strip(">").strip("<")) for x in exon_coord[1:].split("]")[0].split(":")])
				direction = exon_coord.split("(")[1].split(")")[0]
				all_starts.append(start)
				all_ends.append(end)
				all_directions.append(direction)
				all_coords.append([start, end, direction])
			start = min(all_starts)
			end = max(all_ends)
			assert len(set(all_directions)) == 1
			direction = all_directions[0]
		elif "{" not in location_string:
			start = (min([int(x.strip(">").strip("<")) for x in location_string[1:].split("]")[0].split(":")])+ 1)
			end = max([int(x.strip(">").strip("<")) for x in location_string[1:].split("]")[0].split(":")])
			direction = location_string.split("(")[1].split(")")[0]
			all_starts.append(start)
			all_ends.append(end)
			all_coords.append([start, end, direction])
		else:
			msg = 'Error: There appears to be a location operator that is neither "join" nor "order". This is currently not supported in lsaBGC-Pan.'
			sys.stderr.write(msg + "\n")
			raise RuntimeError(f"Error processing location string {location_string}")

		is_multi_part = False
		if len(all_coords) > 1:
			is_multi_part = True
			# Biopython lists parts in biological order: ascending for + strand features and descending for - strand
			# features. A violation of this order indicates the feature spans the origin of a circular record.
			listed_starts = [c[0] for c in all_coords]
			if direction == '-':
				spans_origin = any(listed_starts[i] > listed_starts[i - 1] for i in range(1, len(listed_starts)))
			else:
				spans_origin = any(listed_starts[i] < listed_starts[i - 1] for i in range(1, len(listed_starts)))
			ascending_coords = sorted(all_coords, key=itemgetter(0))
			if spans_origin:
				# order parts along the genome across the origin: the largest gap between consecutive parts (in 
				# ascending order) is the region not covered by the feature, so parts after it come first.
				gaps = [ascending_coords[i + 1][0] - ascending_coords[i][1] for i in range(len(ascending_coords) - 1)]
				k = gaps.index(max(gaps)) + 1
				all_coords = ascending_coords[k:] + ascending_coords[:k]
				start = all_coords[0][0]
				if record_length != None:
					end = all_coords[-1][1] + record_length
				else:
					end = max([c[1] for c in ascending_coords[k:]])
			else:
				all_coords = ascending_coords
		return [all_coords, start, end, direction, is_multi_part] # type: ignore
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		raise RuntimeError(
			f"Error processing location string {location_string}"
		) from e

def cleanUpSampleName(original_name):
	"""
	Description:
	Function to clean up a string corresponding to some name to make it easier to write files referencing.
	********************************************************************************************************************
	Parameters:
	- name: A string corresponding to some name.
	********************************************************************************************************************
	Returns:
	-  cleaned_up_name: A cleaned up string where troublesome characters for having as part of filenames are replaced 
	                    with less problematic characters.
	********************************************************************************************************************
	"""
	return original_name.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ',
																										'_').replace(
		':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(',
																													'').replace(
		')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')

def parseGECCOGBKForFunction(bgc_gbk, logObject):
	"""
	Description:
	Function to parse GECCO BGC GenBank for product function
	********************************************************************************************************************
	Parameters:
	- bgc_gbk: The GECCO BGC GenBank file path.
	- logObject: a python logging object
	********************************************************************************************************************
	Returns:
	-  product: The predicted product/type of BGC by GECCO.
	********************************************************************************************************************
	"""	
	try:
		rec = SeqIO.read(bgc_gbk, 'genbank')
		product = 'unknown'
		try:
			product = rec.annotations['structured_comment']['GECCO-Data']['biosyn_class']
		except Exception as e:
			try:
				product = rec.annotations['structured_comment']['GECCO-Data']['cluster_type']
			except Exception as e:
				pass
		if product == 'Unknown':
			product = 'unknown'
		return(product)
	except Exception as e:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		logObject.error(traceback.format_exc())
		raise RuntimeError('Issues parsing BGC Genbank %s' % bgc_gbk) from e

def parseAntiSMASHGBKForFunction(bgc_gbk, logObject, compress_multi=True):
	"""
	Description:
	Function to parse antiSMASH BGC GenBank for product function.
	********************************************************************************************************************
	Parameters:
	- bgc_gbk: The antiSMASH BGC GenBank file path.
	- logObject: a python logging object.
	********************************************************************************************************************
	Returns:
	-  product: The predicted product(s)/type(s) of BGC by antiSMASH.
	********************************************************************************************************************
	"""	
	product = 'unknown'
	try:
		products = set([])
		with open(bgc_gbk) as obg:
			for rec in SeqIO.parse(obg, 'genbank'):
				for feat in rec.features:
					if feat.type == 'protocluster':
						try:
							products.add(feat.qualifiers.get('product')[0])
						except Exception as e:
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
	except Exception as e:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		logObject.error(traceback.format_exc())
		raise RuntimeError('Issues parsing BGC Genbank %s' % bgc_gbk) from e
	return(product)

def parseOrthoFinderMatrix(orthofinder_matrix_file, relevant_gene_lts, all_primary=True):
	"""
	Description:	
	Function to parse and return information from OrthoFinderV2 de novo homolog group identification.
	********************************************************************************************************************
	Parameters:
	- orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file
	- relevant_gene_lts: set of all the relevant gene locus tag identifiers found in BGC Genbanks
	********************************************************************************************************************
	Return:
	- gene_to_hg: dictionary mapping gene locus tags to homolog group
	- return hg_genes: dictionary with set of gene locus tags for each homolog group
	- return hg_median_gene_counts: median copy count for each homolog group
	- return hg_multicopy_proportion: proportion of samples with homolog group which have multiple (paralogous) genes in the homolog group.
	********************************************************************************************************************
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

				try:
					hg_multicopy_proportion[hg] = float(sum([1 for x in gene_counts if x > 1])) / sum(
						[1 for x in gene_counts if x > 0])
					hg_median_gene_counts[hg] = statistics.median(gene_counts)
				except Exception as e:
					print(line)
	return ([gene_to_hg, hg_genes, hg_median_gene_counts, hg_multicopy_proportion])

def run_cmd(cmd, logObject, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	"""
	Description:
	Void simple function to run a single command through subprocess with logging.
	********************************************************************************************************************
	Parameters:
	- cmd: A list corresponding to a command.
	- logObject: A python logging object.
	********************************************************************************************************************
	"""

	logObject.info('Running the following command: %s' % ' '.join(cmd))
	try:
		return_code = subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr, executable='/bin/bash')
		if return_code != 0:
			raise RuntimeError('Command exited with non-zero return code %d' % return_code)
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd)) from e

def multiProcessNoLogWithTimeout(input):
	"""
	Description:
	Void, generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list). Similar to more general function, only this one times out after a while.
	********************************************************************************************************************
	Parameters:
	- input:
		- input_cmd: Up to last item the input should correspond to the command to run in subprocess. 
		- TIMEOUT: The last item in the list should correspond to the timeout in seconds.
	********************************************************************************************************************
	"""
	input_cmd = input[:-1]
	TIMEOUT = input[-1] 
	try:
		return_code = subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash', timeout=TIMEOUT)
		if return_code != 0:
			sys.stderr.write('Warning: command exited with non-zero return code %d: %s\n' % (return_code, ' '.join(input_cmd)))
	except subprocess.TimeoutExpired:
		sys.stderr.write('Warning: command timed out after %s seconds: %s\n' % (str(TIMEOUT), ' '.join(input_cmd)))
		return

def multiProcessNoLog(input_cmd):
	"""
	Description:
	Void, generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list). Same as the generalizable function, only this one doesn't take a log object
	as part of the input.
	********************************************************************************************************************
	Parameters:
	- input_cmd: A list corresponding to an input command.
	"""
	try:
		# Check if command contains output redirection
		cmd_str = ' '.join(input_cmd)
		if '>' in cmd_str:
			# Use shell=True to handle redirection properly
			return_code = subprocess.call(cmd_str, shell=True, executable='/bin/bash')
		else:
			# For commands without redirection, use the original method
			return_code = subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
		if return_code != 0:
			sys.stderr.write('Warning: command exited with non-zero return code %d: %s\n' % (return_code, cmd_str))
	except Exception as e:
		sys.stderr.write(traceback.format_exc())

def multiProcess(input):
	"""
	Description:
	Void, generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	********************************************************************************************************************
	Parameters: 
	- input: 
		- input_cmd: Up to last item the input should correspond to the command to run in subprocess.
		- logObject: The last item in the list should correspond to a python logging object.
	********************************************************************************************************************
	"""
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		return_code = subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if return_code != 0:
			msg = 'Warning: command exited with non-zero return code %d: %s' % (return_code, ' '.join(input_cmd))
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')
		else:
			logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except Exception as e:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))
		logObject.warning(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())

def addLocusTagsToGBKs(inputs):
	"""
	Description:
	Void function add locus tags to a sample's genome and BGC GenBank files.
	********************************************************************************************************************
	Parameters:
	- inputs:
		- sample: The sample name.
		- locus_tag_prefix: The locus tag prefix to use.
		- resdir: The results directory to write the updated GenBank files to.
		- genome_gbk: The genome GenBank file to process.
		- region_gbks: The BGC region GenBank files to process.
		- keep_ids_flag: A flag variable on whether to attempt to keep the original locus tags.
		- logObject: a python logging object.
	********************************************************************************************************************
	"""
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
		with open(updated_genome_gbk, 'w') as ugg_handle, open(genome_gbk) as ogg:
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
						except Exception as e:
							#try:
							#	cds_lt = feat.qualifiers.get('protein_id')[0]
							#except Exception as e:
							msg = f'Issue finding a locus_tag associated with CDS feature at { str(feat.location)} of scaffold {scaff}'
							sys.stderr.write(msg + '\n')
							logObject.info(msg)
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

		for bgc_gbk in region_gbks:
			bgc_starts = []
			with open(bgc_gbk) as ogbf:
				for line in ogbf:
					line = line.strip()
					if line.startswith('Orig. start') and '::' in line:
						bgc_starts.append(int(line.split()[-1].replace('>', '').replace('<', '')))			
			assert(len(bgc_starts) == 1)
			bgc_start = bgc_starts[0]

			updated_bgc_gbk = samp_resdir + bgc_gbk.split('/')[-1]	
			with open(updated_bgc_gbk, 'w') as ubg_handle, open(bgc_gbk) as ogbf:
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

	except Exception as e:
		msg = "Problem processing one of the GenBank files for sample %s to add CDS locus tags." % sample
		logObject.warning(msg)
		logObject.warning(traceback.format_exc())
		sys.stderr.write(msg)
		raise RuntimeError(msg) from e

def checkCDSHaveLocusTags(inputs):
	"""
	Description:
	Void function to check if all CDS features in a GenBank file have locus_tag qualifiers.
	********************************************************************************************************************
	Parameters:
	- inputs:
		- sample: The sample name.
		- gbk_type: The type of the GenBank file.
		- gbk: The input GenBank file.
		- outf: A file to write stats on CDS count and whether they all had locus tags. 
		- logObject: a python logging object.
	******************************************************************************************************************** 
	"""
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
					except Exception as e:
						pass

					cds_count += 1
					if locus_tag == None:
						all_cds_have_lt = False

		with open(outf, 'w') as outh:
			outh.write(gbk_type + '\t' + str(sample) + '\t' + gbk + '\t' + str(cds_count) + '\t' + str(all_cds_have_lt) + '\n')
	except Exception as e:
		logObject.warning("Problem processing GenBank file %s for sample %s" % (gbk, sample))
		logObject.warning(traceback.format_exc())
		raise RuntimeWarning("Problem processing GenBank file %s for sample %s" % (gbk, sample)) from e

def findAntiSMASHBGCInFullGenbank(inputs):
	"""
	Description:
	Void function to parse antiSMASH BGC GenBank from antiSMASH to determine the location/coordinates of it in the 
	genome.
	********************************************************************************************************************
	Parameters:
	- inputs:
		- full_gbk: The full genome GenBank file.
		- bgc_gbk: The antiSMASH BGC GenBank file.
		- outf: The output file where to write information on BGC location to.
	********************************************************************************************************************
	"""
	try:
		sample, genome_gbk, sample_bgc_info_tuples = inputs
		full_genome_index = SeqIO.index(genome_gbk, 'genbank')

		for bgc_gbk, outf in sample_bgc_info_tuples:
			bgc_starts = []
			with open(bgc_gbk) as ogbf:
				for line in ogbf:
					line = line.strip()
					if line.startswith('Orig. start') and '::' in line:
						bgc_starts.append(int(line.split()[-1].replace('>', '').replace('<', '')))			

			bgc_scaffold = None
			bgc_length = 0
			with open(bgc_gbk) as obg:
				for rec in SeqIO.parse(obg, 'genbank'):
					bgc_scaffold = rec.id
					bgc_length = len(str(rec.seq))

			full_scaff_length = None
			try:
				scaff_rec_with_bgc = full_genome_index[rec.id]
				full_scaff_length = len(str(scaff_rec_with_bgc.seq))		
				assert(full_scaff_length != None)				
			except Exception as e2:
				msg = f'Issue finding the location of the BGC {bgc_gbk} for sample {sample}.'
				raise RuntimeWarning(msg) from e2
				return

			with open(outf, 'w') as outf_handle:
				for start_coord in bgc_starts:
					start_coord += 1
					end_coord = start_coord + bgc_length - 1
					outf_handle.write(bgc_scaffold + '\t' + str(start_coord) + '\t' + str(end_coord) + '\t' + str(min([start_coord-1, full_scaff_length-end_coord])) + '\n')

	except Exception as e:
		raise RuntimeWarning('Issue determining the location of antiSMASH BGC(s) in the full genome GenBank.') from e

def processAntiSMASHBGCtoGenomeMappingResults(process_data, logObject):
	"""
	Definition:
	Function to process antiSMASH BGC to genome mapping results and return the info as a list.
	********************************************************************************************************************
	Parameters:
	- process_data: A list of lists containing paths to files with BGC to genome mapping information.
	- logObject: A python logging object.
	********************************************************************************************************************
	Returns:
	- bgc_locations: A list of lists where the first item in a sublist is the bgc identifier/path followed by which 
	                 sample it belongs to, the scaffold it is on, the start and end coordinates, the length of the bgc,
					 and the distance of the BGC start/end to the nearest scaffold edge.
	********************************************************************************************************************
	"""
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
				msg = f'Warning: no location or more than 2 potential locations for BGC region {bgc}. Skipping incorporation.'
				logObject.warning(msg)
				sys.stderr.write(msg + '\n')
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
		raise RuntimeError(msg) from e

def extractProteinsFromGenBank(inputs):
	"""
	Definition:
	Void function to extract protein sequences from CDS features in a GenBank file and output in FASTA format.
	********************************************************************************************************************
	Parameters:
	- inputs:
		- input_genbank: The input GenBank file with CDS features and translation + locus_tag qualifiers 
		- output_proteome:	The output FASTA proteome file
		- logObject: A python logging object
	********************************************************************************************************************
	"""
	input_genbank, output_proteome, logObject = inputs
	try:
		with open(output_proteome, 'w') as op_handle, open(input_genbank) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						try:
							prot_lt = feature.qualifiers.get('locus_tag')[0]
							prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
							op_handle.write('>' + prot_lt + '\n' + prot_seq + '\n')
						except Exception as e:
							msg = f'Warning: no locus_tag or protein_id for feature {feature.type} in GenBank file {input_genbank}. Skipping.'
							sys.stderr.write(msg + '\n')
							logObject.warning(msg)
							pass
	except Exception as e:
		logObject.error("Issues with parsing out protein sequences from GenBank file %s." % input_genbank)
		logObject.error(traceback.format_exc())
		raise RuntimeError("Issues with parsing out protein sequences from GenBank file %s." % input_genbank) from e

def determineAsofName(asof_index):
	"""
	Simple function to make integers into strings of equal length by appending 0s as prefixes.
	"""
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
		asof_name = asof_index_str
	assert(asof_name != None)
	return(asof_name)

def runOrthoFinder2Full(prot_directory, orthofinder_outdir, run_msa, logObject, threads=1):
	"""
	Definition:
	Wrapper function to run OrthoFinder2 for bacteria (involves some extra processing of hierarchical orthogroups).
	********************************************************************************************************************
	Parameters:
	- prot_directory: The input directory of proteomes in FASTA format to provide to OrthoFinder2 as input.
	- orthofinder_outdir: The workspace/output directory for OrthoFinder2 analysis.
	- run_msa: Flag on whether to infer hierarchical orthogroups using MSA instead of DendroBlast.
	- logObject: A python logging object.
	- threads (optional): The number of threads. Default is 1. 
	********************************************************************************************************************
	Returns:
	- result_file: The resulting orthogroup by sample carriage matrix.
	********************************************************************************************************************
	"""
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

		shutil.move(tmp_orthofinder_dir, orthofinder_outdir)
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		n0_file = orthofinder_outdir + 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
		putative_xenologs_dir = orthofinder_outdir + 'Putative_Xenologs/'
		phylo_misplaced_genes_dir = orthofinder_outdir + 'Phylogenetically_Misplaced_Genes/'
		
		with open(result_file, 'w') as result_handle:

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
				hog_missing_to_add[top_hog][fg_to_genome[fg]].add(fg)

		genomes = []
		genome_genes_accounted = defaultdict(set)
		with open(result_file, 'w') as result_handle:
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

		assert (genomes == genomes_mf and genomes == genomes_sf)
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd)) from e
	return result_file


def runOrthoFinder2FullFungal(prot_directory, orthofinder_outdir, run_msa, logObject, threads=1):
	"""
	Definition:
	Wrapper function to run OrthoFinder2 for fungi.
	********************************************************************************************************************
	Parameters:
	- prot_directory: The input directory of proteomes in FASTA format to provide to OrthoFinder2 as input.
	- orthofinder_outdir: The workspace/output directory for OrthoFinder2 analysis.
	- run_msa: Flag on whether to infer hierarchical orthogroups using MSA instead of DendroBlast.
	- logObject: A python logging object.
	- threads (optional): The number of threads. Default is 1.
	********************************************************************************************************************
	Returns:
	- result_file: The resulting orthogroup by sample carriage matrix.
	********************************************************************************************************************
	"""

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

		shutil.move(tmp_orthofinder_dir, orthofinder_outdir)
		n0_file = orthofinder_outdir + 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
		
		with open(result_file, 'w') as result_handle:

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
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd)) from e
	return result_file


def runOrthoFinder2(prot_directory, orthofinder_outdir, logObject, threads=1):
	"""
	Definition:
	Wrapper function to run OrthoFinder2 up until the coarse orthogroups are determined (basically the end of 
	OrthoFinder 1).
	********************************************************************************************************************
	Parameters:
	- prot_directory: The input directory of proteomes in FASTA format to provide to OrthoFinder2 as input.
	- orthofinder_outdir: The workspace/output directory for OrthoFinder2 analysis.
	- logObject: A python logging object.
	- threads (optional): The number of threads. Default is 1. 
	********************************************************************************************************************
	Returns:
	- result_file: The resulting orthogroup by sample carriage matrix.
	********************************************************************************************************************
	"""
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
		shutil.move(tmp_orthofinder_dir, orthofinder_outdir)
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		with open(result_file, 'w') as result_handle:
			with open(main_file) as omf:
				for line in omf:
					result_handle.write(line)
			with open(singletons_file) as osf:
				for i, line in enumerate(osf):
					if i == 0: continue
					result_handle.write(line)
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd)) from e
	return result_file

def runPanaroo(detailed_BGC_listing_file, panaroo_input_dir, results_directory, logObject, threads=1, panaroo_options='--clean-mode moderate --remove-invalid-genes'):
	"""
	Definition:
	Wrapper function to run Panaroo.
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_file: A detailed BGC listing file with positional/location information for each BGC.
	- panaroo_input_dir: The directory with input GFF files (in Prokka format) to provide to Panaroo as input. 
	- results_directory: The workspace/output directory for Panaroo analysis.
	- logObject: A python logging object.
	- threads (optional): . Default is 1. 
	- panaroo_options (optional): Additional options for running Panaroo. Default is "--clean-mode moderate 
	                              --remove-invalid-genes". 
	********************************************************************************************************************
	Returns:
	- result_file: The resulting orthogroup by sample carriage matrix.
	********************************************************************************************************************
	"""
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
				submitted_jobs = [executor.submit(multiProcess, arg) for arg in reformat_cmds]
				for future in tqdm(concurrent.futures.as_completed(submitted_jobs), total=len(submitted_jobs), desc="Creating Prokka-like GFF files for Panaroo"):
					future.result()
		except Exception as e:
			msg = 'Issues creating Prokka-like GFF files for Panaroo.'
			sys.stderr.write(msg + '\n')
			sys.stderr.write(traceback.format_exc() + '\n')
			logObject.error(msg)
			logObject.error(traceback.format_exc())
			sys.exit(1)

		panaroo_cmd = ['panaroo', '-t', str(threads), panaroo_options, '-i', ' '.join(panaroo_inputs), '-o', results_directory]

		try:
			logObject.info('Running the following command: %s' % ' '.join(panaroo_cmd))
			subprocess.call(' '.join(panaroo_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran Panaroo!')
		except Exception as e:
			msg = 'Issue running Panaroo with the following command: %s' % ' '.join(panaroo_cmd)
			sys.stderr.write(msg + '\n')
			sys.stderr.write(traceback.format_exc() + '\n')
			logObject.error(msg)
			logObject.error(traceback.format_exc())
			sys.exit(1)

		main_ortho_file = results_directory + 'gene_presence_absence.csv'
		
		with open(result_file, 'w') as result_handle:
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
							try:
								lt = feat.qualifiers.get('locus_tag')[0]
							except Exception as e:
								msg = f'Warning: no locus tag for feature {feat.type} in GenBank file {gbk}. Skipping.'
								sys.stderr.write(msg + '\n')
								continue

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

		assert (os.path.isfile(result_file))
	except Exception as e:
		msg = "Problem with determining orthogroups using Panaroo."
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)
	return result_file

def runBofasa(detailed_BGC_listing_file, bofasa_input_dir, bofasa_prep_dir, results_directory, logObject, threads=1, ignore_upperbound_limit=False, rooting_seeds=1):
	"""
	Definition:
	Wrapper function to run bofasa (https://github.com/raufs/bofasa) for multi-species bacterial orthology inference.
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_file: A detailed BGC listing file with positional/location information for each BGC.
	- bofasa_input_dir: Directory where symlinks to the full-genome GenBank files (named by sample) will be created.
	- bofasa_prep_dir: The output directory for "bofasa prep".
	- results_directory: The output directory for "bofasa run".
	- logObject: A python logging object.
	- threads (optional): The number of threads to use. Default is 1.
	- ignore_upperbound_limit (optional): Pass --ignore-upperbound-limit to bofasa prep to allow > 200 genomes.
	- rooting_seeds (optional): The maximum number of nodes bofasa run tries for rooting (-rs). Default is 1.
	********************************************************************************************************************
	Returns:
	- result_file: The resulting orthogroup by sample carriage matrix.
	********************************************************************************************************************
	"""
	result_file = results_directory + 'Final_Orthogroups.tsv'
	try:
		sample_genomes = {}
		with open(detailed_BGC_listing_file) as odlf:
			for i, line in enumerate(odlf):
				if i == 0: continue
				line = line.strip()
				sample, _, genome_gbk = line.split('\t')[:3]
				sample_genomes[sample] = genome_gbk

		# bofasa uses the GenBank file basename (minus extension) as the sample name, so create symlinks named by
		# sample to make sure columns in the resulting matrix match lsaBGC-Pan sample names.
		bofasa_inputs = []
		for sample in sorted(sample_genomes):
			link_path = bofasa_input_dir + sample + '.gbk'
			os.symlink(os.path.abspath(sample_genomes[sample]), link_path)
			bofasa_inputs.append(link_path)

		# locus tags in the lsaBGC-Pan GenBank files are retained since --rename-locus-tags is not specified.
		prep_cmd = ['bofasa', 'prep', '-i', ' '.join(bofasa_inputs), '-o', bofasa_prep_dir, '-c', str(threads), '-y']
		if ignore_upperbound_limit:
			# bofasa prep has its own upper limit of 200 genomes
			prep_cmd.append('--ignore-upperbound-limit')
		runCmdViaSubprocess(prep_cmd, logObject=logObject, check_directories=[bofasa_prep_dir],
							stdout=bofasa_prep_dir.rstrip('/') + '.stdout.log', stderr=bofasa_prep_dir.rstrip('/') + '.stderr.log')

		bofasa_ogs_file = results_directory + 'Final_Results/Protein_Ortholog_Groups.tsv'
		run_cmd = ['bofasa', 'run', '-i', bofasa_prep_dir, '-o', results_directory, '-c', str(threads), '-rs', str(rooting_seeds), '-y']
		runCmdViaSubprocess(run_cmd, logObject=logObject, check_files=[bofasa_ogs_file],
							stdout=results_directory.rstrip('/') + '.stdout.log', stderr=results_directory.rstrip('/') + '.stderr.log')

		# Reformat to the standard lsaBGC-Pan orthogroup by sample matrix. bofasa should report singletons, but
		# any CDS features not accounted for (e.g. very short proteins filtered by bofasa) are added as singletons.
		with open(result_file, 'w') as result_handle:
			clustered_lts = set([])
			genomes = []
			cds_iter = 1
			with open(bofasa_ogs_file) as obof:
				for i, line in enumerate(obof):
					line = line.rstrip('\n')
					ls = line.split('\t')
					if i == 0:
						genomes = ls[1:]
						result_handle.write('\t'.join(['Orthogroup'] + genomes) + '\n')
					else:
						for lts in ls[1:]:
							for lt in lts.split(', '):
								lt = lt.strip()
								if lt != '': clustered_lts.add(lt)
						result_handle.write('\t'.join(['OG' + determineAsofName(cds_iter)] + ls[1:]) + '\n')
						cds_iter += 1

			for sample in sample_genomes:
				if not sample in genomes:
					raise RuntimeError(f'Sample {sample} not found as a column in bofasa results: {bofasa_ogs_file}')
				gbk = sample_genomes[sample]
				with open(gbk) as ogbk:
					for rec in SeqIO.parse(ogbk, 'genbank'):
						for feat in rec.features:
							if not feat.type == 'CDS': continue
							try:
								lt = feat.qualifiers.get('locus_tag')[0]
							except Exception as e:
								msg = f'Warning: no locus tag for feature {feat.type} in GenBank file {gbk}. Skipping.'
								sys.stderr.write(msg + '\n')
								continue

							if not lt in clustered_lts:
								printlist = ['OG' + determineAsofName(cds_iter)]
								cds_iter += 1
								for genome in genomes:
									if genome == sample:
										printlist.append(lt)
									else:
										printlist.append('')
								result_handle.write('\t'.join(printlist) + '\n')

		assert (os.path.isfile(result_file))
	except Exception as e:
		msg = "Problem with determining orthogroups using bofasa."
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)
	return result_file

def processPrecomputedGECCOResults(gecco_precomputed_dir, sample_genomes, gecco_results_dir, logObject):
	"""
	Definition:
	Function to process precomputed GECCO results (e.g. from nf-core/funcscan) for use in lsaBGC-Pan. Only the BGC
	coordinates and types (from the *.clusters.tsv files) and the Pfam domain annotations (from the cluster GenBank
	files) are used. BGC GenBank files are recreated from the lsaBGC-Pan gene-called full genome GenBank files so that
	CDS features and locus tags are consistent with the rest of the analysis.

	The precomputed GECCO directory is expected to have one subdirectory per sample, named the same as the genome
	FASTA file (minus extension) provided to lsaBGC-Pan, each containing GECCO "run" results for that same FASTA.
	********************************************************************************************************************
	Parameters:
	- gecco_precomputed_dir: Directory with per-sample subdirectories of precomputed GECCO results.
	- sample_genomes: Dictionary mapping sample names to lsaBGC-Pan's gene-called full genome GenBank files.
	- gecco_results_dir: Directory where recreated BGC GenBank files will be written (in per-sample subdirectories).
	- logObject: A python logging object.
	********************************************************************************************************************
	Returns:
	- gecco_bgc_locations: List of BGC location information lists: [bgc_file, sample, scaffold, start, end,
	                       bgc_length, dist_to_edge, 'gecco'].
	********************************************************************************************************************
	"""
	gecco_bgc_locations = []
	try:
		# map samples to precomputed GECCO subdirectories
		sample_to_gecco_dir = {}
		for sd in os.listdir(gecco_precomputed_dir):
			subdir = gecco_precomputed_dir + sd + '/'
			if not os.path.isdir(subdir): continue
			# mirror the sample naming used by lsaBGC-Pan for genome FASTA files with troublesome characters.
			sample = sd
			for tc in [' ', '\t', ',', '|', '\\']:
				sample = sample.replace(tc, '_')
			sample = cleanUpSampleName(sample)
			if sample in sample_genomes:
				sample_to_gecco_dir[sample] = subdir
			else:
				msg = 'Warning: precomputed GECCO results subdirectory %s does not match any input genome sample name - ignoring it.' % subdir
				sys.stderr.write(msg + '\n')
				logObject.warning(msg)

		missing_samples = sorted([s for s in sample_genomes if not s in sample_to_gecco_dir])
		if len(missing_samples) > 0:
			raise RuntimeError('No precomputed GECCO results subdirectory found for the following samples (subdirectory names must match genome FASTA file names minus extension): %s' % ', '.join(missing_samples))

		for sample in sorted(sample_genomes):
			gecco_dir = sample_to_gecco_dir[sample]
			cluster_tsvs = [f for f in os.listdir(gecco_dir) if f.endswith('.clusters.tsv')]
			if len(cluster_tsvs) != 1:
				raise RuntimeError('Expected exactly one *.clusters.tsv file in precomputed GECCO results subdirectory %s, found %d.' % (gecco_dir, len(cluster_tsvs)))
			cluster_tsv = gecco_dir + cluster_tsvs[0]

			scaffold_recs = {}
			with open(sample_genomes[sample]) as ogg:
				for rec in SeqIO.parse(ogg, 'genbank'):
					scaffold_recs[rec.id] = rec

			sample_outdir = gecco_results_dir + sample + '/'
			setupReadyDirectory([sample_outdir])

			with open(cluster_tsv) as octf:
				header = []
				for i, line in enumerate(octf):
					ls = line.rstrip('\n').split('\t')
					if i == 0:
						header = ls
						continue
					row = dict(zip(header, ls))
					scaff, bgc_id = ls[0], ls[1]
					start_coord, end_coord = sorted([int(ls[2]), int(ls[3])])
					bgc_type = row.get('type', 'Unknown')
					if not scaff in scaffold_recs:
						raise RuntimeError('Scaffold %s (from %s) not found in genome for sample %s - please make sure GECCO was run on the same genome FASTA file provided to lsaBGC-Pan.' % (scaff, cluster_tsv, sample))
					genome_rec = scaffold_recs[scaff]
					scaff_seq = str(genome_rec.seq).upper()

					# determine Pfam domain annotations and confirm sequence agreement using GECCO cluster GenBank
					domain_features = []
					gecco_gbk = gecco_dir + bgc_id + '.gbk'
					if os.path.isfile(gecco_gbk):
						gecco_rec = SeqIO.read(gecco_gbk, 'genbank')
						gecco_seq = str(gecco_rec.seq).upper()
						window_start = max(0, start_coord - 1 - 10)
						idx = scaff_seq.find(gecco_seq, window_start, min(len(scaff_seq), end_coord + 10))
						if idx == -1:
							raise RuntimeError('Sequence of GECCO BGC %s does not match the expected region (%s:%d-%d) in the genome for sample %s - please make sure GECCO was run on the same genome FASTA file provided to lsaBGC-Pan.' % (gecco_gbk, scaff, start_coord, end_coord, sample))
						start_coord = idx + 1
						end_coord = idx + len(gecco_seq)
						for feat in gecco_rec.features:
							if feat.type == 'misc_feature':
								domain_features.append(feat)
					else:
						msg = 'Warning: GECCO BGC GenBank %s not found - Pfam domain annotations will not be available for this BGC.' % gecco_gbk
						sys.stderr.write(msg + '\n')
						logObject.warning(msg)

					# extend boundaries to fully encompass CDS features from lsaBGC-Pan gene calling which overlap the
					# region partially (gene calls can differ slightly from GECCO's internal gene calling).
					cds_features = []
					region_start, region_end = start_coord - 1, end_coord
					for feat in genome_rec.features:
						if feat.type != 'CDS': continue
						fstart, fend = int(feat.location.start), int(feat.location.end)
						if fend <= start_coord - 1 or fstart >= end_coord: continue
						cds_features.append(feat)
						region_start = min(region_start, fstart)
						region_end = max(region_end, fend)

					if len(cds_features) == 0:
						msg = 'Warning: no CDS features from lsaBGC-Pan gene calling found within GECCO BGC %s for sample %s - skipping it.' % (bgc_id, sample)
						sys.stderr.write(msg + '\n')
						logObject.warning(msg)
						continue

					domain_shift = (start_coord - 1) - region_start
					new_rec = SeqRecord(genome_rec.seq[region_start:region_end], id=scaff, name=scaff[:16], description='GECCO BGC %s recreated by lsaBGC-Pan' % bgc_id)
					new_rec.annotations['molecule_type'] = 'DNA'
					new_rec.annotations['structured_comment'] = OrderedDict([('GECCO-Data', OrderedDict([('cluster_type', bgc_type), ('source', 'precomputed GECCO results processed by lsaBGC-Pan')]))])
					for feat in sorted(cds_features, key=lambda x: int(x.location.start)):
						new_rec.features.append(feat._shift(-region_start))
					for feat in domain_features:
						new_rec.features.append(feat._shift(domain_shift))

					bgc_file = sample_outdir + bgc_id + '.gbk'
					SeqIO.write(new_rec, bgc_file, 'genbank')

					bgc_start, bgc_end = region_start + 1, region_end
					bgc_length = bgc_end - bgc_start + 1
					dist_to_edge = min([bgc_start, len(scaff_seq) - bgc_end])
					gecco_bgc_locations.append([bgc_file, sample, scaff, bgc_start, bgc_end, bgc_length, dist_to_edge, 'gecco'])
	except Exception as e:
		msg = 'Issues processing precomputed GECCO results.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.exit(1)
	return gecco_bgc_locations

def setupReadyDirectory(directories):
	"""
	Simple function to (re-)create directories provided as a list.  
	"""
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				shutil.rmtree(d)
			os.makedirs(d)
	except Exception as e:
		raise RuntimeError('Issue setting up directories: %s' % str(directories)) from e


def p_adjust_bh(p):
	"""
	Benjamini-Hochberg p-value correction for multiple hypothesis testing.
	"""
	p = np.asarray(p, dtype=float)
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
	except Exception as e:
		return False


def is_fastq(fastq):
	"""
	Function to validate if FASTQ file is correctly formatted.
	"""
	try:
		with open(fastq) as of:
			return next(SeqIO.parse(of, 'fastq'), None) is not None
	except Exception as e:
		return False


def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				return next(SeqIO.parse(ogf, 'fasta'), None) is not None
		else:
			with open(fasta) as of:
				return next(SeqIO.parse(of, 'fasta'), None) is not None
	except Exception as e:
		return False


def is_genbank(gbk):
	"""
	Function to check if Genbank file is correctly formatted.
	"""
	try:
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				return next(SeqIO.parse(ogf, 'genbank'), None) is not None
		else:
			with open(gbk) as of:
				return next(SeqIO.parse(of, 'genbank'), None) is not None
	except Exception as e:
		return False
	
def createLoggerObject(log_file):
	"""
	Description:
	Function which creates logger object.
	********************************************************************************************************************	
	Parameters:
	- log_file: path to log file.
	********************************************************************************************************************
	Returns:
	- logObject: logging logger object.
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
	Description:
	Void function which closes/terminates a logger object.
	********************************************************************************************************************	
	Parameters:
	- logObject: A python logger object to close.
	********************************************************************************************************************	
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
	with open(parameter_file, 'w') as parameter_handle:
		for i, pv in enumerate(parameter_values):
			pn = parameter_names[i]
			parameter_handle.write(pn + ': ' + str(pv) + '\n')


def logParametersToObject(logObject, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))

def determineSeqSimProteinAlignment(protein_alignment_file, use_only_core=True):
	"""
	Definition:	
	Function to determine sequence similarity between pairs of proteins in a protein multiple sequence alignment.
	Used in popstrat.
	********************************************************************************************************************
	Parameters:
	- protein_alignment_file: The input protein multiple sequence alignment in FASTA format.
	- use_only_core (optional): Use only core sites shared by both sequences (neither sequence has an ambiguous allele). 
	                            Default is True.
	********************************************************************************************************************
	Returns:
	- pair_seq_matching: A dictionary showing the sequence similarity between proteins in the alignment.
	********************************************************************************************************************
	"""
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

def parseVersion():
	"""
	Returns suite version. 
	"""
	return(package_version)

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
	except Exception as e:
		return float('nan')

def loadTableInPandaDataFrame(input_file, numeric_columns):
	"""
	Description:
	This function reads a TSV file and stores it as a pandas dataframe.
	********************************************************************************************************************
	Parameters:
	- input_file: The input TSV file, with first row corresponding to the header.
	- numeric_columns: Set of column names which should have numeric data.
	********************************************************************************************************************
	Returns:
	- panda_df: A pandas DataFrame object representation of the input TSV file.
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
	"""
	********************************************************************************************************************
	Definition:
	Void, major function to create the final consolidated spreadsheet that is the major result of lsaBGC-Pan.
	********************************************************************************************************************
	Parameters:
	- detailed_BGC_listing_with_pop_and_GCF_map_file: The detailed BGC listing file produced in lsaBGC-Pan with much of 
	                                                  the information for sample mapping, GCF mapping, and context of 
													  each BGC.
	- zol_results_dir:
	- zol_high_qual_flag: A flag for whether zol analyses were run with a more comprehensive/high-quality preset of 
	                      options
	- mibig_dir: The directory with results from lsaBGC-MIBiGMapper for all GCFs.
	- recon_result_file: The major result file from lsaBGC-Recon showing a table showing how BGC-associated orthogroups
	                     are distributed in different populations and GCF contexts.
	- recon_og_result_file: The secondary major result file from lsaBGC-Recon which shows a matrix of orthogroup 
	                        carriage across samples.
	- recon_pop_color_file: The population coloring determined in lsaBGC-Recon to eventually color the spreadsheet 
	                        using. Currently is not used!
	- sociate_result_file: The major result file from lsaBGC-Sociate showing the results of finding associated/
	                       dissociated orthogroups and GCFs with focal GCFs.
	- final_spreadsheet_xlsx: The final consolidated spreadsheet path; to be made in the function.
	- scratch_dir: A scratch space to write temporary/intermediate files.
	- logObject: A python logging object.
	********************************************************************************************************************
	"""

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

		zol_sheet_header = [
			'GCF ID',
			'Ortholog Group (OG) ID',
			'OG is Single Copy?',
			'Proportion of Total Gene Clusters with OG',
			'Proportion of Complete Gene Cluster Instances with OG',
			'OG Median Length (bp)',
			'OG Consensus Order',
			'OG Consensus Direction',
			'PGAP Annotation (E-value)',
			'Pfam Domains',
			'Average Pairwise Nucleotide Identity',
			'Average Pairwise Protein Identity',
			'Std Dev of Pairwise Protein Identity',
			'Protein to Nucleotide Identity Ratio',
			'Tajima\'s D',
			'Proportion of Filtered Codon Alignment is Segregating Sites',
			'Entropy',
			'Upstream Region Entropy',
			'Median Beta-RD-gc',
			'Max Beta-RD-gc',
			'Proportion of sites which are highly ambiguous in codon alignment',
			'Proportion of sites which are highly ambiguous in trimmed codon alignment',
			'Median GC',
			'Median GC Skew',
			'BGC score (GECCO weights)',
			'Viral score (V-Score)',
			'Hydrophobicity Mean',
			'Hydrophobicity Std Dev',
			'Aliphatic Index Mean',
			'Aliphatic Index Std Dev',
			'm/z Mean',
			'm/z Std Dev',
			'KO Annotation (E-value)',
			'PaperBLAST Annotation (E-value)',
			'CARD Annotation (E-value)',
			'IS Finder (E-value)',
			'MIBiG Annotation (E-value)',
			'VOG Annotation (E-value)',
			'VFDB Annotation (E-value)',
			'DefenseFinder Annotation (E-value)'
		]

		# Optional HyPhy fields when zol is run with HyPhy
		if zol_high_qual_flag:
			zol_sheet_header += [
				'GARD Partitions Based on Recombination Breakpoints',
				'Gene-wide dN/dS by MG94',
				'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
				'Average delta(Beta, Alpha) by FUBAR across sites',
				'Gene-wide dN/dS by FUBAR (Sum of Beta / Sum of Alpha)',
				'Proportion of Sites Under Selection which are Positive',
				'P-value for gene-wide episodic selection by BUSTED'
			]

		zol_sheet_header += [
			'CDS Locus Tags',
			'OG Consensus Sequence'
		]

		zol_full_dir = zol_results_dir + 'Comprehensive/'
		gcf_comp_cons_dir = zol_results_dir + 'Complete_Instances/'	

		comp_cons = defaultdict(lambda: defaultdict(lambda: 'NA'))
		for f in os.listdir(gcf_comp_cons_dir):
			gcf_comp_cons_file = gcf_comp_cons_dir + f
			gcf = f.split('.txt')[0]
			if not os.path.isfile(gcf_comp_cons_file): continue			
			with open(gcf_comp_cons_file) as ogrf:
				for line in ogrf:
					line = line.strip()
					if not line: continue
					ls = line.split('\t')
					og = ls[0]
					cons = ls[1]
					comp_cons[gcf][og] = cons

		with open(zol_combined_tsv_file, 'w') as zctf_handle:
			zctf_handle.write('\t'.join(zol_sheet_header) + '\n')
			num_rows = 1
			# Dynamically parse zol headers and rows to avoid hard-coded indices
			zol_header_from_report = None
			for gcf in os.listdir(zol_full_dir):
				gcf_result_file = zol_full_dir + gcf + '/Final_Results/Consolidated_Report.tsv'
				if not os.path.isfile(gcf_result_file):
					continue
				with open(gcf_result_file) as ogrf:
					for i, line in enumerate(ogrf):
						line = line.strip('\n')
						if i == 0:
							zol_header_from_report = line.split('\t')
							continue
						ls = line.split('\t')
						col_to_val = {}
						for ci, cv in enumerate(ls):
							if zol_header_from_report is not None and ci < len(zol_header_from_report):
								col_to_val[zol_header_from_report[ci]] = cv

						og_id = col_to_val.get('Ortholog Group (OG) ID', 'NA')
						row_vals = []
						for h in zol_sheet_header:
							if h == 'GCF ID':
								row_vals.append(gcf)
							elif h == 'Ortholog Group (OG) ID':
								row_vals.append(og_id)
							elif h == 'Proportion of Complete Gene Cluster Instances with OG':
								row_vals.append(comp_cons[gcf][og_id])
							else:
								row_vals.append(col_to_val.get(h, 'NA'))
						zctf_handle.write('\t'.join([str(x) for x in row_vals]) + '\n')
						num_rows += 1

		zr_numeric_columns = set([
			'Proportion of Total Gene Clusters with OG',
			'Proportion of Complete Gene Cluster Instances with OG',
			'OG Median Length (bp)',
			'OG Consensus Order',
			'Average Pairwise Nucleotide Identity',
			'Average Pairwise Protein Identity',
			'Std Dev of Pairwise Protein Identity',
			'Protein to Nucleotide Identity Ratio',
			'Tajima\'s D',
			'GARD Partitions Based on Recombination Breakpoints',
			'Gene-wide dN/dS by MG94',
			'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
			'Average delta(Beta, Alpha) by FUBAR across sites',
			'Gene-wide dN/dS by FUBAR (Sum of Beta / Sum of Alpha)',
			'Proportion of Sites Under Selection which are Positive',
			'P-value for gene-wide episodic selection by BUSTED',
			'Proportion of Filtered Codon Alignment is Segregating Sites',
			'Entropy',
			'Upstream Region Entropy',
			'Median Beta-RD-gc',
			'Max Beta-RD-gc',
			'Proportion of sites which are highly ambiguous in codon alignment',
			'Proportion of sites which are highly ambiguous in trimmed codon alignment',
			'Median GC',
			'Median GC Skew',
			'BGC score (GECCO weights)',
			'Viral score (V-Score)',
			'Hydrophobicity Mean',
			'Hydrophobicity Std Dev',
			'Aliphatic Index Mean',
			'Aliphatic Index Std Dev',
			'm/z Mean',
			'm/z Std Dev'
		])
		
		zr_data = loadTableInPandaDataFrame(zol_combined_tsv_file, zr_numeric_columns)
		zr_data.to_excel(writer, sheet_name='zol Results', index=False, na_rep="NA")
		zr_sheet =  writer.sheets['zol Results']

		# Apply borders to all cells similar to zol styling
		border_format = workbook.add_format({'border': 1, 'border_color': '#DCDCDC'})
		zr_sheet.set_column(0, len(zr_data.columns) - 1, None, border_format)

		# Header/NA/Warning formatting
		zr_sheet.conditional_format('C2:C' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"False"', 'format': warn_format})
		zr_sheet.conditional_format('A2:CA' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"NA"', 'format': na_format})
		zr_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# Build dynamic column letter map based on current dataframe columns
		def _excel_col_letters(n):
			cols = []
			for i in range(n):
				div, mod = divmod(i, 26)
				letter = chr(65 + mod)
				while div:
					div, mod = divmod(div - 1, 26)
					letter = chr(65 + mod) + letter
				cols.append(letter)
			return cols

		excel_cols = _excel_col_letters(len(zr_data.columns))
		col_map = {col: excel_cols[i] for i, col in enumerate(list(zr_data.columns))}

		# Helper to conditionally apply formatting if column exists
		def _fmt_if(col_name, fmt_dict):
			if col_name in col_map:
				cell_range = f"{col_map[col_name]}2:{col_map[col_name]}{num_rows}"
				zr_sheet.conditional_format(cell_range, fmt_dict)

		# Apply coloring similar to zol v1.8.0 using dynamic mapping
		_fmt_if('Proportion of Total Gene Clusters with OG', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#CCCCCC", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Proportion of Complete Gene Cluster Instances with OG', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#CCCCCC", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('OG Median Length (bp)', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#5A8AC6", "min_value": 100, "max_value": 2500, "min_type": "num", "max_type": "num"})

		# Pairwise identity statistics (added in zol v1.7.0; percentages) - coloring matches zol
		for identity_col in ['Average Pairwise Nucleotide Identity', 'Average Pairwise Protein Identity']:
			_fmt_if(identity_col, {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#A5D6E3", "min_value": 40.0, "max_value": 90.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Std Dev of Pairwise Protein Identity', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F86B6B", "min_value": 0.0, "max_value": 100.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Protein to Nucleotide Identity Ratio', {"type": "3_color_scale", "min_color": "#E6B0AA", "mid_color": "#FFFFFF", "max_color": "#B8CCE4", "min_value": 0.5, "mid_value": 1.0, "max_value": 1.5, "min_type": "num", "mid_type": "num", "max_type": "num"})

		_fmt_if('Tajima\'s D', {"type": "3_color_scale", "min_color": "#E6B0AA", "mid_color": "#FFFFFF", "max_color": "#B8CCE4", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0, "min_type": "num", "mid_type": "num", "max_type": "num"})
		_fmt_if('Proportion of Filtered Codon Alignment is Segregating Sites', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#A37BA8", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Entropy', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F86B6B", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Upstream Region Entropy', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F86B6B", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})

		# Beta-RD: apply across both columns if both present, else individually
		if 'Median Beta-RD-gc' in col_map and 'Max Beta-RD-gc' in col_map:
			cell_range = f"{col_map['Median Beta-RD-gc']}2:{col_map['Max Beta-RD-gc']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "3_color_scale", "min_color": "#F8696B", "mid_color": "#FFEB84", "max_color": "#63BE7B", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25, "min_type": "num", "mid_type": "num", "max_type": "num"})
		else:
			_fmt_if('Median Beta-RD-gc', {"type": "3_color_scale", "min_color": "#F8696B", "mid_color": "#FFEB84", "max_color": "#63BE7B", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25, "min_type": "num", "mid_type": "num", "max_type": "num"})
			_fmt_if('Max Beta-RD-gc', {"type": "3_color_scale", "min_color": "#F8696B", "mid_color": "#FFEB84", "max_color": "#63BE7B", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25, "min_type": "num", "mid_type": "num", "max_type": "num"})

		# HyPhy (zol high-quality preset) columns coloring to match zol v1.8.0
		if 'GARD Partitions Based on Recombination Breakpoints' in col_map:
			cell_range = f"{col_map['GARD Partitions Based on Recombination Breakpoints']}2:{col_map['GARD Partitions Based on Recombination Breakpoints']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F8696B", "min_value": 1, "max_value": 5, "min_type": "num", "max_type": "num"})
		if 'Number of Sites Identified as Under Positive or Negative Selection by FUBAR' in col_map:
			cell_range = f"{col_map['Number of Sites Identified as Under Positive or Negative Selection by FUBAR']}2:{col_map['Number of Sites Identified as Under Positive or Negative Selection by FUBAR']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#B8D4B8", "min_value": 0, "max_value": 10, "min_type": "num", "max_type": "num"})
		# Gene-wide dN/dS (added in zol v1.8.0) - coloring matches zol
		for dnds_col in ['Gene-wide dN/dS by MG94', 'Gene-wide dN/dS by FUBAR (Sum of Beta / Sum of Alpha)']:
			_fmt_if(dnds_col, {"type": "3_color_scale", "min_color": "#E6B0AA", "mid_color": "#FFFFFF", "max_color": "#B8D4B8", "min_value": 0, "mid_value": 1, "max_value": 2, "min_type": "num", "mid_type": "num", "max_type": "num"})
		if 'Average delta(Beta, Alpha) by FUBAR across sites' in col_map:
			cell_range = f"{col_map['Average delta(Beta, Alpha) by FUBAR across sites']}2:{col_map['Average delta(Beta, Alpha) by FUBAR across sites']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "3_color_scale", "min_color": "#E6B0AA", "mid_color": "#FFFFFF", "max_color": "#B8D4B8", "min_value": -5, "mid_value": 0, "max_value": 5, "min_type": "num", "mid_type": "num", "max_type": "num"})
		if 'Proportion of Sites Under Selection which are Positive' in col_map:
			cell_range = f"{col_map['Proportion of Sites Under Selection which are Positive']}2:{col_map['Proportion of Sites Under Selection which are Positive']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "3_color_scale", "min_color": "#E6B0AA", "mid_color": "#FFFFFF", "max_color": "#B8D4B8", "min_value": 0, "mid_value": 0.5, "max_value": 1, "min_type": "num", "mid_type": "num", "max_type": "num"})
		if 'P-value for gene-wide episodic selection by BUSTED' in col_map:
			cell_range = f"{col_map['P-value for gene-wide episodic selection by BUSTED']}2:{col_map['P-value for gene-wide episodic selection by BUSTED']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "2_color_scale", "min_color": "#63BE7B", "max_color": "#FFFFFF", "min_value": 0, "max_value": 0.05, "min_type": "num", "max_type": "num"})

		_fmt_if('Proportion of sites which are highly ambiguous in codon alignment', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#E6B0AA", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Proportion of sites which are highly ambiguous in trimmed codon alignment', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#E6B0AA", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		_fmt_if('Median GC', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#63BE7B", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"})
		# Median GC Skew: 3-color scale centered at 0
		if 'Median GC Skew' in col_map:
			cell_range = f"{col_map['Median GC Skew']}2:{col_map['Median GC Skew']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "3_color_scale", "min_color": "#E8D5D5", "mid_color": "#FFFFFF", "max_color": "#D5E8F5", "min_value": -0.2, "mid_value": 0.0, "max_value": 0.2, "min_type": "num", "mid_type": "num", "max_type": "num"})
		# BGC score (GECCO weights): highlight > 2
		gecco_format = workbook.add_format({"bg_color": "#d5abde", "border": 1, "border_color": "#DCDCDC"})
		if 'BGC score (GECCO weights)' in col_map:
			cell_range = f"{col_map['BGC score (GECCO weights)']}2:{col_map['BGC score (GECCO weights)']}{num_rows}"
			zr_sheet.conditional_format(cell_range, {"type": "cell", "criteria": ">", "value": 2, "format": gecco_format})
		# Viral score
		_fmt_if('Viral score (V-Score)', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#E6B0AA", "min_value": 0, "max_value": 4, "min_type": "num", "max_type": "num"})
		
		_fmt_if('Hydrophobicity Mean', {"type": "3_color_scale", "min_color": "#5A8AC6", "mid_color": "#FFFFFF", "max_color": "#F8696B", "min_value": -2.5, "mid_value": 0, "max_value": 2.5, "min_type": "num", "mid_type": "num", "max_type": "num"})
		_fmt_if('Hydrophobicity Std Dev', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F8696B", "min_value": 0, "max_value": 2, "min_type": "num", "max_type": "num"})
		_fmt_if('Aliphatic Index Mean', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#63BE7B", "min_value": 50, "max_value": 150, "min_type": "num", "max_type": "num"})
		_fmt_if('Aliphatic Index Std Dev', {"type": "2_color_scale", "min_color": "#FFFFFF", "max_color": "#F8696B", "min_value": 0, "max_value": 30, "min_type": "num", "max_type": "num"})
		_fmt_if('m/z Mean', {"type": "2_color_scale", "min_color": "#dcedde", "max_color": "#87bba2", "min_value": 5000, "max_value": 75000, "min_type": "num", "max_type": "num"})
		_fmt_if('m/z Std Dev', {"type": "2_color_scale", "min_color": "#dcedde", "max_color": "#87bba2", "min_value": 0, "max_value": 5000, "min_type": "num", "max_type": "num"})

		# Add autofilter across header row and default filter for Proportion >= 0.1 if present
		zr_sheet.autofilter(f"A1:{excel_cols[len(zr_data.columns) - 1]}{num_rows}")
		if 'Proportion of Total Gene Clusters with OG' in col_map:
			zr_sheet.filter_column(col_map['Proportion of Total Gene Clusters with OG'], "x >= 0.1")

		# create MIBiG mapping spreadsheet
		mibig_json_tar_url = 'https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz'
		mibig_json_tar_dir = scratch_dir + 'mibig_json_4.0/'
		mibig_json_tar_file = scratch_dir + 'mibig_json_4.0.tar.gz'

		if not os.path.isdir(mibig_json_tar_dir):
			if os.path.isfile(mibig_json_tar_file):
				os.remove(mibig_json_tar_file)
			wget_cmd = ['curl', '-L', mibig_json_tar_url, '-o', mibig_json_tar_file]

			try:
				subprocess.call(' '.join(wget_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
				assert (os.path.isfile(mibig_json_tar_file))
			except Exception as e:
				sys.stderr.write('Had an issue running: %s\n' % ' '.join(wget_cmd))
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)

			uncompress_cmd = ['tar', '-zxvf', mibig_json_tar_file, '-C', scratch_dir]
			try:
				os.mkdir(mibig_json_tar_dir)
				subprocess.call(' '.join(uncompress_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
				assert (os.path.isdir(mibig_json_tar_dir))
			except Exception as e:
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
					for comp in data['compounds']:
						for key in comp:
							if key == 'compound':
								mibig_bgc_compounds[bgc].append(comp[key])


		mb_combined_tsv_file = scratch_dir + 'mibig_map_results.tsv'
		mb_header = ['GCF ID', 'MIBiG BGC ID', 'GCF OG ID', 'MIBiG Protein Matching', 'MIBiG Compound(s)']
		with open(mb_combined_tsv_file, 'w') as mb_outf:
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

		num_rows = 0
		with open(recon_og_result_file) as oorf:
			for i, line in enumerate(oorf):
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
					  'phylogenetically corrected pvalue', 'beta', 'beta-std-err', 'variant_h2'])
		so_data = loadTableInPandaDataFrame(sociate_result_file, so_numeric_columns)
		so_data.to_excel(writer, sheet_name='lsaBGC-Sociate Results', index=False, na_rep="NA")
		so_sheet =  writer.sheets['lsaBGC-Sociate Results']
		so_sheet.conditional_format('A1:CA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		writer.close()
		
	except Exception as e:
		logObject.error("Issues creating consolidated results files.")
		logObject.error(e)
		logObject.error(traceback.format_exc())
		sys.exit(1)
