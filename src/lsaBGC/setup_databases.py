"""
setup_databases: implementation of "lsaBGC-Pan setup".

Sets up the annotation databases used by zol (via zol's setup_annotation_dbs.py) and then makes the databases bofasa
needs within lsaBGC-Pan (Pfam and ISfinder - geNomad is not used by lsaBGC-Pan) available in the bofasa database
directory (BOFASA_DB_PATH) as symbolic links to zol's copies, so the same databases are not downloaded/stored twice.
bofasa's own database manifest (database_location_paths.txt) is then written, since it has a different layout than
zol's manifest.
"""

import os
import sys
import shutil
import argparse
import subprocess
import traceback

# ISfinder protein sequences - the same pinned release used by "bofasa setup" (only needed if zol was set up in
# minimal mode, which does not include ISfinder).
ISFINDER_URL = ("https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/"
				"bedbeb5c16618eca3e16ef74adf865d6b88cafb5/IS.faa")

BOFASA_MANIFEST_HEADER = 'Database\tType\tPath\tAdditional_Info\n'


def create_parser(argv):
	parser = argparse.ArgumentParser(prog='lsaBGC-Pan setup', formatter_class=argparse.RawTextHelpFormatter,
									 description="""
	Program: lsaBGC-Pan setup
	Description: Sets up the annotation databases used by lsaBGC-Pan.

	1. Runs zol's setup_annotation_dbs.py to download the annotation databases used by zol (into ZOL_DATA_PATH).
	2. Links the databases needed by bofasa within lsaBGC-Pan (Pfam and ISfinder; geNomad is not needed) from zol's
	   database directory into bofasa's database directory (BOFASA_DB_PATH) to avoid storing duplicate databases.

	ZOL_DATA_PATH and BOFASA_DB_PATH are set automatically when the lsaBGC-Pan conda environment is activated, they
	can also be set manually to use other locations (they should be separate, dedicated directories).

	Note, by default zol's comprehensive set of annotation databases is downloaded, which requires considerable disk
	space. Use --minimal to only download what lsaBGC-Pan needs (Pfam, PGAP, and MIBiG for zol - ISfinder is then
	set up separately for bofasa).
	""")
	parser.add_argument('-ld', '--lsabgc-minimal', '-m', '--minimal', dest='minimal', action='store_true', help="Only download the minimal databases for lsaBGC-Pan (Pfam, PGAP,\nand MIBiG) - passes -ld/--lsabgc-minimal to zol's\nsetup_annotation_dbs.py. ISfinder is then set up directly in the\nbofasa database directory.", required=False, default=False)
	parser.add_argument('-l', '--link-only', action='store_true', help="Do not (re-)run zol's database setup, only link the existing zol\ndatabases into the bofasa database directory.", required=False, default=False)
	parser.add_argument('-c', '--threads', type=int, help="Number of threads to use [Default is 4].", required=False, default=4)
	parser.add_argument('-ds', '--download-segments', type=int, help="Number of parallel connections zol uses to download each large file\n(if the server supports it; passed to setup_annotation_dbs.py -ds)\n[Default is 8; use 1 to disable].", required=False, default=8)
	parser.add_argument('-y', '--yes', action='store_true', help="Do not ask for confirmation before zol's setup clears and\nrecreates ZOL_DATA_PATH.", required=False, default=False)
	return parser.parse_args(argv)


def _resolve_dir(env_var, default_subdir):
	"""Get a database directory from an environment variable, falling back to <environment prefix>/share/..."""
	path = os.getenv(env_var, '').strip()
	if not path:
		prefix = os.getenv('CONDA_PREFIX', '').strip() or sys.prefix
		path = os.path.join(prefix, 'share', default_subdir)
		sys.stdout.write('%s is not set, using: %s\n' % (env_var, path))
	return os.path.abspath(path)


def _is_within(path, parent):
	path, parent = os.path.realpath(path), os.path.realpath(parent)
	return path == parent or path.startswith(parent.rstrip(os.sep) + os.sep)


def read_zol_manifest(zol_db_dir):
	"""Read zol's database_location_paths.txt: name -> [descriptions_file, database_path, additional_info]."""
	manifest = os.path.join(zol_db_dir, 'database_location_paths.txt')
	entries = {}
	if not os.path.isfile(manifest):
		return entries
	with open(manifest) as ozm:
		for line in ozm:
			ls = line.rstrip('\n').split('\t')
			if len(ls) >= 3 and ls[0]:
				entries[ls[0]] = ls[1:]
	return entries


def read_bofasa_manifest(bofasa_db_dir):
	"""Read bofasa's database_location_paths.txt: name -> full line fields (Database, Type, Path, Additional_Info)."""
	manifest = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
	entries = {}
	if not os.path.isfile(manifest):
		return entries
	with open(manifest) as obm:
		for i, line in enumerate(obm):
			if i == 0:
				continue
			ls = line.rstrip('\n').split('\t')
			if len(ls) >= 3 and ls[0]:
				entries[ls[0]] = ls
	return entries


def _replace_with_symlink(link_path, target_path):
	"""Create (or replace an existing file/link at) link_path with a symbolic link to target_path."""
	os.makedirs(os.path.dirname(link_path), exist_ok=True)
	if os.path.islink(link_path) or os.path.isfile(link_path):
		if not os.path.islink(link_path):
			sys.stdout.write('Replacing existing (duplicate) database file with a link: %s\n' % link_path)
		os.remove(link_path)
	elif os.path.isdir(link_path):
		raise RuntimeError('Expected a file but found a directory at: %s' % link_path)
	os.symlink(target_path, link_path)


def _count_hmms(hmm_file):
	z = 0
	with open(hmm_file) as ohf:
		for line in ohf:
			if line.startswith('NAME'):
				z += 1
	return z


def setup_isfinder_for_bofasa(isfinder_dir, threads):
	"""Download ISfinder and create its DIAMOND database (as done by "bofasa setup"). Returns the .dmnd path."""
	from Bio import SeqIO
	os.makedirs(isfinder_dir, exist_ok=True)
	isfinder_faa = os.path.join(isfinder_dir, 'ISfinder.faa')
	isfinder_proc_faa = os.path.join(isfinder_dir, 'ISfinder_proc.faa')
	isfinder_dmnd = os.path.join(isfinder_dir, 'ISfinder.dmnd')
	if os.path.islink(isfinder_dmnd):
		os.remove(isfinder_dmnd)

	sys.stdout.write('Downloading ISfinder protein sequences from %s ...\n' % ISFINDER_URL)
	subprocess.run(['curl', '-fsSL', '-o', isfinder_faa, ISFINDER_URL], check=True)
	with open(isfinder_faa, errors='replace') as oif:
		first_line = next((x.strip() for x in oif if x.strip()), '')
	if not first_line.startswith('>'):
		raise RuntimeError('Downloaded ISfinder file does not look like a FASTA file: %s' % isfinder_faa)

	with open(isfinder_faa) as oif, open(isfinder_proc_faa, 'w') as oipf:
		for rec in SeqIO.parse(oif, 'fasta'):
			if '~~~Passenger Gene~~~' in rec.description:
				continue
			oipf.write('>%s\n%s\n' % (rec.description, str(rec.seq)))
	subprocess.run(['diamond', 'makedb', '--in', isfinder_proc_faa, '--db', isfinder_dmnd, '--threads', str(threads)],
				   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if not os.path.isfile(isfinder_dmnd):
		raise RuntimeError('DIAMOND database for ISfinder was not created: %s' % isfinder_dmnd)
	# only the DIAMOND database is used by bofasa
	for f in [isfinder_faa, isfinder_proc_faa]:
		os.remove(f)
	return isfinder_dmnd


def link_zol_databases_for_bofasa(zol_db_dir, bofasa_db_dir, threads=4):
	"""
	Make the Pfam and ISfinder databases from zol's database directory available to bofasa (via symbolic links) and
	write bofasa's database manifest. If zol's setup did not include ISfinder (minimal mode), ISfinder is set up
	directly in the bofasa database directory instead. Returns True if successful.
	"""
	zol_entries = read_zol_manifest(zol_db_dir)
	if 'pfam' not in zol_entries:
		sys.stderr.write('Error: no Pfam database listed in zol\'s database manifest (%s) - please (re-)run '
						 'lsaBGC-Pan setup without --link-only.\n' % os.path.join(zol_db_dir, 'database_location_paths.txt'))
		return False

	os.makedirs(bofasa_db_dir, exist_ok=True)
	previous_entries = read_bofasa_manifest(bofasa_db_dir)

	# Pfam
	zol_pfam_hmm = os.path.abspath(zol_entries['pfam'][1])
	if not os.path.isfile(zol_pfam_hmm):
		sys.stderr.write('Error: zol\'s Pfam database file was not found: %s\n' % zol_pfam_hmm)
		return False
	pfam_z = zol_entries['pfam'][2].strip() if len(zol_entries['pfam']) > 2 else ''
	if not pfam_z.isdigit() or int(pfam_z) <= 0:
		pfam_z = str(_count_hmms(zol_pfam_hmm))
	if int(pfam_z) <= 0:
		sys.stderr.write('Error: no HMMs found in zol\'s Pfam database file: %s\n' % zol_pfam_hmm)
		return False
	bofasa_pfam_hmm = os.path.join(bofasa_db_dir, 'pfam_db', 'Pfam-A.hmm')
	_replace_with_symlink(bofasa_pfam_hmm, zol_pfam_hmm)
	sys.stdout.write('Linked Pfam (%s HMMs): %s -> %s\n' % (pfam_z, bofasa_pfam_hmm, zol_pfam_hmm))

	# ISfinder
	isfinder_dir = os.path.join(bofasa_db_dir, 'isfinder_db')
	bofasa_isfinder_dmnd = os.path.join(isfinder_dir, 'ISfinder.dmnd')
	zol_isfinder_dmnd = os.path.abspath(zol_entries['isfinder'][1]) if 'isfinder' in zol_entries else None
	if zol_isfinder_dmnd and os.path.isfile(zol_isfinder_dmnd):
		_replace_with_symlink(bofasa_isfinder_dmnd, zol_isfinder_dmnd)
		# remove intermediate files from any earlier "bofasa setup" which are no longer needed
		for f in ['ISfinder.faa', 'ISfinder_proc.faa']:
			if os.path.isfile(os.path.join(isfinder_dir, f)) and not os.path.islink(os.path.join(isfinder_dir, f)):
				os.remove(os.path.join(isfinder_dir, f))
		sys.stdout.write('Linked ISfinder: %s -> %s\n' % (bofasa_isfinder_dmnd, zol_isfinder_dmnd))
	else:
		sys.stdout.write('ISfinder is not part of zol\'s database setup (minimal mode) - setting it up for bofasa ...\n')
		try:
			setup_isfinder_for_bofasa(isfinder_dir, threads)
		except Exception:
			sys.stderr.write('Error: issues setting up the ISfinder database for bofasa.\n')
			sys.stderr.write(traceback.format_exc())
			return False
		sys.stdout.write('ISfinder database created: %s\n' % bofasa_isfinder_dmnd)

	# write bofasa's manifest, retaining any previously set up geNomad database (not used by lsaBGC-Pan)
	manifest = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
	with open(manifest, 'w') as obm:
		obm.write(BOFASA_MANIFEST_HEADER)
		if 'genomad' in previous_entries and os.path.isdir(previous_entries['genomad'][2]):
			obm.write('\t'.join(previous_entries['genomad']) + '\n')
		obm.write('isfinder\tdiamond\t%s\t\n' % bofasa_isfinder_dmnd)
		obm.write('pfam\thmm\t%s\t%s\n' % (bofasa_pfam_hmm, pfam_z))
	sys.stdout.write('Wrote bofasa database manifest: %s\n' % manifest)
	return True


def check_bofasa_databases(bofasa_db_dir):
	"""Validate the bofasa database directory the same way bofasa does (Pfam & ISfinder; geNomad is not needed)."""
	entries = read_bofasa_manifest(bofasa_db_dir)
	for db in ['pfam', 'isfinder']:
		if db not in entries or not os.path.isfile(entries[db][2]):
			return False, 'the %s database is missing or its file was not found' % db
	pfam = entries['pfam']
	if len(pfam) < 4 or not pfam[3].strip().isdigit() or int(pfam[3].strip()) <= 0:
		return False, 'the number of Pfam HMMs is missing from the manifest'
	return True, 'bofasa databases are set up'


def main(argv):
	myargs = create_parser(argv)

	zol_db_dir = _resolve_dir('ZOL_DATA_PATH', 'zol/db')
	bofasa_db_dir = _resolve_dir('BOFASA_DB_PATH', 'bofasa/db')
	if _is_within(bofasa_db_dir, zol_db_dir) or _is_within(zol_db_dir, bofasa_db_dir):
		sys.stderr.write('Error: ZOL_DATA_PATH (%s) and BOFASA_DB_PATH (%s) must be separate directories (zol\'s '
						 'setup clears ZOL_DATA_PATH).\n' % (zol_db_dir, bofasa_db_dir))
		return 1

	if myargs.download_segments < 1:
		sys.stderr.write('Error: --download-segments must be at least 1.\n')
		return 1

	sys.stdout.write('zol databases directory: %s\n' % zol_db_dir)
	sys.stdout.write('bofasa databases directory: %s\n' % bofasa_db_dir)

	if not myargs.link_only:
		if shutil.which('setup_annotation_dbs.py') is None:
			sys.stderr.write('Error: zol\'s setup_annotation_dbs.py was not found - is zol installed in this environment?\n')
			return 1
		os.makedirs(zol_db_dir, exist_ok=True)
		zol_cmd = ['setup_annotation_dbs.py', '-c', str(myargs.threads), '-ds', str(myargs.download_segments)]
		if myargs.minimal:
			zol_cmd.append('-ld')
		if myargs.yes:
			zol_cmd.append('-y')
		env = dict(os.environ, ZOL_DATA_PATH=zol_db_dir)
		sys.stdout.write('Running: %s\n' % ' '.join(zol_cmd))
		result = subprocess.run(zol_cmd, env=env)
		if result.returncode != 0:
			sys.stderr.write('Error: zol\'s database setup exited with return code %d.\n' % result.returncode)
			return 1
		issues_file = os.path.join(zol_db_dir, 'issues.txt')
		if os.path.isfile(issues_file) and os.path.getsize(issues_file) > 0:
			with open(issues_file) as oif:
				sys.stderr.write('Warning: zol reported issues setting up some databases:\n' + oif.read())

	if not link_zol_databases_for_bofasa(zol_db_dir, bofasa_db_dir, threads=myargs.threads):
		return 1

	ok, message = check_bofasa_databases(bofasa_db_dir)
	if not ok:
		sys.stderr.write('Error: bofasa database validation failed - %s.\n' % message)
		return 1
	sys.stdout.write('lsaBGC-Pan database setup completed successfully!\n')
	return 0
