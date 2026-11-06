#!/usr/bin/env python3
#------------------------------------------------------------------------------
#  python version of Stntuple/scripts/catalog_stntuples (was ruby)
#
#  usage:
#  ------
#  catalog_stntuples.py -b book -d dataset -m mode -v verbose \
#                       [-i install_directory] \
#                       -u user -H host -D directory -p pattern \
#                       [-n nev_fileset] \
#                       [-o offline_version]
#
#  starts from already created file-level catalog file and assigns fileset
#  numbers. Then creates fileset-level catalog
#
#  access modes:
#  -------------
#  'bluearc' : file in NFS  (default)
#  'xrootd'  : file in PNFS
#
#  1) if pattern is not specified, it is assumed that all filenames contain #{dataset} in them
#  2) the fileset size can be specified, #{nev_fileset}: number of events per fileset, default:200000
#
#  N.B. it is assumed that stntuple files have '.stn' extension
#
#  examples:
#  ---------
#  catalog_stntuples.py -H diskpool -d jbmmgh -b ewk/wtaunu -p jbmm -i /web/sites/mu2e.fnal.gov/htdocs/atwork/computing/Stntuple/cafdfc
#------------------------------------------------------------------------------
import argparse
import os
import re
import shutil
import subprocess
import sys


def to_int(s):
    """Ruby-style String#to_i: leading integer or 0."""
    m = re.match(r'\s*[+-]?\d+', s or '')
    return int(m.group(0)) if m else 0


def to_float(s):
    """Ruby-style String#to_f: leading float or 0.0."""
    m = re.match(r'\s*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?', s or '')
    return float(m.group(0)) if m else 0.0


def field(parts, i):
    return parts[i] if i < len(parts) else None


def shell(cmd):
    """run a shell command and return its stdout"""
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


class Fileset:
    def __init__(self, nevents, lorun, loevt, hirun, hievt):
        self.nevents = nevents
        self.lorun   = lorun
        self.loevt   = loevt
        self.hirun   = hirun
        self.hievt   = hievt
        print(f'new fileset: {nevents} {lorun}  {loevt} {hirun} {hievt}')


def parse_args(argv):
    p = argparse.ArgumentParser(description='catalog stntuple files into filesets', add_help=False)
    p.add_argument('-h', '--help',           action='store_true',    help='print usage and exit')
    p.add_argument('-a', '--access',         default='bluearc',      help='bluearc (NFS, default) or xrootd (PNFS)')
    p.add_argument('-b', '--book',           default='stntuple/dev_242')
    p.add_argument('-c', '--cleanup',        default='yes',          help="'no' keeps the temporary files")
    p.add_argument('-d', '--dataset',        default='',             help='dataset name, sets the default pattern')
    p.add_argument('-D', '--Directory',      default=None,           help='input directory')
    p.add_argument('-H', '--host',           default=None,           help='host (default: this host)')
    p.add_argument('-i', '--install',        nargs='?', const='',    default=None,
                   help='install the catalog; optional install directory')
    p.add_argument('-m', '--mode',           default='1')
    p.add_argument('-n', '--nev_fileset',    type=int, default=200000, help='events per fileset')
    p.add_argument('-o', '--offline_version', default='vXX_XX_XX')
    p.add_argument('-p', '--pattern',        default=None,           help='file name pattern (default: dataset)')
    p.add_argument('-u', '--user',           default=os.environ.get('USER'))
    p.add_argument('-v', '--verbose',        action='store_true')
    return p.parse_args(argv)


def resolve_parameters(a):
    """the host, directory, pattern and protocol defaults of the ruby version"""
    local_host = shell('hostname -f').strip()
    r = {
        'access': a.access, 'book': a.book, 'cleanup': a.cleanup, 'dataset': a.dataset,
        'directory': a.Directory, 'host': a.host or local_host, 'user': a.user,
        'pattern': a.pattern, 'mode': a.mode, 'nev_fileset': a.nev_fileset,
        'offline_version': a.offline_version, 'verbose': a.verbose,
        'install': a.install is not None, 'install_dir': a.install or
                   '/web/sites/mu2e.fnal.gov/htdocs/atwork/computing/Stntuple/cafdfc',
        'diskpool': None, 'protocol': 'rsh', 'local_host': local_host,
    }
    if r['pattern'] is None:
        r['pattern'] = r['dataset']
    if r['pattern'] is None or r['pattern'] == '':
        r['pattern'] = r['dataset']

    host, user, directory, dataset = r['host'], r['user'], r['directory'], r['dataset']
    if host == 'ncdf131':
        r['user'], r['host'] = 'murat', 'ncdf131.fnal.gov'
        if not directory:
            r['directory'] = f'/cdf/data37a/s1/datasets/cdfpewk/{dataset}'
    elif host == 'fcdfdata122':
        r['host'] = 'fcdfdata122.fnal.gov'
        if not directory:
            r['directory'] = f'/export/data2/ewk/datasets/cdfpewk/{dataset}'
    elif r['access'] == 'xrootd':
        if not directory:
            r['directory'] = f'/pnfs/mu2e/scratch/users/{os.environ.get("USER")}/datasets/{dataset}'
    elif r['access'] == 'bluearc':
        if not directory:
            r['directory'] = f'/mu2e/data/users/{os.environ.get("USER")}/datasets/{dataset}'
    elif host == 'diskpool':
        r['diskpool'] = 'cdf://analysis_pool'
        r['protocol'] = 'dcap'
        r['user']     = 'cdfdata'
        r['host']     = 'fcdfrdc3.fnal.gov'
        if not directory:
            r['directory'] = f'/pnfs/diskpool/{r["book"]}/{dataset}'
    elif not directory:
        # Mu2e shortcut: /mu2e/data/users/$USER/datasets/$dataset
        r['directory'] = f'/mu2e/data/users/{os.environ.get("USER")}/datasets/{dataset}'
    return r


def write_request_file(r, request_file):
    """list of the input files, in the format read by dh.C"""
    with open(request_file, 'w') as f:
        if r['host'] == r['local_host']:
            f.write('# DATA_SERVER   .            \n')
            f.write(f'# INPUT_DIR     {r["directory"]}\n')
            f.write(f'# BOOK          {r["book"]}     \n')
            f.write(f'# DATASET       {r["dataset"]}  \n')
            f.write(f'# OUTPUT_DIR    {r["directory"]}\n')
            names = sorted(n for n in os.listdir(r['directory'])
                           if n.endswith('.stn') and r['pattern'] in n)
            for n in names:
                f.write(f'{r["directory"]}/{n}\n')
        else:
            f.write(f'# DATA_SERVER root://{r["host"]}  \n')
            f.write(f'# INPUT_DIR     {r["directory"]}\n')
            f.write(f'# BOOK          {r["book"]}     \n')
            f.write(f'# DATASET       {r["dataset"]}  \n')
            f.write(f'# OUTPUT_DIR    {r["user"]}@{r["host"]}{r["directory"]}\n')
            # the ruby version printed "/filename": awk's dir variable was never set
            out = shell(f'ssh -l {r["user"]} {r["host"]} ls -l {r["directory"]}/ | '
                        f'grep "\\.stn$" | grep {r["pattern"]} | '
                        f"awk '{{print \"{r['directory']}/\"$9}}'")
            f.write(out)


def run_root(r, request_file, pid):
    """run the root macro that writes the file-level catalog into a log file"""
    tmp_script = f'script.{pid}.C'
    logfile    = f'catalog_stntuples.{pid}.log'
    with open(tmp_script, 'w') as f:
        f.write('{\n')
        f.write('  gInterpreter->LoadMacro("Stntuple/scripts/dh.C");\n')
        f.write(f'  gInterpreter->ProcessLine("dh::catalog_list_of_files(\\"{request_file}\\",{r["mode"]})");\n')
        f.write('}\n')
    cmd = f'root.exe -b -q {tmp_script} | grep {r["pattern"]} >> {logfile}  2>&1'
    if r['verbose']:
        print(f'::: {cmd}')
    subprocess.run(cmd, shell=True)
    return tmp_script, logfile


def process_logfile(r, logfile, pid):
    """sort the file-level catalog, assign fileset numbers, write the file catalog"""
    with open(logfile) as f:
        lines = [l for l in f.read().splitlines() if l.strip()]

    # sort by run, then event (the ruby comparator); comment lines sort first
    def key(line):
        p = line.split()
        return (to_int(field(p, 6)), to_int(field(p, 7)))
    lines.sort(key=key)

    file_catalog_fn = r['dataset'].replace('/', '_') + '.files.' + str(pid)
    filesets        = []
    fs = 1
    nevents = 0
    lo_run, lo_evt = 1000000, 1000000
    hi_run, hi_evt = -1, -1

    with open(file_catalog_fn, 'w') as fc:
        for line in lines:
            if line[0:1] != '#':
                p = line.split()
                nev   = to_int(field(p, 5))
                lorun = to_int(field(p, 6))
                loevt = to_int(field(p, 7))
                hirun = to_int(field(p, 8))
                hievt = to_int(field(p, 9))
                fc.write('%06d %s %10.3f  %s %s %6d %7d %9d %7d %9d\n' % (
                    fs, field(p, 1) or '', to_float(field(p, 2)),
                    field(p, 3) or '', field(p, 4) or '', nev, lorun, loevt, hirun, hievt))

                # fileset-level information
                nevents += nev
                if lorun < lo_run:
                    lo_run, lo_evt = lorun, loevt
                elif lorun == lo_run and loevt < lo_evt:
                    lo_evt = loevt
                if hirun > hi_run:
                    hi_run, hi_evt = hirun, hievt
                elif hirun == hi_run and hievt > hi_evt:
                    hi_evt = hievt

                if nevents > r['nev_fileset']:
                    filesets.append(Fileset(nevents, lo_run, lo_evt, hi_run, hi_evt))
                    nevents = 0
                    lo_run, lo_evt = 1000000, 1000000
                    hi_run, hi_evt = -1, -1
                    fs += 1
            else:
                fc.write(line + '\n')

    # last, not complete, fileset
    if hi_run != -1:
        filesets.append(Fileset(nevents, lo_run, lo_evt, hi_run, hi_evt))
    return file_catalog_fn, filesets


def write_fileset_catalog(r, filesets, pid):
    """the html-style fileset-level catalog"""
    fileset_catalog_fn = f'AAA_CATALOG.html.{pid}'
    comment_line = ('#' + '-' * 57 + '-' * 43 + '-' * 43)

    if not r['diskpool']:
        if r['access'] == 'xrootd':
            host = 'xroot://fndca1.fnal.gov'
        elif r['access'] == 'bluearc':
            host = 'bluearc://localhost.fnal.gov'
        else:
            host = r['host']
        directory = r['directory']
    else:
        host = r['diskpool']
        directory = f'/{r["book"]}/{r["dataset"]}'

    nev_total = sum(fs.nevents for fs in filesets)
    ldir = len(directory) if filesets else -1

    with open(fileset_catalog_fn, 'w') as fc:
        fc.write('<pre>\n')
        fc.write(comment_line + '\n')
        fc.write(f'# <b>dataset: {r["dataset"]}</b> \n')
        fc.write('# :DESCRIPTION: \n')
        fc.write('# fs       server                   directory   ')
        fc.write(' ' * max(ldir - 16, 0) + ' nevents\n')
        fc.write(comment_line + '\n')
        for i, fs in enumerate(filesets, start=1):
            fc.write('%06d  %-20s  %-45s %6d %7d %9d %7d %9d\n' % (
                i, host, directory, fs.nevents, fs.lorun, fs.loevt, fs.hirun, fs.hievt))

        width = ldir + 20
        fc.write(comment_line + '\n')
        fc.write(f'# :NTOTAL:          {nev_total:>{width}}\n')
        fc.write(f'# :NGENERATED:      {"":>{width}}\n')
        fc.write(f'# :MC_FLAG:         {1:>{width}}\n')
        fc.write(f'# :PROCESS_CODE:    {"":>{width}}\n')
        fc.write(f'# :PDG_CODE:        {"":>{width}}\n')
        fc.write(f'# :OFFLINE_VERSION: {r["offline_version"]:>{width}}\n')
        fc.write(comment_line + '\n')
        fc.write('</pre>\n')
    return fileset_catalog_fn


def install(r, fileset_catalog_fn, file_catalog_fn):
    cdir = os.path.join(r['install_dir'], r['book'], r['dataset'])
    os.makedirs(cdir, exist_ok=True)
    shutil.move(fileset_catalog_fn, os.path.join(cdir, 'AAA_CATALOG.html'))
    shutil.move(file_catalog_fn,    os.path.join(cdir, 'AAA_FILES.txt'))


def main(argv):
    a = parse_args(argv)
    if a.help or not argv:
        print('usage: catalog_stntuples.py [-v] [-u user] [-a access] -d dataset -H host '
              '-D directory [-b book] -p pattern [-m mode]')
        return 0

    r   = resolve_parameters(a)
    pid = os.getpid()
    if r['verbose']:
        print(f'host = {r["host"]}, hostname={r["local_host"]} directory={r["directory"]} '
              f'pattern={r["pattern"]}')

    request_file = f'aaa.request.{pid}'
    write_request_file(r, request_file)

    tmp_script, logfile = run_root(r, request_file, pid)

    file_catalog_fn, filesets = process_logfile(r, logfile, pid)
    fileset_catalog_fn = write_fileset_catalog(r, filesets, pid)

    if r['install']:
        install(r, fileset_catalog_fn, file_catalog_fn)

    if r['cleanup'] != 'no':
        for fn in (logfile, tmp_script, request_file):
            if os.path.exists(fn):
                os.remove(fn)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
