#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ported by nickp60 on Pi Day 2017
# version 0.1
import sys
import os
import platform
import argparse
import subprocess
import operator
import shutil
import logging
import pkg_resources


# global variables
OPSYS = platform.system().lower()
EXE = __file__
VERSION = "0.0.4"
DESC = "rapid ribosomal RNA prediction"
AUTHOR = 'Torsten Seemann <torsten.seemann@gmail.com>'
URL = 'https://github.com/Victorian-Bioinformatics-Consortium/barrnap'

#DBDIR = os.path.join(os.path.dirname(os.path.dirname(EXE)),"db")
# our db of HMMs should be installed in site-packages somewhere when the
# package was installed.  we set this below, to avoid issues when testing
resource_package = pkg_resources.Requirement.parse("barrnap")
try:
    DBDIR = pkg_resources.resource_filename(
        resource_package, 'barrnap/db')
except Exception as e:
    DBDIR = ""

NHMMER = shutil.which("nhmmer")

LENG = {
    "5S_rRNA": 119,
    "16S_rRNA": 1585,
    "23S_rRNA": 3232,
    "5_8S_rRNA": 156,
    "18S_rRNA": 1869,
    "28S_rRNA": 2912,
    "12S_rRNA": 954}
MAXLEN = {k: 1.2 * v for (k, v) in LENG.items()}


def showCitation():
    return("""
If you use Barrnap in your work, please cite:
    Seemann T (2013)
    {0} {1} : {2}
    {3}

Thank you.
    """.format(os.path.basename(EXE), VERSION, DESC, URL))


def showUsage():
    return(str("{0} [options] <chromosomes.fasta>\n".format(EXE) +
           "Synopsis:\n  {0} {1} - {2}\n".format(EXE, VERSION, DESC)) +
           "Author:\n  {0}\n".format(AUTHOR))


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(
        description="barrnap ported to python3",
        usage=showUsage(),
        add_help=False)  # to allow for custom help
    parser.add_argument("fasta", action="store")

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-k", "--kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-t", "--threads",
                          action="store", type=float,
                          help="Number of threads/cores/CPUs to use;" +
                          "default: %(default)s", default=8,
                          )
    optional.add_argument("-e", "--evalue",
                          default=.000001,
                          help="Similarity e-value cut-off; " +
                          "default: %(default)s",
                          action="store", type=float)
    optional.add_argument("-l", "--lencutoff",
                          default=0.8,
                          help="Proportional length threshold to label as " +
                          "partial; default: %(default)s",
                          action="store", type=float)
    optional.add_argument("-r", "--reject",
                          help="Proportional length threshold to reject " +
                          "prediction; default: %(default)s",
                          default=0.5,
                          action="store", type=float)
    optional.add_argument("-i", "--incseq",
                          default=False,
                          help="Include FASTA input sequences in GFF3 output",
                          action="store_true")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="This help")
    optional.add_argument("-v", "--version",
                          action="version", version=VERSION,
                          help="Print version and exit")
    optional.add_argument("--citation",
                          action="version", version=showCitation(),
                          help="Print citation for referencing barrnap")
    args = parser.parse_args()
    return args


def main(args, logger=None):
    if logger is None:
        logger = logging.getLogger('barrnap')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    if args is None:
        args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    logger.info("This is {0} {1}".format(EXE, VERSION))
    logger.info("Written by {0}".format(AUTHOR))
    logger.info("Obtained from {0}".format(URL))

    logger.debug("Detected operating system: {0}".format(platform.system()))
    # if OPSYS not in os.listdir(os.path.join(os.path.dirname(os.path.dirname(EXE)),"binaries")):
    #     if shutil.which(NHMMER) is None:
    #         logger.error("No binary for your OS '{0}' is included. If you have one, copy it to $NHMMER. Exiting...".format(OPSYS))
    #         sys.exit(1)
    #     else:
    #         NHMMER = shutil.which(nhmmer)
    # else:
    #     pass
    #  check threads
    if not args.threads > 0:
        logger.error("Invalid --threads $threads")
        sys.exit(1)
    else:
        logger.debug("Will use {0} threads".format(args.threads))

    #  check lencutoff
    if not args.lencutoff > 0:
        logger.error("Invalid --lencutoff {0}")
        sys.exit(1)
    else:
        logger.debug("Will tag genes  < {0} of expected length ".format(args.lencutoff))

    #  check evalue
    if not args.evalue > 0:
        logger.error("Invalid --evalue {0}")
        sys.exit(1)
    else:
        logger.debug("Setting evalue cutoff to {0}".format(args.evalue))

    #  check reject
    if not args.reject > 0:
        logger.error("Invalid --reject cutoff")
        sys.exit(1)
    else:
        logger.debug("Will reject genes < {0} of expected length.".format(args.reject))
    # no need to check kingdom, argparse handles that
    hmmdb = os.path.join(DBDIR, "{0}.hmm".format(args.kingdom))
    if not os.path.exists(hmmdb):
        logger.warning("Can't find database: {0}".format(hmmdb))
        logger.warning(
            "attempting to find databases in bin/db, in case " +
            " this was installed via conda alongside regular" +
            " barrnap")
        hmmdb = os.path.join(
            os.path.dirname(__file__),
            "db",
            "{0}.hmm".format(args.kingdom))
        if not os.path.exists(hmmdb):
            logger.error("Nope, can't find database: {0}".format(hmmdb))
            sys.exit(1)
    else:
        logger.debug("Using database: {0}".format(hmmdb))


    logger.info("Scanning {0} for {1} rRNA genes... please wait".format(args.fasta, args.kingdom))
    cmd = str("{0} --cpu {1} -E {2} --w_length {3}  -o /dev/null --tblout /dev/stdout {4} {5}").format(
        NHMMER,
        args.threads,
        args.evalue,
        int(max(MAXLEN.values())),  # must be an int
        hmmdb,
        args.fasta)
    logger.debug(cmd)
    hmmer_results = subprocess.run(cmd, shell=sys.platform != "win32",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, check=True)
    results_list = [x for x in
                    hmmer_results.stdout.decode("utf-8").split("\n")]
    feats = []
    for idx, line in enumerate(results_list):
        if line.startswith("#") or len(line.replace(" ", "")) == 0:
            logger.debug("skipping line %i" % idx)
            continue
        x = [x for x in line.split(" ") if len(x) != 0]
        if int(x[6]) < int(x[7]):
            begin, end, strand = int(x[6]), int(x[7]), '+'
        else:
            begin, end, strand = int(x[7]), int(x[6]), '-'
        seqid, gene, prod = x[0], x[2], x[2]
        try:
            score = x[12]
        except:
            score = "."
        #  check length against refs
        if gene not in LENG.keys():
            logger.error("Detected unknown gene '$gene' in scan, aborting.")
            sys.exit(1)
        # assumes NAME field in .HMM file is of form "16s_rRNA" etc
        prod = prod.replace("_r", " ribosomal ")  # convert the short ID to an English string
        prod = prod.replace("5_8", "5.8")  # correct naming for 5.8S

        # check for incomplete alignments
        note = ''
        length = end - begin + 1
        if length < int(args.reject * LENG[gene]):
            logger.info("Rejecting short {0} nt predicted $gene. Adjust via --reject option.".format(length))
            continue
        elif length < int(args.lencutoff * LENG[gene]):
            note = ";note=aligned only {0} percent of the $prod".format(
                100 * length / LENG[gene])  # sic
            prod = prod + " (partial)"

        logger.info("Found: {0} {1} L={2} {3}..{4} {5} {6}".format(
            gene,
            seqid,
            length / LENG[gene],
            begin,
            end,
            strand,
            prod))
        tags = "Name={0};product={1}{2}".format(gene, prod, note)
        feats.append([seqid, "{0}:{1}".format(EXE, VERSION),
                      begin, end, score, strand, tags])

    logger.info("Found %i ribosomal RNA features." % len(feats))
    logger.info("Sorting features and outputting GFF3...")
    sys.stdout.write("##gff-version 3\n")
    seqids = set([i[0] for i in feats])
    for seqid in seqids:
        for feat in sorted(feats, key=operator.itemgetter(2)):
            if feat[0] != seqid:
                continue
            sys.stdout.write("\t".join([str(j) for j in feat]) + "\n")
    if args.incseq:
        sys.stdout.write("##FASTA\n")
        with open(args.fasta, "r") as fa:
            for line in fa:
                sys.stdout.write(line)

    sys.exit(0)


if __name__ == "__main__":
    args = get_args()
    main(args)
