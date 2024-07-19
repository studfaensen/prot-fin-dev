#! python3
"""
A tool to find a protein and its relatives based on chemical-physical
features of their amino acid sequences
"""

import argparse
from tools import DBConfig
from actions import create_db, find_matches, match_family
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)  # fixes weird python error, look: https://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html

DB_DEFAULT = "database.pickle"


def main():
    """
    This is where the command line interface is defined to interpret
    the passed arguments to execute the respective functions.
    """

    cli = get_cli()

    args = cli.parse_args()
    args.func(args)


def get_cli():
    parser = argparse.ArgumentParser(
        prog='prot-fin',
        description=__doc__,
    )
    parser.set_defaults(func=lambda _: parser.print_help())
    sub_commands = parser.add_subparsers()

    # protfin.py create-db [-p] <fasta-file>
    create_db_parser = sub_commands.add_parser("create-db", help="Create Database")
    create_db_parser.add_argument("fasta-file")
    create_db_parser.add_argument("-p", "--path", default=DB_DEFAULT)
    create_db_parser.add_argument("-c", "--cpu", default=1, type=int)
    create_db_parser, dbconfig = cli_dbconfig(create_db_parser)
    create_db_parser.set_defaults(func=lambda args:
                                  create_db(
                                      getattr(args, "fasta-file"),
                                      db_out=args.path,
                                      cpu_count=args.cpu,
                                      **dbconfig(args)
                                  ))

    # protfin.py find-matches [-d] <fasta-file>
    find_match_parser = sub_commands.add_parser("find-matches", help="Find Matches for Proteins")
    find_match_parser.add_argument("fasta-file")
    find_match_parser.add_argument("-d", "--database", default=DB_DEFAULT)
    find_match_parser.add_argument("-f", "--filter", default=1., type=float)
    find_match_parser.set_defaults(func=lambda args:
                                   find_matches(
                                       getattr(args, "fasta-file"),
                                       db_in=args.database,
                                       filter_quantile=args.filter
                                   ))

    # protfin.py match-family [-d] <fasta-file>
    find_match_parser = sub_commands.add_parser("match-family", help="Find Matches for Proteins")
    find_match_parser.add_argument("family-file")
    find_match_parser.add_argument("-d", "--database", default=DB_DEFAULT)
    find_match_parser.add_argument("-f", "--filter", default=1., type=float)
    find_match_parser.set_defaults(func=lambda args:
                                   match_family(
                                       getattr(args, "family-file"),
                                       db_in=args.database,
                                       filter_quantile=args.filter
                                   ))

    return parser


def cli_dbconfig(parser: argparse.ArgumentParser):
    default_config = DBConfig()
    parser.add_argument("-w", "--window-size", default=default_config.window_size, type=int)
    parser.add_argument("-o", "--overlap", default=default_config.overlap, type=int)
    parser.add_argument("-n", "--npeaks", default=default_config.n_peaks, type=int)
    parser.add_argument("-s", "--significance", default=default_config.significance, type=float)
    parser.add_argument("-m", "--selection-method", default=default_config.selection_method, choices=("none", "absolute", "deviation"), type=str)
    parser.add_argument("-k", "--skip-first-k-freqs", default=default_config.skip_first_k_freqs, type=int)
    return parser, lambda args: {
        "window_size": args.window_size,
        "overlap": args.overlap,
        "n_peaks": args.npeaks,
        "selection_method": args.selection_method,
        "significance": args.significance,
        "skip_first_k_freqs": args.skip_first_k_freqs
    }


if __name__ == '__main__':
    main()
