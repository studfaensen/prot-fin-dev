"""
This module is to do some statistic analysis on the protfin results
"""

from actions import evaluate_protfin, select_samples, print_hash_counts, plot_frequencies
import argparse
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)  # fixes weird python error, look: https://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html


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
        prog='evaluation',
        description=__doc__,
    )
    parser.set_defaults(func=lambda _: parser.print_help())
    sub_commands = parser.add_subparsers()

    # evaluation.py eval <protfin-out-file>
    eval_parser = sub_commands.add_parser("eval", help="Evaluate protfin find-matches output")
    eval_parser.add_argument("protfin-out-file")
    eval_parser.set_defaults(func=lambda args: evaluate_protfin(getattr(args, "protfin-out-file")))

    # evaluation.py select-samples [-s] <mapman-file> <protein-file>
    eval_parser = sub_commands.add_parser("select-samples", help="Select samples from reference")
    eval_parser.add_argument("mapman-file")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("-s", "--samples-per-family", default=1, type=int)
    eval_parser.set_defaults(func=lambda args:
                             select_samples(
                                getattr(args, "mapman-file"),
                                getattr(args, "protein-file"),
                                args.samples_per_family
                             ))

    # evaluation.py print-hash-counts <database-file>
    eval_parser = sub_commands.add_parser("print-hash-counts", help="Print the calculated hash counts per sequence as csv")
    eval_parser.add_argument("database")
    eval_parser.set_defaults(func=lambda args: print_hash_counts(args.database))

    # evaluation.py print-hash-counts <database-file>
    eval_parser = sub_commands.add_parser("plot-frequencies", help="Plot the selected frequencies")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("out-file")
    eval_parser.set_defaults(func=lambda args: plot_frequencies(getattr(args, "protein-file"), getattr(args, "out-file")))

    return parser


if __name__ == '__main__':
    main()
