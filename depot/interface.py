#!/usr/bin/env python3

"""
Usage:
  R2V2 [<module>] [<args>...] [-h] [-v]

Modules:
  prepare         prepare libraries for otu clustering
  dereplicate     dereplicate and perform orf filtering
  cluster         cluster raw sequence data into OTUs

Options
  -h, --help      show this
  -v, --version   show version number

  See 'R2V2 <command> --help' for more information on a specific command.

"""

import sys
from docopt import docopt
from depot import prepare, dereplicate, cluster

def main():
    """Modules"""
    args = docopt(__doc__,version='1.0.0', options_first=True)
    if args['<module>'] == 'prepare':
        prepare.main()
    elif args['<module>'] == 'dereplicate':
        dereplicate.main()
    elif args['<module>'] == 'cluster':
        cluster.main()
    else:
        sys.exit(f"{args['<module>']} is not an R2V2 module. See 'R2V2 -h'.")
