# python -m chibreak <run> [--save]
import sys

from .experiments import main

sys.exit(main(sys.argv[1:]))
