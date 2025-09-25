import sys
from . import printcode

def main():
    args = sys.argv[1:]
    if len(args)!=1 or args[0] not in printcode.languages:
        raise SystemExit("ERROR: Please specify language (one of %s)"%(' '.join(printcode.languages)))
    lang = args[0]
    printcode.doit(lang)
