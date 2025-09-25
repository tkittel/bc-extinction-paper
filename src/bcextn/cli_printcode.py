import sys
from . import printcode

def main():
    args = sys.argv[1:]
    def usage():
        import os
        bn = os.path.basename(sys.argv[0])
        print("Usage:\n\n")
        print('  %s (%s) [OUTFILE]\n'%(bn,'|'.join(printcode.languages)))
        raise SystemExit(1)
    if len(args) not in (1,2) or args[0] not in printcode.languages:
        usage()
    outfile = args[1] if len(args)>1 else None
    if outfile is not None and (not outfile or not '.' in outfile or outfile.startswith('-')):
        usage()
    printcode.doit( lang=args[0],
                    outfile=outfile )
