
import sys
from . import plot

def modes():
    return sorted( m[5:] for m in dir(plot)
                   if m.startswith('main_') )

def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit('Please provide plot mode. One of "%s"'%('", " '.join(modes())))
    mode = args[0]
    args = args[1:]
    if mode=='list':
        print( '\n'.join(modes() ) + '\n' )
    fctname = f'main_{mode}'
    fct = getattr(plot,fctname,None)
    if not fct:
        raise SystemExit('Invalid mode: "%s"'%mode)
    fct( args )
