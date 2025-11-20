#!/usr/bin/env python

import pathlib

srcdir = pathlib.Path(__file__).parent
outdir = srcdir.joinpath('submit')
maintex = srcdir.joinpath('bcextnupdate.tex')

def main():
    assert not outdir.exists()
    outdir.mkdir()
    main_merged = []
    ifig = 1
    for line in maintex.read_text().splitlines():
        if line.startswith(r'\input{'):
            n=line[len(r'\input{'):].strip()
            assert n.endswith('}')
            n = n[:-1]
            f = srcdir.joinpath(*n.split('/'))
            print(f"Including {n}")
            for e in f.read_text().splitlines():
                main_merged.append(e.rstrip())
        elif r'\includegraphics' in line:
            assert line.count('{')==1
            assert line.count('}')==1
            n = line.split('{')[1].split('}')[0]
            assert n.startswith('generated')
            assert '{%s}'%n in line
            f = srcdir.joinpath(*n.split('/'))
            assert f.name.endswith('.pdf')
            assert f.exists()
            print(f"Using graphics {n}")
            newname = 'fig_%02i_%s'%(ifig,f.name)
            fo = outdir.joinpath(newname)
            fo.write_bytes(f.read_bytes())
            ifig += 1
            main_merged.append(line.replace('{%s}'%n,'{%s}'%fo.name))
        else:
            main_merged.append(line.rstrip())

    main_merged.append('')
    outdir.joinpath('paper.tex').write_text('\n'.join(main_merged))
    for n in ['iucrjournals.cls','harvard.sty',
              'iucr.bst','iucr.bib']:
        f = srcdir.joinpath(n)
        assert f.is_file()
        outdir.joinpath(f.name).write_bytes(f.read_bytes())



if __name__=='__main__':
    main()
