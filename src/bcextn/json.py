import json
import pathlib
import gzip

def load_json( f ):
    f = pathlib.Path(f)
    if f.name.endswith('.gz'):
        with gzip.open(f, 'r') as fh:
            return json.loads(fh.read().decode('utf8'))
    else:
        with f.open('rt') as fh:
            return json.load(fh)

def save_json( f, data, force = False ):
    f = pathlib.Path(f)
    if f.is_file():
        if not force:
            raise RuntimeError(f'File exists: {f}')
        else:
            f.unlink()
    assert not f.is_file()
    assert f.parent.is_dir()
    with f.open('wt') as fh:
        json.dump(data, fh, sort_keys=True)
    print(f"Wrote {f.name}")
