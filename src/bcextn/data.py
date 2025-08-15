
def datadir():
    import pathlib
    dd= pathlib.Path(__file__).parent.parent.parent.joinpath('data')
    assert dd.is_dir()
    return dd

def data_file( name, must_exist = True ):
    f = datadir().joinpath(name)
    if must_exist and not f.is_file():
        raise RuntimeError(f'File not found: {f}')
    return f

def load_json_data( name ):
    from .json import load_json
    return load_json(data_file(name))
