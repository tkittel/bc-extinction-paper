import pathlib
import numpy as np

def format_table1_heatmap( table_orig, table_new, out_filename, mode='html',do_print = False ):
    is_html = out_filename.endswith('.html')
    is_tex = out_filename.endswith('.tex')
    assert int(is_html) + int(is_tex) == 1
    table_diff = (table_orig/table_new - 1.0)*100.0
    if do_print:
        print(table_diff)
    vmin, vmax = table_diff.min().min(), table_diff.max().max()
    table_diff.columns = ['%g'%e for e in table_diff.columns]
    table_diff.index = ['%g'%e for e in table_diff.index]

    print( table_diff )

    def my_scale( val ):
        ##Any value diverging less than 1% we clip to neutral. Otherwise we use a colourbar with extreme values
        #if abs(val)<1.0:
        #    return 0.5
        x = 0.5 + 0.5 * val/table_diff.max().max()
        assert -1 <= x <= 1
        return x
    cmap = create_custom_cmap( table_diff,
                               scale_fct = my_scale,
                               orig_name = 'custom' )
    td = table_diff.style.background_gradient(cmap=cmap,vmin = vmin,vmax = vmax)
    td = td.format("{:.2f}")

    if is_html:
        output = td.to_html()
    else:
        output = td.to_latex(convert_css=True)

    with pathlib.Path(out_filename).open('wt') as fh:
        fh.write(output)
        print(f"Wrote: {out_filename}")

def create_custom_cmap( table,
                        scale_fct,
                        orig_name = 'coolwarm' ):
    import matplotlib.pyplot as plt
    if orig_name == 'custom':
        import matplotlib.colors as mcol
        orig_cmap = mcol.LinearSegmentedColormap.from_list( name='customcoolwarm',
                                                            colors =[(0, 0, 1),
                                                                     (1, 1., 1),
                                                                     (1, 0, 0)],
                                                            N=256,gamma=1.0 )
    else:
        orig_cmap = plt.get_cmap(orig_name)
    amin, amax = float(table.min().min()), float(table.max().max())
    def x_to_unscaled_a( x ):
        assert -0.000001 < x < 1.000001
        return amin + x * ( amax-amin )
    def new_cmap(x):
        a = x_to_unscaled_a( x )
        newx = scale_fct( a )
        return orig_cmap(newx)
    from matplotlib.colors import Colormap

    class MyCMap(Colormap):
        def __init__(self, name='custom_cmap'):
            super().__init__(name)
        def __call__(self, x, alpha=1.0, **kwargs):
            assert alpha==1.0
            x = np.clip(x, 0.0, 1.0)
            return np.array([new_cmap(val) for val in x])
    return MyCMap()
