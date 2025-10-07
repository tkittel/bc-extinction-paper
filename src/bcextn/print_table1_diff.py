
#FIXME rename file

import pathlib
import numpy as np

def is_html( out_filename ):
    is_html = pathlib.Path(out_filename).name.endswith('.html')
    is_tex = pathlib.Path(out_filename).name.endswith('.tex')
    assert int(is_html) + int(is_tex) == 1
    return is_html

def fix_table_fmt( table, out_filename ):
    table.columns = ['%g'%e for e in table.columns]
    table.index = ['%g'%e for e in table.index]
    if is_html(out_filename):
        table.columns.name = 'sin(\u03B8)'
        table.index.name = 'x'
    else:
        table.columns.name = r'$\sin(\theta)$'
        table.index.name = '$x$'

def table_output( styled_table, out_filename, table_id ):
    assert '"' not in table_id
    #styled_table = styled_table.set_table_attributes('id="%s"'%table_id)
    if is_html(out_filename):
        output = styled_table.to_html(uuid=table_id)
    else:
        output = styled_table.to_latex(convert_css=True)
    with pathlib.Path(out_filename).open('wt') as fh:
        fh.write(output)
        print(f"Wrote: {out_filename}")

def write_updated_table1( table1_updated,
                          out_filename,
                          *,
                          tablename,
                          do_print = False,
                          ndigits = 5 ):
    t = table1_updated
    fix_table_fmt(t,out_filename)
    #Follow style of BC1974 Table 1, multiply by 1e{ndigits} and show now decimals:
    t *= 10**ndigits
    t = table1_updated.style.format("{:0%i.0f}"%ndigits)
    if do_print:
        print(t)#fixme does not work
    table_id = 'Upd'+tablename.capitalize()
    if ndigits!=4:
        table_id+='Digits%i'%ndigits
    table_output(t,out_filename,table_id=table_id)

def write_table1_diff_heatmap( table_orig,
                               table_new,
                               out_filename,
                               *,
                               tablename,
                               do_print = False,
                               norm_is_min_y_1minusy = True ):
    if not norm_is_min_y_1minusy:
        table_diff = (table_orig/table_new - 1.0)*100.0
    else:
        table_diff = (table_orig-table_new)*100.0 / np.minimum(table_new.values, 1 - table_new.values)
    vmin, vmax = table_diff.min().min(), table_diff.max().max()
    fix_table_fmt(table_diff,out_filename)
    if do_print:
        print(table_diff)
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
    table_id = 'Diff'+tablename.capitalize()
    if not norm_is_min_y_1minusy:
        table_id+='RawNorm'
    table_output(td,out_filename,table_id =table_id)

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
