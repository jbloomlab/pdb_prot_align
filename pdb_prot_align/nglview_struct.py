"""
===============
nglview_struct
===============

Utilities for `nglview <https://github.com/arose/nglview>`_ structures.

"""


import os
import re
import tempfile

import IPython.display

import nglview

import pandas as pd

import pdb_prot_align.colorschemes


def colorscheme_by_site(colorscheme_name,
                        sites_df,
                        color_by,
                        ):
    """Add a color scheme to the color registry.

    The scheme can then be added to a `ngview.widget.NGLWidget` as described
    `here <https://github.com/dwhswenson/contact_map/pull/62>`_. For instance::

        view.add_cartoon(color=colorscheme_name)


    Parameters
    ----------
    colorscheme_name : str
        Name of the color scheme.
    sites_df : pandas.DataFrame or str
        Information on how to color sites. Can either be data frame or name
        of CSV file with data frame. Must have columns named 'pdb_chain'
        and 'pdb_site' as well as the column specified by `color_by`.
    color_by : str or 2-tuple
        How to color the sites. Can either specify as a str the name of
        a column in `sites_df` that has the name of a color for each site,
        or can be the 2-tuple `(val_col, color_map)`. In this case,
        `val_col` is name of column with numerical values, and `color_map`
        is a :class:`pdb_prot_align.colorschemes.ValueToColorMap` that maps
        the numbers in this column to colors. If colors are specified as
        str and are hex, then they need to be like this '#25828e'.

    """
    site_col = 'pdb_site'
    chain_col = 'pdb_chain'

    if isinstance(sites_df, str):
        sites_df = pd.read_csv(sites_df)
    elif not isinstance(sites_df, pd.DataFrame):
        raise ValueError('`sites_df` must be data frame or name of CSV file')

    if isinstance(color_by, str):
        color_col = color_by
    elif (len(color_by) == 2 and isinstance(color_by[0], str) and
          isinstance(color_by[1],
                     pdb_prot_align.colorschemes.ValueToColorMap)):
        val_col = color_by[0]
        if val_col not in sites_df.columns:
            raise ValueError(f"`sites_df` lacks column {val_col}")
        color_map = color_by[1]
        color_col = 'color'
        if color_col in sites_df.columns:
            raise ValueError(f"`sites_df` can not have column {color_col} "
                             'if `color_by` is a 2-tuple')
        sites_df = (sites_df
                    .assign(
                        color=lambda x: x[val_col].map(color_map.val_to_color)
                        )
                    )

    cols = [site_col, chain_col, color_col]
    for col in cols:
        if col not in sites_df.columns:
            raise ValueError(f"`sites_df` lacks column {col}")

    # Drop duplicate and NaN rows, which could be case if data frame
    # is tidy and has amino acid identity, and ensure site is integer.
    sites_df = (
        sites_df
        [cols]
        .drop_duplicates()
        .dropna()
        .assign(**{site_col: lambda x: x[site_col].astype('int')})
        )
    # make sure just one color per site / chain
    dups = (sites_df
            .groupby([chain_col, site_col])
            .aggregate(nrows=pd.NamedAgg(color_col, 'count'))
            .query('nrows > 1')
            .reset_index()
            [[chain_col, site_col]]
            )
    if len(dups):
        raise ValueError('non-unique colors for some sites:\n' + str(dups))

    # Do the coloring; details on selection schemes:
    # https://github.com/arose/ngl/blob/master/doc/usage/selection-language.md
    colorscheme = []
    for tup in sites_df.itertuples():
        chain = getattr(tup, chain_col)
        resi = getattr(tup, site_col)
        if isinstance(resi, float):
            if resi != int(resi):
                raise ValueError(f"non-integer residue {resi}")
            resi = int(resi)
        sel_str = f":{chain} and {resi}"
        color = getattr(tup, color_col)
        colorscheme.append([color, sel_str])

    nglview.color.ColormakerRegistry.add_scheme(colorscheme_name, colorscheme)


def render_html(view,
                *,
                html_file=None,
                orientation=None,
                remove_widget_view=False,
                ):
    """Render widget to HTML display.

    Parameters
    ----------
    view : ngview.widget.NGLWidget
        The structure widget to render.
    html_file : None or str
        If you want to also save to permanent HTML file, provide name here.
    orientation : list
        Set to this camera orientation (list of 16 numbers), fixing this bug:
        https://github.com/dwhswenson/contact_map/pull/62#issuecomment-583788933
        You can get the desired orientation by manually manipulating the widget
        in a Jupyter notebook and then calling `view._camera_orientation`.
    remove_widget_view : bool
        Remove the widget view lines, so the HTML just gives the widget state.
        Helpful if you want to embed widgets in HTML rendering without showing
        another time.

    Returns
    -------
    IPython.display.DisplayHandle
        A handle that displays the HTML structure in a Jupyter notebook.

    """
    with tempfile.TemporaryFile(mode='w+', suffix='.html') as f:
        nglview.write_html(f, view)
        f.flush()
        f.seek(0)
        html_text = f.read()

    if orientation:
        if len(orientation) != 16:
            raise ValueError('`orientation` must be list of 16 numbers')
        html_text = re.sub(r'"_camera_orientation":\s+\[[^\]]*\]',
                           '"_camera_orientation": [' +
                           ', '.join(map(str, orientation)) + ']',
                           html_text)

    if remove_widget_view:
        widget_view_regex = (
                r'<script type="application/vnd\.jupyter\.widget\-view\+json">'
                r'\n.*?\n'
                '</script>'
                )
        html_text = re.sub(widget_view_regex, '', html_text)

    if html_file:
        if os.path.splitext(html_file)[1] != '.html':
            raise ValueError(f"`html_file` needs extension .html: {html_file}")
        with open(html_file, 'w') as f_html:
            f_html.write(html_text)

    return IPython.display.display(IPython.display.HTML(data=html_text))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
