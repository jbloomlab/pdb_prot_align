"""
===============
nglview_struct
===============

Utilities for `nglview <http://nglviewer.org/nglview/latest/index.html>`_
structures.

"""


import pandas as pd

import pdb_prot_align.colorschemes


def color_by_site(cmd,
                  sites_df,
                  color_by,
                  pymol_object=None,
                  ):
    """Color residues in a ``pymol`` object by site.

    Parameters
    ----------
    cmd : pymol.cmd module
        Name of `pymol.cmd` module for current session, typically imported
        with ``import pymol.cmd as cmd``. The structure(s) must already be
        loaded into this session.
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
        str and are hex, then they need to be like this '0x25828e', not like
        this '#25828e'
    pymol_object : None or str
        Name of pymol object to which chains and sites belong. If `None`,
        color all residues in current session with specified chains / sites.

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

    # do the coloring
    if pymol_object:
        object_sel = f"{pymol_object} and "
    else:
        object_sel = ''
    for tup in sites_df.itertuples():
        chain = getattr(tup, chain_col)
        resi = getattr(tup, site_col)
        sel_str = f"{object_sel}chain {chain} and resi {resi}"
        color = getattr(tup, color_col).replace('#', '0x')
        cmd.color(color, sel_str)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
