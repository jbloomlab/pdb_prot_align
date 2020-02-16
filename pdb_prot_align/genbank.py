"""
=======
genbank
=======

Extract information from Genbank-format files.

"""


import collections
import os

import Bio.SeqFeature  # noqa: F401
import Bio.SeqIO
import Bio.SeqRecord  # noqa: F401

import pandas as pd


def genbank_to_feature_df(gb,
                          multi_features,
                          *,
                          feature_color_map=None,
                          featureless_sites='NaN',
                          ):
    """Convert Genbank with annotated features to map of sites to features.

    Parameters
    ----------
    gb : BioPython.SeqRecord.SeqRecord or name of Genbank file
        Contains features annotated by their `type`.
    multi_features : {'last', 'first', 'all', 'error'}
        What to do if a site is in multiple features? Assign it to 'last' or
        'first' feature (as ordered in `gb`) that it is part of, assign it
        to 'all' features it is part of, or raise an 'error'.
    feature_color_map : None or dict
        If not `None`, add a column mappping each feature to a color. The
        dict must have a key for each feature type in `gb`.
    featureless_sites : {'NaN', 'drop'}
        How to handle sites without features. Keep them with features of
        `NaN`, or drop them from returned data frame.

    Returns
    -------
    pandas.DataFrame
        Columns are 'isite' (1, ... numbering), 'amino_acid', 'feature', and
        optionally 'color' if `feature_color_map` is not `None`. If mapping
        sites to multiple features, the returned data frame is tidy.

    Example
    -------
    >>> gb = Bio.SeqRecord.SeqRecord(
    ...         seq='MGKLLIT',
    ...         id='example',
    ...         name='example',
    ...         description='example',
    ...         features=[Bio.SeqFeature.SeqFeature(
    ...                     type='1to6',
    ...                     location=Bio.SeqFeature.FeatureLocation(0, 6)),
    ...                   Bio.SeqFeature.SeqFeature(
    ...                     type='1to2',
    ...                     location=Bio.SeqFeature.FeatureLocation(0, 2)),
    ...                   Bio.SeqFeature.SeqFeature(
    ...                     type='4to6',
    ...                     location=Bio.SeqFeature.FeatureLocation(3, 6)),
    ...                   ],
    ...         )

    >>> genbank_to_feature_df(gb, 'last')
       isite amino_acid feature
    0      1          M    1to2
    1      2          G    1to2
    2      3          K    1to6
    3      4          L    4to6
    4      5          L    4to6
    5      6          I    4to6
    6      7          T     NaN

    >>> genbank_to_feature_df(gb, 'last', featureless_sites='drop')
       isite amino_acid feature
    0      1          M    1to2
    1      2          G    1to2
    2      3          K    1to6
    3      4          L    4to6
    4      5          L    4to6
    5      6          I    4to6

    >>> genbank_to_feature_df(gb, 'first')
       isite amino_acid feature
    0      1          M    1to6
    1      2          G    1to6
    2      3          K    1to6
    3      4          L    1to6
    4      5          L    1to6
    5      6          I    1to6
    6      7          T     NaN

    >>> genbank_to_feature_df(gb, 'all')
        isite amino_acid feature
    0       1          M    1to6
    1       1          M    1to2
    2       2          G    1to6
    3       2          G    1to2
    4       3          K    1to6
    5       4          L    1to6
    6       4          L    4to6
    7       5          L    1to6
    8       5          L    4to6
    9       6          I    1to6
    10      6          I    4to6
    11      7          T     NaN

    >>> genbank_to_feature_df(gb, 'error')
    Traceback (most recent call last):
        ...
    ValueError: site 1 in multiple features

    >>> feature_color_map = {'1to6': 'red', '1to2': 'blue', '4to6': 'green'}
    >>> genbank_to_feature_df(gb, 'last', feature_color_map=feature_color_map)
       isite amino_acid feature  color
    0      1          M    1to2   blue
    1      2          G    1to2   blue
    2      3          K    1to6    red
    3      4          L    4to6  green
    4      5          L    4to6  green
    5      6          I    4to6  green
    6      7          T     NaN    NaN

    """
    if isinstance(gb, str):
        if not os.path.isfile(gb):
            raise IOError(f"no `gb` file:\n{gb}")
        gb = Bio.SeqIO.read(gb, 'genbank')
    elif not isinstance(gb, Bio.SeqRecord.SeqRecord):
        raise ValueError(f"`gb` not str or SeqRecord:\n{gb}")

    if multi_features == 'last':
        featurelist = list(reversed(gb.features))
    else:
        featurelist = gb.features

    site = []  # 0-indexed sites
    feature = []
    nfeatures_for_site = collections.defaultdict(int)
    for r in range(len(gb)):
        for feat in featurelist:
            if r in feat:
                site.append(r)
                feature.append(feat.type)
                if multi_features in {'last', 'first'}:
                    break
                elif multi_features == 'error':
                    nfeatures_for_site[r] += 1
                    if nfeatures_for_site[r] > 1:
                        raise ValueError(f"site {r + 1} in multiple features")
                elif multi_features != 'all':
                    raise ValueError(f"bad `multi_features` {multi_features}")

    if featureless_sites == 'NaN':
        site_set = set(site)
        for r in range(len(gb)):
            if r not in site_set:
                site.append(r)
                feature.append(float('nan'))
    elif featureless_sites != 'drop':
        raise ValueError(f"invalid `featureless_sites` {featureless_sites}")

    df = (pd.DataFrame({'site': site,
                        'feature': feature,
                        })
          .assign(isite=lambda x: x['site'] + 1,
                  amino_acid=lambda x: x['site'].map(dict(enumerate(gb.seq))),
                  )
          [['isite', 'amino_acid', 'feature']]
          .sort_values('isite')
          )

    if feature_color_map is not None:
        missing = [feat.type for feat in featurelist
                   if feat.type not in feature_color_map]
        if missing:
            raise ValueError(f"`feature_color_map` lacks features:\n{missing}")
        df = df.assign(color=lambda x: x['feature'].map(feature_color_map))

    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
