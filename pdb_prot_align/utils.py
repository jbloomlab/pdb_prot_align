"""Utilities for aligning proteins to PDB."""


import Bio.PDB
import Bio.SeqUtils

import natsort

import pandas as pd


def pdb_seq_to_number(pdbfile, chain_ids, chain_identity='union'):
    """Get sequence and number of chain(s) for some protein in PDB file.

    Parameters
    ----------
    pdbfile : str
        Path to existing PDB file.
    chain_ids : list
        List of chains in PDB file. All these chains must correspond to the
        same underlying molecule (i.e., be monomers in homo-oligomer).
    chain_identity : {'union', 'intersection', 'require_same'}
        How to parse chains. They are required to share the same wildtype
        at all sites they have in common. If all sites are not shared
        between all chains, take the union, the intersection, or raise
        an error if they are not exactly the same.

    Returns
    -------
    pandas.DataFrame
        Columns are:
            - 'sequential' : sequential 1, 2, ... numbering of residues
            - 'pdb_site' : PDB file residue number
            - 'wildtype' : wildtype residue identity (1-letter code)

    """
    structure = Bio.PDB.PDBParser().get_structure('_', pdbfile)
    if len(structure) != 1:
        raise ValueError(f"{pdbfile} has multiple models")
    else:
        model = list(structure)[0]

    aa_3to1 = {three.upper(): one for three, one in
               Bio.SeqUtils.IUPACData.protein_letters_3to1.items()}

    df = pd.DataFrame({}, columns=['pdb_site', 'wildtype', 'chain'])
    for chain_id in chain_ids:
        chain = model[chain_id]
        pdb_sites = []
        wildtypes = []
        for res in chain.child_list:
            heteroflag, resnum, insertcode = res.get_id()
            if not heteroflag.strip():  # ignore hetero-residues
                pdb_sites.append(f"{resnum}{insertcode.strip()}")
                wildtypes.append(aa_3to1[res.resname])
        if pdb_sites != natsort.natsorted(pdb_sites):
            raise ValueError(f"residues in {pdbfile} chain {chain_id} not "
                             'naturally sorted')
        df = df.append(pd.DataFrame({'pdb_site': pdb_sites,
                                     'wildtype': wildtypes,
                                     'chain': chain_id}))

    # check all chains have same wildtype at each position
    mismatched_wt = (
        df
        .groupby('pdb_site', sort=False)
        .aggregate({'wildtype': 'nunique'})
        .reset_index()
        .query('wildtype != 1')
        ['pdb_site']
        .tolist()
        )
    if mismatched_wt:
        raise ValueError(f"{pdbfile} chains {', '.join(chain_ids)} differ in "
                         f"wildtype at sites:\n{', '.join(mismatched_wt)}")

    if chain_identity == 'require_same':
        not_exact = (
            df
            .groupby(['pdb_site', 'wildtype'], sort=False)
            .aggregate({'chain': 'count'})
            .reset_index()
            .query(f"chain != {len(chain_ids)}")
            ['pdb_site']
            .tolist()
            )
        if not_exact:
            raise ValueError(f"`chain_identity` is {chain_identity}, but "
                             f"{pdbfile} chains {', '.join(chain_ids)} "
                             f"don't all have sites:\n{', '.join(not_exact)}")
        else:
            df = df[['pdb_site', 'wildtype']].drop_duplicates()
    elif chain_identity == 'union':
        df = (df
              [['pdb_site', 'wildtype']]
              .drop_duplicates()
              )
    elif chain_identity == 'intersection':
        df = (df
              .groupby(['pdb_site', 'wildtype'], sort=False)
              .aggregate({'chain': 'count'})
              .reset_index()
              .query(f"chain == {len(chain_ids)}")
              )
    else:
        raise ValueError(f"invalid `chain_identity` {chain_identity}")

    # natural sort on site as here: https://stackoverflow.com/a/29582718
    df = df.set_index('pdb_site')
    df = (df
          .reindex(index=natsort.natsorted(df.index))
          .reset_index()
          .assign(sequential=lambda x: x.reset_index().index + 1)
          [['sequential', 'pdb_site', 'wildtype']]
          )
    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
