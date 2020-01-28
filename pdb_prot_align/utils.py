"""Utilities for aligning proteins to PDB."""


import re
import subprocess
import tempfile
import warnings

import Bio.AlignIO
import Bio.Data.IUPACData
import Bio.PDB
import Bio.SeqIO
import Bio.SeqUtils

import natsort

import pandas as pd


def align_prots_mafft(prots, *, mafft='mafft'):
    """Align protein sequences with ``mafft``.

    Parameters
    ----------
    prots : list
        List of sequences as `Biopython.SeqRecord.SeqRecord` objects.
    alignmentfile : str
        Name of created FASTA alignment.
    mafft : str
        Command that resolves to ``mafft`` executable.

    Returns
    -------
    Bio.Align.MultipleSeqAlignment
        The protein alignment.

    """
    for p in prots:
        if not re.fullmatch(
                f"[{Bio.Data.IUPACData.extended_protein_letters}]+",
                str(p.seq)):
            raise ValueError(f"invalid amino acids for {p.description}\n" +
                             str(p.seq))

    try:
        _ = subprocess.run([mafft, '--version'],
                           check=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL,
                           )
    except FileNotFoundError:
        raise ValueError(f"cannot call `mafft` executable with: {mafft}")

    with tempfile.NamedTemporaryFile('w') as prots_in:
        Bio.SeqIO.write(prots, prots_in, 'fasta')
        prots_in.flush()
        res = subprocess.run([mafft, '--amino', prots_in.name],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL,
                             universal_newlines=True,
                             )
        alignment_text = res.stdout

    with tempfile.TemporaryFile('w+') as f:
        f.write(alignment_text)
        f.flush()
        f.seek(0)
        aln = Bio.AlignIO.read(f, 'fasta')

    return aln


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
    with warnings.catch_warnings():
        warnings.simplefilter(
                'ignore',
                category=Bio.PDB.PDBExceptions.PDBConstructionWarning)
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
