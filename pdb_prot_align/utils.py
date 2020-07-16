"""
=====
utils
=====

Utilities for aligning proteins to PDB.

"""


import collections
import re
import subprocess
import tempfile
import textwrap  # noqa: F401
import warnings

import Bio.AlignIO
import Bio.Data.IUPACData
import Bio.PDB
import Bio.SeqIO
import Bio.SeqUtils

import natsort

import numpy

import pandas as pd


def alignment_to_count_df(fasta_alignment,
                          *,
                          ignore_gaps=False,
                          as_freqs=False,
                          add_entropy_neff=False,
                          ):
    """Convert FASTA alignment to counts data frame.

    Parameters
    ----------
    fasta_alignment : str
        FASTA file with alignment.
    ignore_gaps : bool
        Ignore gap characters ('-') in counts.
    as_freqs : bool
        Return frequencies rather than counts.
    add_entropy_neff : bool
        Add columns with site entropy (bits) and number effective amino acids.
        Requires `as_freqs` to be `True`.

    Returns
    -------
    pandas.DataFrame
        First column ('isite') is 1, 2, ... site numbering. Other columns give
        counts of each observed alignment character at that site.

    Example
    --------
    >>> with tempfile.NamedTemporaryFile('w+') as f:
    ...     _ = f.write(textwrap.dedent(
    ...             '''
    ...             >seq1
    ...             GK-AC-L
    ...             >seq2
    ...             -K-AM-L
    ...             >seq3
    ...             G-NAC-L
    ...             >seq4
    ...             GK-ACSI
    ...             '''
    ...             ))
    ...     f.flush()
    ...     counts = alignment_to_count_df(f.name)
    ...     counts_nogaps = alignment_to_count_df(f.name, ignore_gaps=True)
    ...     freqs_nogaps = alignment_to_count_df(f.name,
    ...                                          ignore_gaps=True,
    ...                                          as_freqs=True)
    ...     freqs_with_ent = alignment_to_count_df(f.name,
    ...                                            ignore_gaps=True,
    ...                                            as_freqs=True,
    ...                                            add_entropy_neff=True,
    ...                                            )

    >>> counts
       isite  -  A  C  G  I  K  L  M  N  S
    0      1  1  0  0  3  0  0  0  0  0  0
    1      2  1  0  0  0  0  3  0  0  0  0
    2      3  3  0  0  0  0  0  0  0  1  0
    3      4  0  4  0  0  0  0  0  0  0  0
    4      5  0  0  3  0  0  0  0  1  0  0
    5      6  3  0  0  0  0  0  0  0  0  1
    6      7  0  0  0  0  1  0  3  0  0  0
    >>> counts_nogaps
       isite  A  C  G  I  K  L  M  N  S
    0      1  0  0  3  0  0  0  0  0  0
    1      2  0  0  0  0  3  0  0  0  0
    2      3  0  0  0  0  0  0  0  1  0
    3      4  4  0  0  0  0  0  0  0  0
    4      5  0  3  0  0  0  0  1  0  0
    5      6  0  0  0  0  0  0  0  0  1
    6      7  0  0  0  1  0  3  0  0  0
    >>> freqs_nogaps.round(2)
       isite    A     C    G     I    K     L     M    N    S
    0      1  0.0  0.00  1.0  0.00  0.0  0.00  0.00  0.0  0.0
    1      2  0.0  0.00  0.0  0.00  1.0  0.00  0.00  0.0  0.0
    2      3  0.0  0.00  0.0  0.00  0.0  0.00  0.00  1.0  0.0
    3      4  1.0  0.00  0.0  0.00  0.0  0.00  0.00  0.0  0.0
    4      5  0.0  0.75  0.0  0.00  0.0  0.00  0.25  0.0  0.0
    5      6  0.0  0.00  0.0  0.00  0.0  0.00  0.00  0.0  1.0
    6      7  0.0  0.00  0.0  0.25  0.0  0.75  0.00  0.0  0.0
    >>> (freqs_with_ent.drop(columns=['entropy', 'n_effective']) ==
    ...  freqs_nogaps).all(None)
    True
    >>> freqs_with_ent[['entropy', 'n_effective']].round(3)
       entropy  n_effective
    0    0.000        1.000
    1    0.000        1.000
    2    0.000        1.000
    3    0.000        1.000
    4    0.811        1.755
    5    0.000        1.000
    6    0.811        1.755

    """
    seqs = [str(s.seq) for s in Bio.AlignIO.read(fasta_alignment, 'fasta')]
    counts = {}
    for i, i_chars in enumerate(zip(*seqs)):
        counts[i + 1] = collections.Counter(i_chars)

    # nested dict to data frame: https://stackoverflow.com/a/54300940
    df = (pd.concat({site: pd.DataFrame(site_counts, index=[0])
                     for site, site_counts in counts.items()},
                    axis=0,
                    ignore_index=True,
                    sort=True)
          .fillna(0)
          .astype(int)
          .assign(isite=lambda x: x.index + 1)
          )
    char_cols = [char for char in df.columns[: -1]
                 if char != '-' or not ignore_gaps]
    df = df[['isite'] + char_cols]

    if as_freqs:
        df = (df
              .set_index('isite')
              .div(df.set_index('isite').sum(axis=1),
                   axis=0)
              .reset_index()
              )
        row_sums = df[char_cols].sum(axis=1)
        if not numpy.allclose(1, row_sums):
            raise RuntimeError('rows do not sum to 1:\n' +
                               str(sorted(row_sums.tolist())))

    if add_entropy_neff:
        if not as_freqs:
            raise ValueError('cannot set `add_entropy_neff` without `as_freq`')
        df = df.assign(
                entropy=lambda x: (x[char_cols]
                                   .apply(lambda r: -(r[r > 0] *
                                                      numpy.log(r[r > 0]) /
                                                      numpy.log(2)
                                                      ).sum(),
                                          axis=1)
                                   ),
                n_effective=lambda x: 2**x['entropy']
                )
        if df['entropy'].any() < 0:
            raise RuntimeError('negative entropy')
        else:
            # get rid of annoying -0 by changing to 0
            df['entropy'] = numpy.abs(df['entropy'])

    return df


def align_prots_mafft(prots, *, mafft='mafft'):
    """Align protein sequences with ``mafft``.

    Parameters
    ----------
    prots : list
        List of sequences as `Biopython.SeqRecord.SeqRecord` objects.
    mafft : str
        Command that resolves to ``mafft`` executable, potentially
        with arbitrary additional options, such as
        ``mafft --reorder``.

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

    mafft_cmds = mafft.split()
    mafft = mafft_cmds[0]
    mafft_args = [arg for arg in mafft_cmds[1:] if arg != '--amino']

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
        res = subprocess.run([mafft, *mafft_args, '--amino', prots_in.name],
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


def aligned_site_map(aln, heads):
    """Get mapping between sequential site numbers of aligned sequences.

    Parameters
    ----------
    aln : Bio.Align.MultipleSeqAlignment
        Alignment of sequences
    heads : list
        Headers of sequences of interest in `aln`, must be at least two.

    Returns
    -------
    pandas.DataFrame
        Columns are sequence names in `heads`. There is a row for alignment
        positions that are non-gapped in at least one sequence in `heads`, and
        the entry in the row gives the number of that site in 1, 2, ...
        numbering of the sequence, or `NaN` if there is no residue for
        sequence at that position in alignment.

    Example
    -------
    Create and parse an alignment:

    >>> with tempfile.TemporaryFile('w+') as f:
    ...     _ = f.write(textwrap.dedent(
    ...             '''
    ...             >seq1
    ...             GK-AC-L
    ...             >seq2
    ...             -K-AM-L
    ...             >seq3
    ...             G-NAC-L
    ...             >seq4
    ...             GK-ACSI
    ...             '''
    ...             ))
    ...     f.flush()
    ...     _ = f.seek(0)
    ...     aln = Bio.AlignIO.read(f, 'fasta')

    Show the alignment site map. The gymnastics with converting to str and
    replacing 'nan' are to make output format the same in both `pandas`
    >= and < 1.0 and earlier versions so doctest passes in both:

    >>> (aligned_site_map(aln, ['seq1', 'seq2', 'seq3'])
    ...  .astype(str).replace('nan', '<NA>'))
       seq1  seq2  seq3
    0     1  <NA>     1
    1     2     1  <NA>
    2  <NA>  <NA>     2
    3     3     2     3
    4     4     3     4
    5     5     4     5

    """
    heads_set = set(heads)
    if len(heads) != len(heads_set):
        raise ValueError(f"non-unique entries in `heads`:\n{heads}")
    if len(heads) < 2:
        raise ValueError('`heads` does not have at least 2 entries')

    aln_heads = [s.description for s in aln]
    aln_heads_set = set(aln_heads)
    if len(aln_heads) != len(aln_heads_set):
        raise ValueError('alignment has non-unique headers')
    missing_heads = heads_set - aln_heads_set
    if missing_heads:
        raise ValueError('following headers in `heads` but not alignment:\n' +
                         '\n'.join(missing_heads))
    select_aln = []
    for s in aln:
        if s.description in heads_set:
            select_aln.append((s.description, str(s.seq)))
    assert len(select_aln) == len(heads)

    aln_length = len(select_aln[0][1])
    if any(aln_length != len(s) for _, s in select_aln):
        raise ValueError('aligned sequences not all of same length')

    site_lists = {head: [] for head in heads}
    indices = {head: 0 for head in heads}

    for i in range(aln_length):
        for head, s in select_aln:
            if s[i] == '-':
                site_lists[head].append(float('nan'))
            else:
                indices[head] += 1
                site_lists[head].append(indices[head])

    df = (pd.DataFrame(site_lists, dtype='Int64')
          [heads]
          .dropna(axis=0, how='all')
          .reset_index(drop=True)
          )

    return df


def pdb_seq_to_number(pdbfile, chain_ids, chain_identity='require_same'):
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


def strip_gaps_to_ref(aln, ref_head):
    """Strip gaps in alignment relative to reference.

    Parameters
    ----------
    aln : Bio.Align.MultipleSeqAlignment
        The alignment of sequences.
    ref_head : str
        Header of reference sequence in `aln`.

    Returns
    -------
    list
        Alignment as list of  2-tuples of strings (`header, sequence`) with
        all gaps relative to reference removed.

    Example
    -------
    >>> with tempfile.TemporaryFile('w+') as f:
    ...     _ = f.write(textwrap.dedent(
    ...             '''
    ...             >seq1
    ...             GK-AC-L
    ...             >seq2
    ...             -K-AM-L
    ...             >seq3
    ...             G-NAC-L
    ...             >seq4
    ...             GK-ACSI
    ...             '''
    ...             ))
    ...     f.flush()
    ...     _ = f.seek(0)
    ...     aln = Bio.AlignIO.read(f, 'fasta')
    >>> strip_gaps_to_ref(aln, 'seq1')  # doctest: +NORMALIZE_WHITESPACE
    [('seq1', 'GKACL'),
     ('seq2', '-KAML'),
     ('seq3', 'G-ACL'),
     ('seq4', 'GKACI')]

    """
    aln = [(s.description, str(s.seq)) for s in aln]
    refseq = [s for head, s in aln if head == ref_head]
    if not refseq:
        raise ValueError(f"alignment lacks `ref_head`:\n{ref_head}")
    elif len(refseq) > 1:
        raise ValueError(f"alignment has {len(refseq)} sequences with "
                         f"`ref_head`:\n{ref_head}")
    else:
        refseq = refseq[0]
    ungapped_indices = [i for i, aa in enumerate(refseq) if aa != '-']
    if not ungapped_indices:
        raise ValueError('no ungapped sites in reference')

    if any(len(refseq) != len(s) for _, s in aln):
        raise ValueError('sequences in `aln` not all same length')

    ungapped_aln = []
    for head, s in aln:
        ungapped_aln.append((head, ''.join(s[i] for i in ungapped_indices)))

    return ungapped_aln


if __name__ == '__main__':
    import doctest
    doctest.testmod()
