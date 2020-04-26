"""
=======
main
=======

Main PDB / protein alignment function.

You can run all of the command-line functionality of ``pdb_prot_align`` via
the :func:`run` function.

"""


import argparse
import os
import re
import sys
import tempfile
import textwrap

import Bio.AlignIO
import Bio.Seq
import Bio.SeqIO

import pandas as pd

import pdb_prot_align
import pdb_prot_align.utils


def _run_main():
    parser = get_parser()
    kwargs = vars(parser.parse_args(sys.argv[1:]))
    run(**kwargs)


def _str2bool(v):
    """Convert string to boolean value."""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# arguments for command-line parser, specified this way so can also
# be used in docs for `run`
_ARGS = {
         'protsfile': {
              'required': True,
              'help': 'input FASTA file of protein sequences',
              'type': str,
              },
         'refprot_regex': {
              'required': True,
              'help': 'regex for reference protein header in `protsfile`',
              'type': str,
              },
         'pdbfile': {
              'required': True,
              'help': 'input PDB file',
              'type': str,
              },
         'chain_ids': {
              'required': True,
              'help': 'chains in PDB file to align; all chains aligning to a '
                      'site must share the same residue number and amino-acid '
                      'or an error will be raised',
              'type': str,
              'nargs': '+',
              },
         'outprefix': {
              'required': True,
              'help': 'prefix for output files (can be / include directory): '
                      '"alignment.fa" (alignment with gaps relative to '
                      'reference stripped); "alignment_unstripped.fa" ('
                      'non-stripped alignment with PDB chains still included);'
                      ' "sites.csv" (sequential sites in reference, PDB sites'
                      ', PDB chains, wildtype in reference, wildtype in PDB, '
                      'site entropy in bits, n effective amino acids at site, '
                      'amino acid, frequency of amino acid)',
              'type': str,
              },
         'ignore_gaps': {
              'required': False,
              'default': True,
              'help': 'ignore gaps (-) when calculating frequencies, number '
                      'effective amino acids, entropy',
              'type': _str2bool,
              },
         'drop_pdb': {
              'required': False,
              'default': True,
              'help': 'drop PDB protein chains from "alignment.fa" and '
                      'computation of stats in "sites.csv" output files',
              'type': _str2bool,
              },
         'mafft': {
             'required': False,
             'default': 'mafft',
             'help': 'path to `mafft`, potentially with additional args '
                     'such as "mafft --reorder" (if multiple args, it all '
                     'needs to be in quotes)',
             'type': str,
             },
         }


def get_parser():
    """Command-line `argparse.ArgumentParser` for ``pdb_prot_align``."""
    parser = argparse.ArgumentParser(
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                    description='Align proteins to reference and PDB.',
                    )

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {0}'.format(
                                pdb_prot_align.__version__))

    for param, param_d in _ARGS.items():
        parser.add_argument(f"--{param}", **param_d)

    return parser


def run(protsfile,
        refprot_regex,
        pdbfile,
        chain_ids,
        outprefix,
        ignore_gaps=True,
        drop_pdb=True,
        mafft='mafft',
        ):
    """Run main function to align proteins to reference and PDB chain(s).

    Note
    ----
    This function implements the full command-line functionality of
    ``pdb_prot_align`` via a Python function.

    Parameters
    ----------
    {_PARAM_DOCS}

    """
    # get output file names and make output directory if needed
    os.makedirs(os.path.dirname(outprefix), exist_ok=True)
    if not os.path.isdir(outprefix):
        outprefix = f"{outprefix}_"
    alignment = f"{outprefix}alignment.fa"
    alignment_unstripped = f"{outprefix}unstripped_alignment.fa"
    csv = f"{outprefix}sites.csv"

    print(f"\nRunning `pdb_prot_align` {pdb_prot_align.__version__}\n")
    if not os.path.isfile(pdbfile):
        raise IOError(f"no `pdbfile` of {pdbfile}")

    # parse PDB chains into dataframe mapping sequential to PDB residues
    print(f"Parsing PDB {pdbfile} chains {' '.join(chain_ids)}")
    pdb_dfs = {}   # dataframe mapping sequential to PDB residue for each chain
    pdb_prots = {}  # SeqRecord for each chain
    pdb_name = os.path.splitext(os.path.basename(pdbfile))[0]
    if len(chain_ids) != len(set(chain_ids)):
        raise ValueError(f"`chain_ids` not unique: {chain_ids}")
    for chain_id in chain_ids:
        pdb_dfs[chain_id] = pdb_prot_align.utils.pdb_seq_to_number(pdbfile,
                                                                   chain_id)
        with tempfile.NamedTemporaryFile('w+') as f:
            f.write(f">PDB_{pdb_name}_{chain_id}\n" +
                    ''.join(pdb_dfs[chain_id]['wildtype']))
            f.flush()
            f.seek(0)
            pdb_prots[chain_id] = Bio.SeqIO.read(f, 'fasta')
        print(f"For chain {chain_id}, parsed {len(pdb_dfs[chain_id])} residues"
              f", ranging from {pdb_dfs[chain_id]['pdb_site'].iloc[0]} to "
              f"{pdb_dfs[chain_id]['pdb_site'].iloc[-1]} in PDB numbering.")
    assert all(len(pdb_dfs[chain_id]) == len(pdb_prots[chain_id])
               for chain_id in chain_ids)

    # read in all the proteins to align
    if not os.path.isfile(protsfile):
        raise IOError(f"no `protsfile` file {protsfile}")
    prots = list(Bio.SeqIO.parse(protsfile, 'fasta'))
    print(f"\nRead {len(prots)} sequences from {protsfile}")

    # ensure no proteins have same header we used for PDB
    pdb_prot_headers = [pdb_prot.description.strip() for pdb_prot
                        in pdb_prots.values()]
    if any(p.description.strip() in pdb_prot_headers for p in prots):
        raise ValueError(f"sequences in {protsfile} cannot have header in "
                         f"{pdb_prot_headers}")

    # get reference protein from set of proteins
    refprot = [s for s in prots if re.search(refprot_regex, s.description)]
    if not refprot:
        raise ValueError(f"no headers in {protsfile} match `refprot_regex` "
                         f"{refprot_regex}")
    elif len(refprot) > 1:
        raise ValueError(f"multiple headers in {protsfile} match "
                         f"`refprot_regex` {refprot_regex}:\n" +
                         '\n'.join(s.description for s in refprot))
    else:
        refprot = refprot[0]
        print(f"Reference protein is of length {len(refprot)} and has "
              f"the following header:\n{refprot.description}\n")

    # make alignment of all proteins and the PDB sequence
    print(f"Using `mafft` to align sequences to {alignment_unstripped}")
    aln = pdb_prot_align.utils.align_prots_mafft(
                                    prots=prots + list(pdb_prots.values()),
                                    mafft=mafft)
    assert len(aln) == len(prots) + len(pdb_prots)
    Bio.AlignIO.write(aln, alignment_unstripped, 'fasta')

    # get mapping of reference protein and PDB sequential sites
    ref_pdb_site_map = pdb_prot_align.utils.aligned_site_map(
                            aln,
                            [refprot.description, *pdb_prot_headers]
                            )
    assert all(len(p) == len(ref_pdb_site_map[p.description].dropna())
               for p in [refprot, *pdb_prots.values()])

    # remove PDB sequences?
    if drop_pdb:
        print('Dropping PDB chains from alignment')
        aln = [s for s in aln if s.description not in pdb_prot_headers]
        assert len(aln) == len(prots)

    # strip gaps relative to reference sequence
    print(f"Stripping gaps relative to reference {refprot.description}")
    aln = pdb_prot_align.utils.strip_gaps_to_ref(aln,
                                                 refprot.description)
    assert all(len(refprot) == len(s) for _, s in aln)
    print(f"Writing gap-stripped alignment to {alignment}\n")
    with open(alignment, 'w') as f:
        f.write('\n'.join(f">{head}\n{s}" for head, s in aln))

    # remove sites gapped in reference from `ref_pdb_site_map`
    ref_pdb_site_map = (ref_pdb_site_map
                        .dropna(axis=0, subset=[refprot.description])
                        )

    # map (pdb_chain, pdb_isite) to pdb_site or pdb_wildtype
    to_pdb = {}
    for valstr in ['pdb_site', 'pdb_wildtype']:
        to_pdb[valstr] = {('NaN', 'NaN'): 'NaN'}
        for chain_id, pdb_df in pdb_dfs.items():
            for isite, val in (pdb_df
                               .rename(columns={'wildtype': 'pdb_wildtype'})
                               .set_index('sequential')
                               [valstr]
                               .to_dict()
                               .items()
                               ):
                to_pdb[valstr][chain_id, isite] = val

    # now create data frame with PDB sites / wildtypes
    df = (
        ref_pdb_site_map
        .rename(columns={old: new for old, new in
                         [(refprot.description, 'isite'),
                          *[(pdb_prot.description, key) for key, pdb_prot
                            in pdb_prots.items()]
                          ]
                         }
                )
        .assign(wildtype=lambda x: (x['isite']
                                    .map({i + 1: wt for i, wt in
                                          enumerate(str(refprot.seq))},
                                         na_action='ignore')
                                    )
                )
        .melt(id_vars=['isite', 'wildtype'],
              var_name='pdb_chain',
              value_name='pdb_isite',
              )
        .assign(pdb_chain=lambda x: (x['pdb_chain']  # chain NaN if site is NaN
                                     .where(pd.notna(x['pdb_isite']), 'NaN')
                                     ),
                pdb_isite=lambda x: x['pdb_isite'].fillna('NaN'),
                )
        .drop_duplicates()
        .assign(
            pdb_site=lambda x: x.apply(lambda r: (to_pdb['pdb_site']
                                                  [r.pdb_chain, r.pdb_isite]
                                                  ),
                                       axis=1),
            pdb_wildtype=lambda x: x.apply(lambda r: (to_pdb['pdb_wildtype']
                                                      [r.pdb_chain,
                                                       r.pdb_isite]
                                                      ),
                                           axis=1),
            nchains=lambda x: x.groupby('isite').pdb_chain.transform('count')
            )
        .query('(pdb_chain != "NaN") or (nchains == 1)')
        .reset_index(drop=True)
        .replace({'NaN': float('nan')})
        [['isite', 'wildtype', 'pdb_chain', 'pdb_site', 'pdb_wildtype']]
        )

    # check that each site aligns to unique PDB wildtype
    nonunique_isites = (
        df
        .drop(columns=['pdb_chain'])
        .drop_duplicates()
        .groupby('isite')
        .size()
        .rename('n')
        .reset_index()
        .query('n > 1')
        ['isite']
        .tolist()
        )
    if nonunique_isites:
        raise ValueError('Sites map to non-unique PDB site / wildtype:\n' +
                         str(df.query('isite in @nonunique_isites')
                             .sort_values(['isite', 'pdb_chain'])))

    # add amino-acid frequencies and entropies / n_effective
    df = (
        df
        .merge(pdb_prot_align.utils.alignment_to_count_df(
                                            alignment,
                                            ignore_gaps=ignore_gaps,
                                            as_freqs=True,
                                            add_entropy_neff=True),
               on='isite',
               validate='many_to_one',
               )
        .melt(id_vars=['isite', 'wildtype', 'pdb_chain', 'pdb_site',
                       'pdb_wildtype', 'entropy', 'n_effective'],
              var_name='amino_acid',
              value_name='frequency',
              )
        .sort_values(['isite', 'pdb_chain', 'amino_acid'])
        )

    print(f"Writing CSV with detailed information to {csv}")
    df.to_csv(csv, index=False, float_format='%.5f')

    print(f"\nProgram complete.\n")


# complete docs for `run` with parameter specs parsed from `_ARGS`
_PARAM_DOCS = []
for param, param_d in _ARGS.items():
    if 'choices' in param_d:
        param_type = f"{{{', '.join(param_d['choices'])}}}"
    elif 'nargs' in param_d and param_d['nargs'] == '+':
        param_type = 'list'
    elif param_d['type'] == _str2bool:
        param_type = 'bool'
    else:
        param_type = param_d['type'].__name__
    param_help = textwrap.indent('\n'.join(textwrap.wrap(param_d['help'], 70)),
                                 '        ')
    _PARAM_DOCS.append(f"    {param} : {param_type}\n{param_help}")
_PARAM_DOCS = '\n'.join(_PARAM_DOCS)
del param, param_d, param_type, param_help

# add parameter specs to docs for `run`
run.__doc__ = run.__doc__.format(_PARAM_DOCS=_PARAM_DOCS.lstrip())


if __name__ == '__main__':
    import doctest
    doctest.testmod()
