"""Main PDB / protein alignment function."""


import argparse
import os
import re
import sys
import tempfile
import textwrap

import Bio.AlignIO
import Bio.Seq
import Bio.SeqIO

import numpy

import pdb_prot_align
import pdb_prot_align.utils


def _run_main():
    parser = get_parser()
    kwargs = vars(parser.parse_args(sys.argv[1:]))
    run(**kwargs)


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
              'help': 'chains in PDB file, must be identical monomers if > 1',
              'type': str,
              'nargs': '+',
              },
         'alignment': {
              'required': True,
              'help': 'created FASTA alignment of proteins with gaps relative '
                      'to reference are stripped; non-stripped alignment with '
                      'PDB protein written to file with suffix "_unstripped"',
              'type': str,
              },
         'csv': {
              'required': True,
              'help': 'created CSV with sequential sites in reference, PDB '
                      'sites, wildtype in reference, wildtype in PDB, amino '
                      'amino acid, frequency of amino acid, site entropy in '
                      'bits, number of effective amino acids at site',
              'type': str,
              },
         'chain_identity': {
              'required': False,
              'default': 'union',
              'choices': ['union', 'intersection', 'require_same'],
              'help': 'if sites differ among chains, get union, intersection, '
                      'or raise error if not same sites',
              },
         'ignore_gaps': {
              'required': False,
              'default': True,
              'help': 'ignore gaps (-) when calculating frequencies, number '
                      'effective amino acids, entropy',
              'type': bool,
              },
         'drop_pdb': {
              'required': False,
              'default': True,
              'help': 'drop protein in `pdbfile` from `alignment` and '
                      '`csv` stats',
              'type': bool,
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
        alignment,
        csv,
        chain_identity='union',
        ignore_gaps=True,
        drop_pdb=True,
        ):
    """Run main function to align proteins to reference and PDB chain(s).

    Parameters
    ----------
    {_PARAM_DOCS}

    """
    print(f"\nRunning `pdb_prot_align` {pdb_prot_align.__version__}\n")
    if not os.path.isfile(pdbfile):
        raise IOError(f"no `pdbfile` of {pdbfile}")

    # parse PDB sequence into dataframe mapping sequential to PDB residues
    print(f"Parsing PDB {pdbfile} chains {' '.join(chain_ids)}")
    pdb_df = pdb_prot_align.utils.pdb_seq_to_number(pdbfile,
                                                    chain_ids,
                                                    chain_identity)
    print(f"Parsed sequence of {len(pdb_df)} residues, ranging "
          f"from {pdb_df['pdb_site'].values[0]} to "
          f"{pdb_df['pdb_site'].values[-1]} in PDB numbering.\n")

    # make a BioPython SeqRecord with PDB protein
    with tempfile.NamedTemporaryFile('w+') as f:
        f.write(f">PDB_{os.path.splitext(os.path.basename(pdbfile))[0]}\n" +
                ''.join(pdb_df['wildtype'].values))
        f.flush()
        f.seek(0)
        pdb_prot = Bio.SeqIO.read(f, 'fasta')
    assert len(pdb_prot) == len(pdb_df)

    # read in all the proteins to align
    if not os.path.isfile(protsfile):
        raise IOError(f"no `protsfile` file {protsfile}")
    prots = list(Bio.SeqIO.parse(protsfile, 'fasta'))
    print(f"Read {len(prots)} sequences from {protsfile}")

    # ensure no proteins have same header we used for PDB
    if any(p.description.strip() == pdb_prot.description.strip()
           for p in prots):
        raise ValueError(f"sequences in {protsfile} cannot have header "
                         f"{pdb_prot.description.strip()}")

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
    base, ext = os.path.splitext(alignment)
    alignment_unstripped = f"{base}_unstripped{ext}"
    print(f"Using `mafft` to align sequences to {alignment_unstripped}")
    aln = pdb_prot_align.utils.align_prots_mafft(prots + [pdb_prot])
    assert len(aln) == len(prots) + 1
    Bio.AlignIO.write(aln, alignment_unstripped, 'fasta')

    # get mapping of reference protein and PDB sequential sites
    ref_pdb_site_map = pdb_prot_align.utils.aligned_site_map(
                            aln,
                            [refprot.description, pdb_prot.description],
                            )
    assert all(len(p) == len(ref_pdb_site_map[p.description].dropna())
               for p in [pdb_prot, refprot])

    # remove PDB sequence?
    if drop_pdb:
        print(f"Dropping PDB sequence {pdb_prot.description} from alignment")
        aln = [s for s in aln if s.description != pdb_prot.description]
        assert len(aln) == len(prots)

    # strip gaps relative to reference sequence
    print(f"Stripping gaps relative to reference {refprot.description}")
    aln = pdb_prot_align.utils.strip_gaps_to_ref(aln,
                                                 refprot.description)
    assert all(len(refprot) == len(s) for _, s in aln)
    print(f"Writing gap-stripped alignment to {alignment}\n")
    with open(alignment, 'w') as f:
        f.write('\n'.join(f">{head}\n{s}" for head, s in aln))

    # Create data frame that will become `csv` output.
    # First get site numbers and wildtype.
    df = (
        ref_pdb_site_map
        .rename(columns={refprot.description: 'isite',
                         pdb_prot.description: 'pdb_isite'})
        .assign(pdb_site=lambda x: (x['pdb_isite']
                                    .map(pdb_df
                                         .set_index('sequential')
                                         ['pdb_site']
                                         .to_dict(),
                                         na_action='ignore')
                                    ),
                wildtype=lambda x: (x['isite']
                                    .map({i + 1: wt for i, wt in
                                          enumerate(str(refprot.seq))},
                                         na_action='ignore')
                                    ),
                pdb_wildtype=lambda x: (x['pdb_isite']
                                        .map(pdb_df
                                             .set_index('sequential')
                                             ['wildtype']
                                             .to_dict(),
                                             na_action='ignore')
                                        ),
                )
        [['isite', 'pdb_site', 'wildtype', 'pdb_wildtype']]
        .merge(pdb_prot_align.utils.alignment_to_count_df(
                                            alignment,
                                            ignore_gaps=ignore_gaps,
                                            as_freqs=True),
               on='isite',
               validate='one_to_one',
               )
        .set_index(['isite', 'pdb_site', 'wildtype', 'pdb_wildtype'])
        .assign(
            entropy=lambda x: x.apply(lambda r: -(r[r > 0] *
                                                  numpy.log(r[r > 0])).sum(),
                                      axis=1),
            n_effective=lambda x: numpy.exp(x['entropy'])
            )
        .reset_index()
        .melt(id_vars=['isite', 'pdb_site', 'wildtype', 'pdb_wildtype',
                       'entropy', 'n_effective'],
              var_name='amino_acid',
              value_name='frequency',
              )
        .sort_values(['isite', 'amino_acid'])
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
