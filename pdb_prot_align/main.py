"""Main PDB / protein alignment function."""


import argparse
import os
import re
import sys
import textwrap

import Bio.AlignIO
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

import pdb_prot_align
import pdb_prot_align.utils


def _run_main():
    parser = get_parser()
    kwargs = vars(parser.parse_args(sys.argv[1:]))
    run(**kwargs)


# arguments for command-line parser and `pdb_prot_align`
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
                      'to reference are stripped, the non-stripped alignment '
                      'written to file with suffix "_unstripped"',
              'type': str,
              },
         'csv': {
              'required': True,
              'help': 'created CSV with sequential sites in reference, PDB '
                      'sites, wildtype in reference, wildtype in PDB, amino '
                      'amino acid, frequencye of amino acid, number of '
                      'effective amino acids at site',
              'type': str,
              },
         'chain_identity': {
              'required': False,
              'default': 'union',
              'choices': ['union', 'intersection', 'require_same'],
              'help': 'if sites differ among chains, get union, intersection, '
                      'or raise error if not same sites',
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


# Parameters doc string created from `_ARGS`
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


def run(protsfile,
        refprot_regex,
        pdbfile,
        chain_ids,
        alignment,
        csv,
        chain_identity='union',
        drop_pdb=True,
        allow_stop=False,
        ):
    """Run main function to align proteins to reference and PDB chain(s).

    Parameters
    ----------
    {_PARAM_DOCS}

    """
    print(f"\nRunning `pdb_prot_align` {pdb_prot_align.__version__}\n")
    if not os.path.isfile(pdbfile):
        raise IOError(f"no `pdbfile` of {pdbfile}")

    print(f"Parsing PDB {pdbfile} chains {' '.join(chain_ids)}")
    pdb_df = pdb_prot_align.utils.pdb_seq_to_number(pdbfile,
                                                    chain_ids,
                                                    chain_identity)
    print(f"Parsed sequence of {len(pdb_df)} residues, ranging "
          f"from {pdb_df['pdb_site'].values[0]} to "
          f"{pdb_df['pdb_site'].values[-1]} in PDB numbering.\n")
    pdb_prot = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(''.join(pdb_df['wildtype'].values)),
            id=f"PDB_{os.path.splitext(os.path.basename(pdbfile))[0]}",
            description='',
            )
    assert len(pdb_prot) == len(pdb_df)

    if not os.path.isfile(protsfile):
        raise IOError(f"no `protsfile` file {protsfile}")
    prots = list(Bio.SeqIO.parse(protsfile, 'fasta'))
    print(f"Read {len(prots)} sequences from {protsfile}")
    if any(p.description.strip() == pdb_prot.description.strip()
           for p in prots):
        raise ValueError(f"sequences in {protsfile} cannot have header "
                         f"{pdb_prot.description.strip()}")

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

    base, ext = os.path.splitext(alignment)
    alignment_unstripped = f"{base}_unstripped{ext}"
    print(f"Using `mafft` to align sequences to {alignment_unstripped}")
    aln = pdb_prot_align.utils.align_prots_mafft(prots + [pdb_prot])
    Bio.AlignIO.write(aln, alignment_unstripped, 'fasta')


run.__doc__ = run.__doc__.format(_PARAM_DOCS=_PARAM_DOCS.lstrip())


if __name__ == '__main__':
    import doctest
    doctest.testmod()
