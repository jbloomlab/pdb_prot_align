"""Main PDB / protein alignment function."""


import argparse
import sys
import textwrap

import pdb_prot_align


def _run_main():
    parser = get_parser()
    kwargs = vars(parser.parse_args(sys.argv[1:]))
    run(**kwargs)


# arguments for command-line parser and `pdb_prot_align`
_ARGS = {
         'prots': {
              'required': True,
              'help': 'Input FASTA file of protein sequences.',
              'type': str,
              },
         'refprot': {
              'required': True,
              'help': 'Unique regex for reference protein in `prots`.',
              'type': str,
              },
         'pdbfile': {
              'required': True,
              'help': 'Input PDB file.',
              'type': str,
              },
         'chain_ids': {
              'required': True,
              'help': 'Chains in PDB file, must be identical monomers if > 1.',
              'type': list,
              'nargs': '+',
              },
         'alignment': {
              'required': True,
              'help': 'Created FASTA alignment of proteins.',
              'type': str,
              },
         'csv': {
              'required': True,
              'help': 'Created tidy CSV with sequential sites ("isite"), '
                      'PDB sites ("pdb_site"), wildtype ("wildtype") in '
                      'refprot, amino acid identity ("amino acid"), '
                      'frequency of amino acid ("frequency"), and '
                      'number effective amino acids at site ("neffective").',
              'type': str,
              },
         'chain_identity': {
              'required': False,
              'default': 'union',
              'choices': ['union', 'intersection', 'exact'],
              'help': 'If sites differ among chains, get union, intersection '
                      'or raise error if not exact same sites.',
              },
         'drop_pdb': {
              'required': False,
              'default': True,
              'help': 'Drop protein in `pdbfile` from `alignment` and '
                      '`csv` stats.',
              'type': bool,
              },
         'allow_stop': {
              'required': False,
              'default': False,
              'help': 'Treat stop ("*") codons as valid amino acid.',
              'type': bool,
              },
         }


def get_parser():
    """Command-line `argparse.ArgumentParser` for ``pdb_prot_align``."""
    parser = argparse.ArgumentParser(
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                    description='Align proteins to reference and PDB chain(s)',
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
    else:
        param_type = param_d['type'].__name__
    param_help = textwrap.indent('\n'.join(textwrap.wrap(param_d['help'], 70)),
                                 '        ')
    _PARAM_DOCS.append(f"    {param} : {param_type}\n{param_help}")
_PARAM_DOCS = '\n'.join(_PARAM_DOCS)
del param, param_d, param_type, param_help


def run(prots,
        refprot,
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
    pass


pdb_prot_align.__doc__ = pdb_prot_align.__doc__.format(
                                _PARAM_DOCS=_PARAM_DOCS.lstrip())


if __name__ == '__main__':
    import doctest
    doctest.testmod()
