``pdb_prot_align`` documentation
===================================================

``pdb_prot_align`` is a program written by `the Bloom lab <https://research.fhcrc.org/bloom/en.html>`_ to align proteins to `Protein Data Bank (PDB) <https://www.rcsb.org/>`_ structures.

Making such alignments can be useful when you want to relate the evolutionary information in a protein alignment to the structure of a protein.
However, the alignments are often difficult because the sequential numbering of sequences used in bioinformatics can differ from the numbering of the protein structure file.
The goal of ``pdb_prot_align`` is to simplify reconciliation of this numbering. 

You may want to simply run ``pdb_prot_align`` via its :ref:`cli`.
However, if you prefer to work within Python you can use the :mod:`pdb_prot_align` package, including the top-level :func:`pdb_prot_align.main.run` that mimics the :ref:`cli`.

See below for more information and examples:

Contents
----------
.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 1

   installation
   cli
   examples
   pdb_prot_align
   package_index
   acknowledgments
