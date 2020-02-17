=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com>`_.

0.3.0
------

Added
++++++
* ``genbank`` module.

Fixed
+++++
* Fixed bug when PDB has insertions relative to reference sequence.

0.2.1
------

Added
+++++
* ``remove_widget_view`` option to ``nglview_struct.render_html`` to enable embedding of structure widgets in HTML without showing them twice.

0.2.0
------

Added
++++++
* Can handle multiple distinct non-overlapping chains; adds a ``pdb_chain`` column to output CSV.

* Changed output to specify a prefix rather than individual ``--alignment`` and ``--csv`` options.

* Added ``add_entropy_neff`` option to ``alignment_to_count_df``

* Added the ``colorschemes`` module.

* Added ``pymol_struct`` module to color ``pymol`` structures.

* Added ``nglview_struct`` module to color ``nglview`` structures.

Fixed
+++++
* Made ``--drop_pdb`` and ``--ignore_gaps`` options work correctly (before they did nothing due to command-line parsing error)

* Allow ``pandas`` >= 1.0

0.1.1
-----

Added
+++++
* ``main.py`` will create the output directory(ies) if they do not exist

0.1.0
-----
Initial release
