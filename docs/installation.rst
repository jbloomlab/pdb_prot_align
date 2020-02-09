Installation
--------------

Installing ``pdb_prot_align`` itself
+++++++++++++++++++++++++++++++++++++

``pdb_prot_align`` requires Python 3.6 or higher.

The easiest way to install ``pdb_prot_align`` is from `PyPI <https://pypi.org/>`_ using `pip <https://pip.pypa.io>`_ with::

    pip install pdb_prot_align

You also need to ensure that `mafft <https://mafft.cbrc.jp/alignment/software/>`_ is installed so it is available at the command line as ``mafft``.

The source code for ``pdb_prot_align`` is available on GitHub at https://github.com/jbloomlab/pdb_prot_align.

Installing PyMol and nglview
++++++++++++++++++++++++++++++
In some cases, you may want to visualize the results on protein structures using pymol_ or nglview_.
This can be done as described in the :doc:`pymol_example` and the :doc:`nglview_example` examples by taking advantage of the :mod:`pdb_prot_align.pymol_struct` and :mod:`pdb_prot_align.nglview_struct` modules.
In order to uses these, you must install pymol_ or nglview_, which is substantially more complex than the simple instructions above to install ``pdb_prot_align``.

pymol_ can be installed a number of ways, as can nglview_.
We suggest doing it via conda_.
The environment.yml_ file in the `pdb_prot_align GitHub repo`_ can be used to install pymol_ and nglview_.
To do this, just follow the `instructions to create a conda environment from a YAML file <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_.
However assuming that you're not building from the source in the `pdb_prot_align GitHub repo`_, you should change the lines in the environment.yml_ file that read::

    - pip:
      - -e .

to instead read::

    - pip:
      - pdb_prot_align

For nglview_, after building and activating the conda_ environment, you then need to run the following command (see `here <https://github.com/jbloomlab/pdb_prot_align>`_)::

    jupyter-nbextension enable nglview --py --sys-prefix

If you are using `mybinder <https://mybinder.org/>`_, this can be done via a `postBuild <https://mybinder.readthedocs.io/en/latest/config_files.html#postbuild-run-code-after-installing-the-environment>`_ file with the `following contents <https://github.com/jbloomlab/pdb_prot_align/blob/master/binder/postBuild>`_ placed in ``./binder/postBuild``.

Note that if you do not want your pymol_ images to have a watermark, you'll also have to manually install a pymol_ license.
Note also note that the current version of nglview_ may only work in a Jupyter **notebook**, if you want to run through Jupyter **lab** you may have to `follow the extra steps here <https://github.com/arose/nglview/issues/801#issuecomment-492744453>`_.

.. _pymol: https://pymol.org/
.. _nglview: https://github.com/arose/nglview
.. _conda: https://docs.conda.io
.. _environment.yml: https://github.com/jbloomlab/pdb_prot_align/blob/master/environment.yml
.. _`pdb_prot_align GitHub repo`: https://github.com/jbloomlab/pdb_prot_align
