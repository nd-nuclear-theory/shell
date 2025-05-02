# Installation instructions for shell #

This repository includes submodules, so be sure to do a recursive clone:

  ~~~~~~~~~~~~~~~~
  % git clone --recursive https://github.com/nd-nuclear-theory/shell.git
  ~~~~~~~~~~~~~~~~

You will also want to clone the `ndconfig` repository:

  ~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/ndconfig.git
  ~~~~~~~~~~~~~~~~

Then please refer to the full installation instructions, which are found in
`ndconfig/INSTALL.md` ("Installation instructions for ND nuclear theory projects
using `ndconfig`").  The `ndconfig` repository is also where you will find
several example `config.mk` files for use in the installation.

*Note:* If you are planning on using `shell` together with `mfdn`, and plan to
be running under the Notre Dame scripting (in `mcscript-ncci`), then you will
want to start by looking at the quickstart guide in the `mcscript-ncci`
repository, under the `docs` directory.  You can find it online here:

  ~~~~~~~~~~~~~~~~
  https://github.com/nd-nuclear-theory/mcscript-ncci/blob/master/docs/nd-mfdn-quickstart-guide.md
  ~~~~~~~~~~~~~~~~
