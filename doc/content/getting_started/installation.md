# Install and Run OpenPronghorn

This guide walks through downloading OpenPronghorn, building it, and running a
small example. You do not need to install MOOSE separately or know how MOOSE is
developed. The commands below take care of the required version for you.

The instructions work on Linux, macOS, and Windows Subsystem for Linux. You
will need an internet connection, several gigabytes of free disk space, and a
terminal. Make sure [Git](https://git-scm.com/downloads) is installed before
starting.

## 1. Install Conda

OpenPronghorn needs a compiler and several supporting libraries. The simplest
way to get them is through Conda, which keeps them together in an isolated
software environment.

First, check whether Conda is already installed:

```bash
conda --version
```

If that command prints a version number, continue to the next section. If it
reports that `conda` was not found, follow the [MOOSE Miniforge installation
instructions](https://mooseframework.inl.gov/getting_started/installation/install_miniconda.html),
then open a new terminal.

If you already had Conda installed, add the package channels used by MOOSE:

```bash
conda config --add channels conda-forge
conda config --add channels https://conda.software.inl.gov/public
```

!alert warning title=Do not use sudo
Conda should be installed in your home directory. Do not use `sudo` with any
of the Conda commands in this guide.

## 2. Prepare the OpenPronghorn environment

Create an environment named `moose`. This downloads the compiler and libraries
needed by OpenPronghorn and may take a few minutes:

!versioner! code
conda create --name moose moose-dev=__VERSIONER_CONDA_VERSION_MOOSE_DEV__
!versioner-end!

When Conda asks whether to proceed, enter `y`. Then activate the environment:

```bash
conda activate moose
```

Your terminal prompt will usually begin with `(moose)` after activation. Run
`conda activate moose` again whenever you open a new terminal to use
OpenPronghorn.

## 3. Download OpenPronghorn

These commands create a `projects` folder in your home directory, download
OpenPronghorn, and download the matching MOOSE source code:

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/idaholab/open_pronghorn.git
cd open_pronghorn
git submodule update --init
```

The last command is important. If it is interrupted, run it again from the
`open_pronghorn` directory.

## 4. Build OpenPronghorn

Turn the downloaded source code into the OpenPronghorn program:

```bash
cd ~/projects/open_pronghorn
make -j 6
```

The build can take a while. The `6` tells the computer to perform up to six
compile tasks at once. On a computer with limited memory, use `make -j 2`
instead. When the build finishes, the program `open_pronghorn-opt` will be in
the `open_pronghorn` directory.

## 5. Check the installation

Run the automated tests to make sure the installation works:

```bash
./run_tests -j 6
```

A successful installation ends with no failed tests. It is normal for some
tests to be skipped when optional software is unavailable.

## 6. Run a first example

From the `open_pronghorn` directory, run the included diffusion example:

```bash
./open_pronghorn-opt -i test/tests/kernels/simple_diffusion/simple_diffusion.i
```

Here, `-i` tells OpenPronghorn which input file to use. The run creates
`simple_diffusion_out.e`, an Exodus results file that can be opened in a
compatible visualization program such as ParaView.

Other OpenPronghorn simulations use the same pattern:

```bash
./open_pronghorn-opt -i path/to/your_input_file.i
```

For everyday use, the only setup normally needed in a new terminal is:

```bash
conda activate moose
cd ~/projects/open_pronghorn
```

## Updating later

To get the latest OpenPronghorn changes, update the source and its matching
MOOSE version:

```bash
conda activate moose
cd ~/projects/open_pronghorn
git pull
git submodule update --init
```

If the update requires a different software environment, update it with:

!versioner! code
conda install moose-dev=__VERSIONER_CONDA_VERSION_MOOSE_DEV__
!versioner-end!

Then rebuild and check the installation:

```bash
make -j 6
./run_tests -j 6
```

## If something goes wrong

- +`conda` is not found:+ Close and reopen the terminal after installing
  Miniforge. If it is still unavailable, return to the linked Miniforge guide
  and complete the `conda init` step.
- +`moose/framework/Makefile` is missing:+ From the `open_pronghorn` directory,
  run `git submodule update --init` again.
- +The build stops with a package version message:+ Run the `conda install`
  command shown in that message, then restart the build.
- +The build is killed or the computer becomes unresponsive:+ Use fewer tasks,
  for example `make -j 2`.
- +Tests fail:+ Search the [OpenPronghorn issue
  tracker](https://github.com/idaholab/open_pronghorn/issues) or ask for help in
  [OpenPronghorn Discussions](https://github.com/idaholab/open_pronghorn/discussions).
  Include the failed test names and error messages when asking for help.

Advanced installations, including high-performance computing systems and
manual compiler setups, are covered by the [MOOSE installation
guide](https://mooseframework.inl.gov/getting_started/installation/index.html).
