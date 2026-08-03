<p align="center">
  <img src="doc/content/media/open_pronghorn_logo.png" alt="OpenPronghorn logo" width="300">
</p>

# OpenPronghorn

[![CIVET status](https://img.shields.io/github/checks-status/idaholab/open_pronghorn/devel?label=CIVET)](https://civet.inl.gov/repo/1471/)

OpenPronghorn is an engineering-scale thermal-hydraulics application built on
the [MOOSE framework](https://mooseframework.inl.gov). It combines MOOSE's
finite-volume flow, heat-transfer, reactor, subchannel, and thermal-hydraulics
capabilities with OpenPronghorn validation cases and application-specific
documentation.

## Getting started

OpenPronghorn is built from source on Linux, macOS, and Windows Subsystem for
Linux. The [installation guide](doc/content/getting_started/installation.md)
walks through the initial setup, downloading and building OpenPronghorn, and
running a first example.

With the development environment active and the repository cloned, the basic
workflow is:

```bash
git submodule update --init
make -j 6
./run_tests -j 6
```

A successful build creates `open_pronghorn-opt` in the repository root.

## Resources

- [Installation guide](doc/content/getting_started/installation.md)
- [OpenPronghorn documentation](https://openpronghorn.inl.gov)
- [MOOSE documentation](https://mooseframework.inl.gov)
- [Questions and discussions](https://github.com/idaholab/open_pronghorn/discussions)
- [Issue tracker](https://github.com/idaholab/open_pronghorn/issues)

OpenPronghorn is licensed under the [GNU Lesser General Public License,
version 2.1](LICENSE). See [NOTICE.txt](NOTICE.txt) for attribution and government
license information.
