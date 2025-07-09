# mrna_petrinet
Project using the PetriNets framework for simulation. This project is originally from https://github.com/lhcnetop/petri_net and is currently being maintained there. But to keep things organized, a refactor will separate the mRNA simulations and code and bring it to this repo, leaving only the main parts of petri net simulation and analysis there.

## Using a local copy of petri_net_core for development

By default, this project depends on the published `petri-net-core` package from PyPI.

If you want to use a local copy of the `petri_net` project (for development or testing changes), run the following command from inside the `mrna_petrinet` directory:

```sh
pip uninstall petri-net-core
pip install -e ../../petri_net
```

This will install your local version in editable mode, so any changes you make to `petri_net` will be reflected immediately in your environment.

To revert to the PyPI version, simply uninstall the local version and reinstall from PyPI:

```sh
pip uninstall petri-net-core
pip install petri-net-core
```
