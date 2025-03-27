## Repo for: PAPER-ID

The structure of the repository is somewhat self-explanatory.

> What else do I want to add here?

### Copernicus data

The copernicus data was obtained using the _Copernicus Marine Service Toolbox
Command Line Interface (CLI)_. For details on how to install it, we refer to
this [webpage](https://pypi.org/project/copernicusmarine/).


The file `north-atlantic.json` specifies the product we are downloading and the
study region we are interested in. The latter is defined by a bounding box. To
download the data, one should run (after installing the CLI):
```
copernicusmarine subset --request-file north-atlantic.json
```
