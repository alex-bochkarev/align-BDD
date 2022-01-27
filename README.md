# Accompanying code for the paper: Aligning BDDs
## Contact info
- Alexey Bochkarev | [contact info](https://www.bochkarev.io/contact/)

## General information and software requirements
The repo is built around the main [makefile](https://en.wikipedia.org/wiki/Makefile) which contains all recipes to build the figures from the paper.

- `TODO:` software requirements. Python packages, using =`venv`:
```bash
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install -r requirements.txt
```

- `TODO:` notes on cluster

## Repo structure
- `TODO:` file structure and a quick description

## Steps to reproduce
To rebuild all the figures from scratch:
```bash
$ make clean && make figures
```

Note that it will almost surely take time, so at least one might be interested in utilizing the 
cluster  or at least using the `make`'s parallel capabilities (`-j` flag).

- `TODO:` PBS instructions

- notes on makefile
- runtimes here?
