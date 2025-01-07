# Contributor Guide

## Environment setup and workflow
- Clone the repo
- Create a fresh python environemt with `.binder/requirements.txt` file
- Open `jupyter lab` or `jupyter notebook` (so that jupytext extension can take effect)
- Then right-click on a MyST markdown notebook, choose `Open With` -> `Notebook`
or `Jupytext Notebook` (more info [here](https://nasa-navo.github.io/navo-workshop/00_SETUP.html#handling-notebooks-in-myst-markdown-format)).

> Note: If you'd like to work with this repository in an editor/IDE of your choice,
make sure that you have relevant extensions to parse MyST markdown notebooks to
interactive python notebook (ipynb) and to reflect the changes made there to markdown
format:
>- VSCode: doesn't seem to have a reliable extension to open md files as notebook, so you may have to manually pair them using jupytext in terminal.
>- Other editors: ...

## Updating an existing notebook (.md format)
- Open it in jupyter notebook or lab, jupytext extension should be installed (it's already listed in `requirements.txt` file mentioned above)
- Right click on a md notebook file and "open with" notebook format. 
- Edit it as you would a regular notebook and once saved the changes should 
reflect in md
- Make sure only content changes appear in md diff, not the frontmatter, before 
comitting
- Make a PR once ready (TODO: explain how like forking etc.)

## Adding a new notebook (from .ipynb format)
To add an ipynb notebook you were working on to this repository, convert that to
MyST markdown using:

```python
jupytext --to md:myst notebook.ipynb 
```

Then add it to `tutorials/` and also add its path in a relevant TOC in `index.md`.