## MERS-CoV recombination: implications about the reservoir and potential for adaptation
### Scripts

The following python scripts were used in the manuscript:
- `Host_association_calculator.py`
- `Hrothvitnir_MERS.py`
- `arlequin_to_fasta.py`
- `cfml_homoplasies.py`

Most figures were made using IPython notebooks:
- `MERS_mapHomoplasies.ipynb`
- `MERS_recombination.ipynb`

### Host_association_calculator.py
`Host_association_calculator.py` is a reworked and rough script that essentially calculates the commonly used linkage disequilibrium metrics to quantify the association between particular alleles in MERS-CoV genomes and whether the name of a sequences has the word "camel" in it.

### Hrothvitnir_MERS.py
`Hrothvitnir_MERS.py` is a descendent of a [script I used for a previous study](http://mbe.oxfordjournals.org/content/32/1/162.full), itself a descendent of a [short linked list script](http://stackoverflow.com/questions/280243/python-linked-list/280286#280286) from StackOverflow.
It is used to pull out posterior distributions of mutations inferred by robust counting from [posterior distributions of trees](https://code.google.com/p/beast-mcmc/).

### arlequin_to_fasta.py
`arlequin_to_fasta.py` is a very quick script I put together to convert fastsimcoal2-simulated sequences from arlequin to fasta format.

### cfml_homoplasies.py
`cfml_homoplasies.py` is a script put together and kindly donated by [Jessica Hedge](https://scholar.google.co.uk/citations?user=u3bsLoQAAAAJ&hl=en&oi=ao) (University of Oxford).
After minor modification this particular version extracts all asymmetric mutations inferred by ClonalFrameML for homoplasy analyses.


## Figures
`MERS_recombination.ipynb` makes most of the figures in the manuscript. **IPython**, **matplotlib** and **numpy** should be the only libraries required to make the figures. You can have a look at the contents of the notebook [here](http://nbviewer.ipython.org/github/evogytis/MERS_recombination/blob/master/scripts/MERS_recombination.ipynb).
`MERS_mapHomoplasies.ipynb` is a specialised IPython notebook that plots phylogenetic trees, maps mutations to individual branches and connects branches with identical mutations inferred by ClonalFrameML. The notebook's contents can be viewed [here](http://nbviewer.ipython.org/github/evogytis/MERS_recombination/blob/master/scripts/MERS_mapHomoplasies.ipynb).


==============
Please keep in mind that these scripts are fairly rough, even though I've tried to comment the code as much as possible. Expect some tedium if you try to run these yourself.