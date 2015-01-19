## Does the Middle East Respiratory Syndrome coronavirus (MERS-CoV) recombine?

Gytis Dudas<sup>1</sup> and Andrew Rambaut<sup>1,2,3</sup>

<sup>1</sup>Institute of Evolutionary Biology, University of Edinburgh, Edinburgh, UK,
<sup>2</sup>Fogarty International Center, National Institutes of Health, Bethesda, MD, USA,
<sup>3</sup>Centre for Immunology, Infection and Evolution at the University of Edinburgh, Edinburgh, UK

### Quick summary
Very likely, yes.
This repo is what I've been up to lately and is nearly finished. Most/all of the work has been done and it shows the presence of the classical recombination triad in the MERS genome: excessive homoplasies, decay of linkage disequilibrium and the presence of alternative topologies.

Done:
- GARD analyses. Show that the MERS genome has well-supported alternative tree topologies and some degree of rate heterogeneity.
- LDhat analyses. Permutation tests for recombination in LDhat appear quite robust to temporal sampling and indicate LD decay. Some other tests don't do so well.
- Homoplasy analyses. Ancestral sequence reconstruction using ClonalFrameML indicates excessive numbers of homoplasies in the MERS genome.
- Host-polymorphism association analysis. Present data do not show the presence of "human" or "camel" alleles, but are also too poorly sampled to address the question in great detail.
- Homoplasy rates integrating over possible tree topologies. Shows almost identical results to ML analysis.
- Marginal likelihood estimates of models including/excluding rate heterogeneity and/or alternative tree topologies. Mostly done and appears to support a model with rate heterogeneity alone, which is not surprising given the overall appearance of recombination and the cost of the number of parameters introduced by another phylogenetic tree.
- Manuscript.
