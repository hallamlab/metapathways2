# MetaPathways v2.0: A master-worker model for environmental Pathway/Genome Database construction on grids and clouds

Niels W. Hanson, Kishori M. Konwar, Shang-Ju Wu, and Steven J. Hallam

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of infor- mation into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. We recently introduced MetaPathways v1.0, a modular annota- tion and analysis pipeline for constructing environmental Path- way/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

Here we present MetaPathways v2.0 incorporating a graphical user interface (GUI) and refined task management methods. The MetaPathways GUI provides an intuitive display for setup and process monitoring and supports interactive data visualization and sub-setting via a custom Knowledge Engine data structure. A master-worker model is adopted for task management allowing users to scavenge computational results from a number of worker grids in an ad hoc, asynchronous, distributed network that dramatically increases fault tolerance. This model facilitates the use of EC2 instances extending ePGDB construction to the Amazon Elastic Cloud.

## Installation

MetaPathways v2.0 requires Python 2.6 or greater and [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) developed by SRI International for full functionality.

The MetaPathways Python codebase as well as the GUI code is self-contained in this GitHub distro.

Please see the [MetaPathways v2.0 wiki](https://github.com/hallamlab/metapathways2/wiki) for more installation and usage information.

A template [MetaPathways_DB](https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip?dl=0) contains starter protein and taxonomic databases
