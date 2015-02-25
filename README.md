# MetaPathways 2: A master-worker model for environmental Pathway/Genome Database construction on grids and clouds

Niels W. Hanson, Kishori M. Konwar, Shang-Ju Wu, and Steven J. Hallam

## Updates

**November 27, 2014**: [MetaPathways v2.5 released](https://github.com/hallamlab/metapathways2/releases/tag/v2.5) with upgrades to the pipeline:
    
* LAST homology searches with BLAST-equivalent output and E-values
* Reads per kilobase per million mapped (RPKM) coverage measure for Contig annotations calculated from raw reads (`.fastq`) or mapping files (`.SAM`) using [bwa](http://bio-bwa.sourceforge.net)
* Addition of the [CAZy sequence database](http://www.cazy.org) as a new compatible functional hierachy
* GUI Keyword-search from annotation subsetting and projection onto different functional hierarcies (KEGG, COG, SEED, MetaCyc, and now CAZy)

See [the release page](https://github.com/hallamlab/metapathways2/releases/tag/v2.5) and [the wiki](https://github.com/hallamlab/metapathways2/wiki) for more information.


## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of infor- mation into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. We recently introduced MetaPathways v1.0, a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

Here we present MetaPathways v2.0 incorporating a graphical user interface (GUI) and refined task management methods. The MetaPathways GUI provides an intuitive display for setup and process monitoring and supports interactive data visualization and sub-setting via a custom Knowledge Engine data structure. A master-worker model is adopted for task management allowing users to scavenge computational results from a number of worker grids in an ad hoc, asynchronous, distributed network that dramatically increases fault tolerance. This model facilitates the use of EC2 instances extending ePGDB construction to the Amazon Elastic Cloud.

## Installation

MetaPathways v2.5 requires Python 2.7 or greater and [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) developed by SRI International for full functionality.

The MetaPathways Python codebase as well as the compiled GUI binaries for Mac OSX and Ubuntu are self-contained in this GitHub distro. GUI source code can be [obtained here](https://github.com/hallamlab/MetaPathwaysGUI).

Please see the [MetaPathways v2.5 wiki](https://github.com/hallamlab/metapathways2/wiki) for more installation details.

A template [MetaPathways_DBs.zip (**Updated: October 2014**)](https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip?dl=0) contains starter protein and taxonomic databases

## Citation

If using MetaPathways v2.0 for reserach work please cite:

Niels W. Hanson, Kishori M. Konwar, Shang-Ju Wu, Steven J. Hallam. *MetaPathways v2.0: A master-worker model for environmental Pathway/Genome Database construction on grids and clouds.* Proceedings of the 2014 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2014), Honolulu, HI, USA, May 21-24, 2014. [doi:10.1109/CIBCB.2014.6845516](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6845516)
