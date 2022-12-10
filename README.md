# JA-graph

Source code for 'Affinity graph with elasticnet regularization for natural image segmentation'.

A list of papers and datasets about natural/color image segmentation can be found in our [repository](https://github.com/Yangzhangcst/Natural-color-image-segmentation).

JA-graph is modified from [the offcial SAS implementation](http://www.ee.columbia.edu/ln/dvmm/SuperPixelSeg/dlform.htm).


### Requirements
The code requires the version of Matlab2018a, Ubuntu 16.04.
(The code can also support the version of Matlab2018a+ in Windows, but the results in the paper may not be available.)

### Data
The original BSD500 dataset can be downloaded from [here](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html). We place the processed dataset in `database/BSDS500.zip`.


### Demo
Run the demo `JA_BSDS500_searchN.m`.
