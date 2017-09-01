# Zero-tree Orthogonal Matching Pursuit

Orthogonal Matching Pursuit inspired by wavelet zero-trees. This repository
contains scripts to test zero-tree OMP and compare it against traditional OMP.

## Content
1. **tools/** contains code for running zero-tree OMP, traditional OMP and other
   scripts to train dictionaries. There are also some mex functions that have
   been compiled against BLAS and LAPACK. If they do not work, please contact
   me for further assistance.
2. **data/** should contain dictionaries and sample video for simulation. See
   instructions below for download the folder
3. **test-video-cs.m** enables you to test compressive sensing for a video.
4. **train-video-dictionaries.m** trains video dictionaries.
5. **train-dictionaries.m** is used as a general training script.

## Data
Dowload the following files and place them in _data/_ folder:
1. [Single scale dictionary](https://cmu.box.com/s/40exxu3cb8lyopm2qba0uwxcicue3xly)
2. [Two scale dictionary](https://cmu.box.com/s/rbvjptnzmtolkw9ysgt8equifduiyxrw)
3. [Test video](https://cmu.box.com/s/vaox6vzzvut201l69gwk7yhekd0v7ik4)

## Notes
1. The folder 'ksvdbox' is by [Ron Rubenstein](http://www.cs.technion.ac.il/~ronrubin/software.html)
2. Mex routines with BLAS/LAPACK will give you accurate timing profiles. However
   if you just want to test, you may equivalently use the matlab scripts.
   Example, use **comp_cs.m** instead of **comp_cs_c**
3. *tools/lib* has certain BLAS and LAPACK libraries compiled for a very
   specific 64bit computer. If your matlab crashes, please compile your own
   BLAS/LAPACK. For help, please contact me.
4. Only sample dictionaries and a sample video is provided. If you want to train
   your own dictionaries, you will need your own data.
5. Test video source is [here](http://high-speed-video.colostate.edu/)

## Citation
Get paper details [here](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=u-xGD2AAAAAJ&citation_for_view=u-xGD2AAAAAJ:YOwf2qJgpHMC)
