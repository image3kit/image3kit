svplot 
========================================================================

Warning: this library needs a lot of clean up and testing.


Easily plot C++ data.  With few lines of C++ code, the svplot C++ library allows efficient plotting of data held in STL containers, like vector and map,  as quality Scalable Vector Graphic (.svg) files, which in turn can be converted in other image formats.


Several 2D plot variants can be produced, and there are a myriad of options to control appearance.

This is used as the standard plotting and output format for two-dimensional (including projections of three-dimensional data) in other modelling and data analysis C++ codes. For instance all plots in <https://doi.org/10.1007/s11242-019-01317-8> are generated using the original version of this library.  


SVG allows embedding of plot data as textual information. This adds the flexibility to save and exchange raw plot data and their representation using a single text (.svg) file.  


## History

The svg_plot project was originally written by Jake Voytko in 2007 as a Boost-sponsored Google Summer of Code project in 2007. It has been maintained and enhanced since then in Boost Sandbox, but is judged unsuitable for a Boost Library, so it was made available on <https://github.com/pabristow/svg_plot>.  

In 2018-209, Ali Q. Raeini, have simplified the svg_plot and renamed it to svplot. In this process, some advanced functionalities, in particular built-in uncertainty handling, and dependency on all other Boost libraries are removed. Therefore, svplot has no external dependency, and only depend on standard template library (STL) and C++11 or newer C++ compiler.   
New functionalities are added, so that additional informations can be visualized by assigning varying sizes or colours to plot data points.  
For simplicity, all the plots are generated as multi-plots, which by default contain a single plot reproducing the original behaviour.


## Usage

For documentation, and a demo example, see the [example] folder.


The documentation is incomplete and does not fully reflect the  changes 
made in version 2018-2019, consider having a look into source code if anything fails.

You can use Inkscape to convert the generated svg files into other image formats.  
For reference, to convert to pdf, you can run the command:

`inkscape --without-gui --export-pdf="FileName.pdf" "FileName.svg"`

## TODO
- The code internals need further simplification and testing.  
- X-axis and Y-axis implementations should be merged into one templated function. 
- Add more tests file and release on Github.
- Update documentation.

## Copyright 
- Copyright Paul A. Bristow 2013, 2014
- Copyright Ali Q. Raeini 2019
- Distributed under the [Boost Software License, Version 1.0](http://www.boost.org/LICENSE_1_0.txt)


