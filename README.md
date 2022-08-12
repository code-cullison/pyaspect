[![DOI](https://zenodo.org/badge/386649496.svg)](https://zenodo.org/badge/latestdoi/386649496)

# PyAspect

PyAspect is an open-source package to assist with automating SPECFEM3D_Cartesian workflows for both waveform inversion and forward modeling.  (This package is still in development and a publication is in progress, so please site this GitHub repository URL if you use PyAspect for your research). 


# Features
* Generate a regular mesh for SPECFEM3D_Cartesian from a 7D (x,y,z,VP,VS,Rho,Q) numpy array
* Generate a project for waveform inversion or forward modeling workflows, as well as tools for generating workflows that take advantage of reciprocity.
* Bash and SLURM scripts to submit a job-array for running multiple inversion iterations or a single forward iteration
* Pythonic tools for creating "SPECFEM-compliant" source and receiver (station) files, including single or moment-tensor sources and multi-component receivers, and for creating headers (Python dictionaries) and header-records via pandas.DataFrame files for tracking both data and related "SPECFEM-compliant" source and receiver files. Because the headers are Python dictionaries, they can have an unlimited number of header fields (words) and lists of headers can be easily transformed into header-records (pandas.DataFrame files) which can be sorted independently of the data.
* Tools for generating reciprocal sources and receivers from lists of standard sources and receivers or from a standard header-record.
* Jupyter Notebook tools for using PyVista for 3D visualization of subsurface models (7D NumPy arrays) and source and receivers plots.
* Note some features in this package are the same or similar to [gnam](https://github.com/code-cullison/gnam) because PyAspect diverged from [gnam](https://github.com/code-cullison/gnam) at an early stage.


### Contributors

- Thomas Cullison:
  - Primary author
  - Tester
- La Ode Marzujriban Masfara:
  - Co-designer of reciprocal moment tensor construction
  - Tester
- Rhys Hawkins:
  - Contributed to design discussions
