# PyAspect

PyAspect is an open-source package the helps setup and automate SPECFEM3D_Cartesian for both inversion and forward workflows.  (This package is still in development and a publication is in progress, so please site this github repository URL if you use PyAspect for your research). 


# Features
* Generate a regular mesh for SPECFEM3D_Cartesian from a 7D (x,y,z,VP,VS,Rho,Q) numpy array
* Generate a project for waveform inversion or foward modeling workflows, as well as tools for generating workflows that take advantange of reciprocity.
* Bash and SLURM scripts to submit a job-array for running multiple inversion iterations or a single forward iteration
* Pythonic tools for creating "SPECFEM-complient" source and receiver (station) files, including single or moment-tensor sources and multi-component receivers, and for creating headers (Python dictionaries) and header-records via pandas.DataFrame files for tracking both data and related "SPECFEM-complient" source and receiver files. Because the headers are Python dictionaries, they can have an unlimited number of header-fields (words) and lists of headers can be easily transformed in to header-records (pandas.DataFrame files) which can be sorted independently of the data.
* Tools for generating reciprocal sources and receivers from lists of standard source and receivers or from a standard header-record.
* Jupyter Notebook tools for using PyVista for 3D visualization of subsurface models (7D numpy arrays) and source and recievers plots
