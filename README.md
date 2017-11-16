# A repository to explore semidefinite programs in quantum information

Right now, it contains two subprojects:

- `diagout`, a pseudo solver interface for CVX, useful to diagnose the problem size and dump the model to disk in SeDuMi format

- `doherty`, a study of three formulations of the Doherty hierarchy.

## Please contribute

Our goal is to collect as many instances of quantum information SDPs. 
Please help by providing additional instances; for that, propose your own folder in a PR.
In the folder, include the `.mat` files in the SeDuMi format, obtained by using `diagout` (CVX) or `export` (YALMIP).
Add a small README with your name, related publications if applicable and most importantly, the expected objective value.
If possible, include self-contained MATLAB code that generated the files in a subfolder of your folder.

In our study, priority will be given to problems satisfying the conditions below:

- Self-contained MATLAB code to generate the problem definitions is provided in a subfolder (see `doherty`). 
  Standard open source libraries such as QETLAB, CVX, YALMIP do not need to be included, but specify in the README
  what should be in the MATLAB path to run the code.

- The contributed problems are part of a series or hierarchy of problems with similar structure, 
  with increasing number of variables/constraints ranging from small to large scale,
  
- The problem has a known analytical solution,

- The coefficients in the problem definition are "simple", i.e. exact forms of the coefficients can be recovered easily (i.e. small fractions, or multiples of small square roots).
