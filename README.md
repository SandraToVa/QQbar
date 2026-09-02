# Hyperfine Splitting scripts
Matlab scripts (quarkonium and hybrids)
- **RubenFiles**: are not modified from the code that Joan gave me from Ruben
- **SandraFiles**: the files that I have created
    - <ins>Spectrum</ins>: in this folder I have all the scripts that are used to compute the spectrum and obtain the correct coeficients for the interpolations
    - <ins>Decays</ins>: in this folder I have all related to the computation of the decay with between Q->Q and H->Q
    - <ins>Practiques</ins>: first scripts, contains the intrepolation script and other practice ones
    - <ins>Spectrum Update</ins>: is the same logic as the <ins>Spectrum</ins> folder but in here we use the new lattice results for the potentials instead of just phenomenological fits. Also, the structure of the codes has changed to imporve understanding.
- **Source**: the files where the computation of the schrödinger equations are done. They compute the actual wavefunctions and spectrum for a general case
- The srcipts withou folder are the actual interesting part. For each quarkonium and hybrid state they compute the spectrum and the wavefunctions using the files from _Source_. The other files all use these scripts
