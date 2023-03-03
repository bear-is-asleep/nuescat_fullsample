# Neutrino Electron Scattering for Flux Constraint

https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=29712
https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=29199
https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=28114
https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=27103

## caf_scripts/mre
This code makes use CAFana macros used to extract data, extract systematics using weights, and produce plots. The common way to execute any of the caf analyzers is to call `cafe -b <CAFanalyzer.C>`. The files in `mre` are meant for debugging using simple cuts and plots. They're also designed to execute quickly. The files in `caf_scripts` are primarily used to for everything interesting. 



Use CAFana for full reconstruction
