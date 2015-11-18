# contact_map.py
A tool to compute and display contact maps between two PDB chains.

Makes a contact map between two chains (Defaults A and B) from a single PDB.
The area is white if there's no contact within the specified cutoff.
Also outputs a modified .pdb file, in case the submitted one has problems.
Outputs the separate models too

Usage: python ./contact_map.py -f file.pdb -o output_name [-d directory, -cmax cutoff_max, -cmin cutoff_min, -fc first_chain, -sc second_chain]
