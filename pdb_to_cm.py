import sys
import math
import argparse
import re
from itertools import combinations


def dist(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)


def read_atoms(file, chain="."):
    pattern = re.compile(chain)

    atoms = []
    for line in file:
        line = line.strip()
        if line.startswith("ATOM"):
            type = line[12:16].strip()
            chain = line[21:22]
            if type == "CA" and re.match(pattern, chain):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append((x, y, z))
    return atoms


def compute_contacts(atoms, threshold):
    contacts = []
    for i in range(len(atoms)-2):
        for j in range(i+2, len(atoms)):
            if dist(atoms[i], atoms[j]) < threshold:
                contacts.append((i+1, j+1))
    return contacts


def write_output(contacts, file):
    for c in contacts:
        file.write("\t".join(map(str, c))+"\n")


def pdb_to_cm(file, threshold, chain="."):
    atoms = read_atoms(file, chain)
    return compute_contacts(atoms, threshold)


def main():
    parser = argparse.ArgumentParser(description="Computes protein contact map from PDB file.")
    parser.add_argument("pdb_file", type=str, help="PDB input file.")
    parser.add_argument("output_file", type=str, help="File to write contact map to.")
    parser.add_argument("-c", "--chain", type=str, required=False, default=".", help="Chain ID to use. Supports regular expression.")
    parser.add_argument("-t", "--threshold", type=float, required=False, default=7.5, help="Contact distance threshold in angstrom.")
    args = parser.parse_args()

    contacts = pdb_to_cm(open(args.pdb_file, "r"), args.threshold, args.chain)
    write_output(contacts, open(args.output_file, "w"))

if __name__ == '__main__':
    main()
