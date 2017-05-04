import sys
import math
import argparse
from itertools import combinations


def dist(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)


def read_atoms(path):
    f = open(path, "r")
    atoms = []
    for line in f:
        line = line.strip()
        if line.startswith("ATOM"):
            type = line[12:16].strip()
            if type == "CA":
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append((x, y, z))
    f.close()
    return atoms


def compute_contacts(atoms, threshold):
    contacts = []
    for i in range(len(atoms)-2):
        for j in range(i+2, len(atoms)):
            if dist(atoms[i], atoms[j]) < threshold:
                contacts.append((i+1, j+1))
    return contacts


def write_output(contacts, path):
    f = open(path, "w")
    for c in contacts:
        f.write("\t".join(map(str, c))+"\n")
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Computes protein contact map from PDB file.")
    parser.add_argument("pdb_file", type=str, help="PDB input file.")
    parser.add_argument("output_file", type=str, help="File to write contact map to.")
    parser.add_argument("-t", "--threshold", type=float, required=False, default=7.5, help="Contact distance threshold in ångström.")
    args = parser.parse_args()

    atoms = read_atoms(args.pdb_file)
    contacts = compute_contacts(atoms, args.threshold)
    write_output(contacts, args.output_file)

if __name__ == '__main__':
    main()
