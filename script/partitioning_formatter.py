"""partitioning_formatter.py -- pretty-print MFDn partitioning files.

Patrick Fasano
University of Notre Dame

- 12/29/17 (pjf): Created."""

import argparse

parser = argparse.ArgumentParser(
    description="Pretty-print MFDn partitioning files.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input',
                    default="mfdn_partitioning.generated",
                    help="input (generated) partitioning file")
parser.add_argument('--orbitals',
                    default="orbitals.dat",
                    help="orbital file")
parser.add_argument('--output',
                    default="mfdn_partitioning.info",
                    help="ouput partitioning file")
args = parser.parse_args()

weights = []
line_number = 1
with open(args.orbitals) as fp:
    for line in fp:
        line = line.rstrip()
        line = line.lstrip()
        if line[0] == '#':
            continue
        line_tokens = line.split()
        if line_number == 1:
            assert int(line_tokens[0]) == 15099
        elif line_number == 2:
            (num_p, num_n) = line_tokens
        else:
            weights.append(line_tokens[-1])
        line_number += 1

partitioning_generated = []
header_line = ""
with open(args.input) as fp:
    header_line = fp.readline()
    partitioning_generated = fp.read().split()
    fp.close()

if len(weights) != len(partitioning_generated):
    raise ValueError

with open(args.output, 'w') as fp:
    fp.write(header_line)
    for i, partition in enumerate(partitioning_generated):
        if weights[i-1] != weights[i]:
            fp.write("\n")
        fp.write("{:>8}".format(partitioning_generated[i]))
    fp.write("\n")
    fp.close()
