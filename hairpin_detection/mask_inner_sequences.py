#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv[1:]) != 2:
        print("Usage: {} <fasta file> <end_length>".format(sys.argv[0]))
        sys.exit()

    file_in = sys.argv[1]
    end_len = int(sys.argv[2])

    header = None

    with open(file_in, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line[0] == ">":
                header = line[1:]
            else:
                if len(line) <= 2*end_len:
                    print(">{}\n{}".format(header, line))
                else:
                    seq_len = len(line)
                    new_seq = line[:end_len] + "N"*len(line[end_len:seq_len-end_len]) + line[seq_len-end_len:]
                    print(">{}\n{}".format(header, new_seq))



if __name__ == "__main__":
    main()