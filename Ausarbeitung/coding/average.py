#!/usr/bin/env python

import argparse

def average(inputfile):
  total = count = 0

  with open(inputfile, 'r') as inp:
    for line in inp:
      try:
        count += 1
        num = float(line)
        total += num
      except ValueError:
	    print('{} is not a number!'.format(line))

    print('Total of all numbers: {}'.format(total))
    print "Avg:", int(total)/count

def main():

  parser = argparse.ArgumentParser()

  parser.add_argument("-f", "--inputfile", help="Input File")

  args = parser.parse_args()

  inputfile = ""
  if args.inputfile:
    average(args.inputfile)
if __name__ == "__main__":
  main()
