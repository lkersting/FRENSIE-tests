#! /usr/bin/env python
from os import path, makedirs
import sys
import argparse as ap

# Add the parent directory to the path
parent_dir=path.dirname(path.dirname(path.abspath(__file__)))
sys.path.insert(1,parent_dir)

import albedo_simulation as simulation

# Set up the argument parser
description = "This script asks for albedo data and run names which "\
              "which it then plots against experimental data."

parser = ap.ArgumentParser(description=description)

parser.add_argument("rendezvous_files", nargs='*', help="the rendezvous file to be processed")

# Parse the user's arguments
user_args = parser.parse_args()

for rendezvous in user_args.rendezvous_files:
  simulation.processData( rendezvous )
