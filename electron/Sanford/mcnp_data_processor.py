#! /usr/bin/env python
import datetime
import os
import shutil
import sys, getopt
from subprocess import call

def main(argv):
    directory = ''
    try:
        opts, args = getopt.getopt(argv,"hd:",["out_dir="])
    except getopt.GetoptError:
        print 'data_processor.py -d <directory>'
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print 'data_processor.py -d <directory>'
            sys.exit(1)
        elif opt in ("-d", "--out_dir"):
            directory = arg

    cell_list = ['10']
    surface_list = ['100', '101']
    angle_list = ['', 'full' ]

    # Get mcnp output file name
    base = "mcnp"
    mcnp_output = base+".o"

    # Check if file exists
    if os.path.isfile(mcnp_output):
        # Check if the output directory exists and make if necessary
        if not os.path.isdir(directory):
            print "Making directory",directory
            os.makedirs(directory)

        # Move file to output directory
        new_name = str(directory)+base
        shutil.move(mcnp_output,new_name+".o")
        shutil.move(base+".m",new_name+".m")
        shutil.move(base+".r",new_name+".r")

        # Move to output data directory
        os.chdir(directory)

        today = datetime.date.today()
        # Read the mcnp data file for surface tallys
        with open(mcnp_output) as data:
            # go through all surface tallies
            for i in cell_list:
                start=" cell  "+i
                name = base+"_cell_flux.txt"
                file = open(name, 'w')
                header = "# Energy   flux \t   Sigma\t"+str(today)+"\n"
                file.write(header)
                # Skips text before the beginning of the interesting block:
                for line in data:
                    if line.startswith(start):
                        data.next()
                        break
                # Reads text until the end of the block:
                for line in data:  # This keeps reading the file
                    if line.startswith('      total'):
                        file.close()
                        break
                    line = line.lstrip()
                    line = line.replace('   ',' ')
                    file.write(line)

        with open(mcnp_output) as data:
            # go through all surface tallies
            for i in surface_list:
                start=" surface  "+i

                # go through the current estimators first angle
                name = base+"_"+i+".txt"
                file = open(name, 'w')
                header = "# Angle     Current     Sigma\t"+str(today)+"\n"
                total_number_of_angles = 18
                number_of_angles = 0
                file.write(header)
                file.write("-1.00000000000000 0.00000E+00 0.0000\n")
                # Skips text before the beginning of the interesting block:
                for line in data:
                    if line.startswith(start):
                        line = data.next().strip()
                        if number_of_angles == 0:
                            line = "0.000000000000000" + line.replace('angle  bin:  -1.          to  0.00000E+00','')
                        elif number_of_angles == 1:
                            line = "0.939692620785908" + line.replace('angle  bin:   0.00000E+00 to  9.39693E-01','')
                        elif number_of_angles == 2:
                            line = "0.965925826289068" + line.replace('angle  bin:   9.39693E-01 to  9.65926E-01','')
                        elif number_of_angles == 3:
                            line = "0.984807753012208" + line.replace('angle  bin:   9.65926E-01 to  9.84808E-01','')
                        elif number_of_angles == 4:
                            line = "0.990268068741570" + line.replace('angle  bin:   9.84808E-01 to  9.90268E-01','')
                        elif number_of_angles == 5:
                            line = "0.994521895368273" + line.replace('angle  bin:   9.90268E-01 to  9.94522E-01','')
                        elif number_of_angles == 6:
                            line = "0.995396198367179" + line.replace('angle  bin:   9.94522E-01 to  9.95396E-01','')
                        elif number_of_angles == 7:
                            line = "0.996194698091746" + line.replace('angle  bin:   9.95396E-01 to  9.96195E-01','')
                        elif number_of_angles == 8:
                            line = "0.996917333733128" + line.replace('angle  bin:   9.96195E-01 to  9.96917E-01','')
                        elif number_of_angles == 9:
                            line = "0.997564050259824" + line.replace('angle  bin:   9.96917E-01 to  9.97564E-01','')
                        elif number_of_angles == 10:
                            line = "0.998134798421867" + line.replace('angle  bin:   9.97564E-01 to  9.98135E-01','')
                        elif number_of_angles == 11:
                            line = "0.998629534754574" + line.replace('angle  bin:   9.98135E-01 to  9.98630E-01','')
                        elif number_of_angles == 12:
                            line = "0.999048221581858" + line.replace('angle  bin:   9.98630E-01 to  9.99048E-01','')
                        elif number_of_angles == 13:
                            line = "0.999390827019096" + line.replace('angle  bin:   9.99048E-01 to  9.99391E-01','')
                        elif number_of_angles == 14:
                            line = "0.999657324975557" + line.replace('angle  bin:   9.99391E-01 to  9.99657E-01','')
                        elif number_of_angles == 15:
                            line = "0.999847695156391" + line.replace('angle  bin:   9.99657E-01 to  9.99848E-01','')
                        elif number_of_angles == 16:
                            line = "0.999961923064171" + line.replace('angle  bin:   9.99848E-01 to  9.99962E-01','')
                        elif number_of_angles == 17:
                            line = "1.000000000000000" + line.replace('angle  bin:   9.99962E-01 to  1.00000E+00','')
                        line = line.replace(' mu',' ')
                        line+=data.next().strip()+'\n'
                        file.write(line)
                        number_of_angles+=1
                        if number_of_angles == total_number_of_angles:
                            break
                file.close()

        # Plot results
#        plot = "../../plot_"+base+".p"
#        call(["gnuplot", plot])

    else:
        print "File ",mcnp_output," does not exist!"

if __name__ == "__main__":
   main(sys.argv[1:])
