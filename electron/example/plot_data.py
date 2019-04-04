#!/usr/bin/python

import sys, getopt
import math as m
import matplotlib.pyplot as plt

def main(argv):
    mcnpinputfile = ''
    aceinputfile = ''
    endlinputfile = ''
    momentsinputfile = ''
    outputname = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:a:e:",["mfile=","afile=","efile=","ifile=","ofile="])
    except getopt.GetoptError:
        print 'plot_data.py -m <mcnpinputfile> -a <aceinputfile> -e <endlinputfile> -i <momentsinputfile> -o <outputname>'
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print 'plot_data.py -m <mcnpinputfile> -a <aceinputfile> -e <endlinputfile> -i <momentsinputfile> -o <outputname>'
            sys.exit(1)
        elif opt in ("-m", "--mfile"):
            mcnpinputfile = arg
        elif opt in ("-a", "--afile"):
            aceinputfile = arg
        elif opt in ("-e", "--efile"):
            endlinputfile = arg
        elif opt in ("-i", "--ifile"):
            momentsinputfile = arg
        elif opt in ("-o", "--ofile"):
            outputname = arg

    outputfile = outputname+'.txt'
    output_c_e = outputname+'_c_e.txt'

    # open the ace file
    f = open(aceinputfile, 'r')

    words = []
    ace_energy = []
    ace_value = []
    ace_error = []
    for line in f:
        words = line.split()
        if words[0] != "#":
            ace_energy.append( float(words[0]) )
            ace_value.append( float(words[1]) )
            ace_error.append( float(words[2]) )

    ace_error.pop(0)

    f.close()


    # open the endl file
    f = open(endlinputfile, 'r')

    words = []
    endl_energy = []
    endl_value = []
    endl_error = []
    for line in f:
        words = line.split()
        if words[0] != "#":
            endl_energy.append( float(words[0]) )
            endl_value.append( float(words[1]) )
            endl_error.append( float(words[2]) )

    endl_error.pop(0)

    f.close()

    words = []
    moments_energy = []
    moments_value = []
    moments_error = []

    print "\n\n",momentsinputfile,"\n\n"
    # open the endl file
    if momentsinputfile != '':
        print "\n\n",momentsinputfile,"\n\n"
        f = open(momentsinputfile, 'r')
        for line in f:
            words = line.split()
            if words[0] != "#":
                moments_energy.append( float(words[0]) )
                moments_value.append( float(words[1]) )
                moments_error.append( float(words[2]) )

        moments_error.pop(0)

        f.close()

    # open the mcnp file
    f = open(mcnpinputfile, 'r')

    words = []
    mcnp_energy = []
    mcnp_value = []
    mcnp_error = []
    for line in f:
        words = line.split()
        if words[0] != '#':
            mcnp_energy.append( float(words[0]) )
            mcnp_value.append( float(words[1]) )
            mcnp_error.append( float(words[2]) )

    mcnp_error.pop(0)

    f.close()

    # Calculate the average value in the bin
    energy = []
    ace_average_values = []
    endl_average_values = []
    moments_average_values = []
    mcnp_average_values = []
    for i in range(1,len(ace_energy)):
        energy.append( (ace_energy[i-1] + ace_energy[i])/2 )
        ace_average_values.append( ace_value[i]/(ace_energy[i]-ace_energy[i-1]) )
        endl_average_values.append( endl_value[i]/(endl_energy[i]-endl_energy[i-1]) )
        mcnp_average_values.append( mcnp_value[i]/(mcnp_energy[i]-mcnp_energy[i-1]) )

        if len(moments_energy) > 0:
            moments_average_values.append( moments_value[i]/(moments_energy[i]-moments_energy[i-1]) )

    # Calculate the propagated uncertainty
    ace_prop_uncert = []
    endl_prop_uncert = []
    moments_prop_uncert = []

    for i in range(0,len(energy)):
        if mcnp_average_values[i] > 0:
            ace_prop_uncert.append( m.sqrt( ((1.0/mcnp_average_values[i])**2)*(ace_average_values[i]*ace_error[i])**2 + ((ace_average_values[i]/mcnp_average_values[i]**2)**2)*(mcnp_average_values[i]*mcnp_error[i])**2 ) )
            endl_prop_uncert.append( m.sqrt( ((1.0/mcnp_average_values[i])**2)*(endl_average_values[i]*endl_error[i])**2 + ((endl_average_values[i]/mcnp_average_values[i]**2)**2)*(mcnp_average_values[i]*mcnp_error[i])**2 ) )
            if len(moments_energy) > 0:
                moments_prop_uncert.append( m.sqrt( ((1.0/mcnp_average_values[i])**2)*(moments_average_values[i]*moments_error[i])**2 + ((moments_average_values[i]/mcnp_average_values[i]**2)**2)*(mcnp_average_values[i]*mcnp_error[i])**2 ) )
        else:
            ace_prop_uncert.append( 0.0 )
            endl_prop_uncert.append( 0.0 )
            if len(moments_energy) > 0:
                moments_prop_uncert.append( 0.0 )

    # Create the output file
    f = open(outputfile, 'w')
    header = ''
    if len(moments_energy) > 0:
        header="#MCNP\tFACEMC-ACE\tFACEMC-ENDL(LinLin)\tFACEMC-Moments\n"
    else:
        header="#MCNP\tFACEMC-ACE\tFACEMC-ENDL(LinLin)\n"
    f.write( header )
    f.write( "#Energy\tValues\tError\n" )

    for i in range(0,len(energy)):
        line = ' '
        line = line + str(energy[i])
        line = line + ' '
        line = line + str(mcnp_average_values[i])
        line = line + ' '
        line = line + str(mcnp_error[i])
        line = line + ' '
        line = line + str(ace_average_values[i])
        line = line + ' '
        line = line + str(ace_error[i])
        line = line + ' '
        line = line + str(endl_average_values[i])
        line = line + ' '
        line = line + str(endl_error[i])
        if len(moments_energy) > 0:
            line = line + ' '
            line = line + str(moments_error[i])
            line = line + ' '
            line = line + str(moments_error[i])
        line = line + '\n'

        f.write( line )

    # Calculate c/e values
    ace_c_e_ratio = [None]*len(energy)
    ace_c_e_uncert = [None]*len(energy)
    endl_c_e_ratio = [None]*len(energy)
    endl_c_e_uncert = [None]*len(energy)
    moments_c_e_ratio = [None]*len(energy)
    moments_c_e_uncert = [None]*len(energy)
    print "----------------------------------------------------------------"
    if len(moments_energy) > 0:
        print "Energy\tFACEMC-ACE\tFACEMC-EDNL(LinLin)\tFACEMC-Moments\tMCNP\tC/E Ratio (ACE)\tC/E Uncertainty(ACE)\tC/E Ratio (ENDL)\tC/E Uncertainty(ENDL)\tC/E Ratio (Moments)\tC/E Uncertainty(Moments)"
    else:
        print "Energy\tFACEMC-ACE\tFACEMC-EDNL(LinLin)\tMCNP\tC/E Ratio (ACE)\tC/E Uncertainty(ACE)\tC/E Ratio (ENDL)\tC/E Uncertainty(ENDL)"
    print "----------------------------------------------------------------"
    for i in range(0,len(energy)):
        if mcnp_average_values[i] > 0:
            if ace_average_values[i] > 0:
                ace_c_e_ratio[i] = ( ace_average_values[i]/mcnp_average_values[i] )
                ace_c_e_uncert[i] = ( ace_prop_uncert[i] )
            elif mcnp_average_values[i] == ace_average_values[i]:
                ace_c_e_ratio[i] = 1.0
                ace_c_e_uncert[i] = ( ace_prop_uncert[i] )
            else:
                ace_c_e_ratio[i] = 0.0
                ace_c_e_uncert[i] = ( ace_prop_uncert[i] )

            if endl_average_values[i] > 0:
                endl_c_e_ratio[i] = ( endl_average_values[i]/mcnp_average_values[i] )
                endl_c_e_uncert[i] = ( endl_prop_uncert[i] )
            elif mcnp_average_values[i] == endl_average_values[i]:
                endl_c_e_ratio[i] = 1.0
                endl_c_e_uncert[i] = ( endl_prop_uncert[i] )
            else:
                endl_c_e_ratio[i] = 0.0
                endl_c_e_uncert[i] = ( endl_prop_uncert[i] )

            if len(moments_energy) > 0:
                if moments_average_values[i] > 0:
                    moments_c_e_ratio[i] = ( moments_average_values[i]/mcnp_average_values[i] )
                    moments_c_e_uncert[i] = ( moments_prop_uncert[i] )
                elif mcnp_average_values[i] == moments_average_values[i]:
                    moments_c_e_ratio[i] = 1.0
                    moments_c_e_uncert[i] = ( moments_prop_uncert[i] )
                else:
                    moments_c_e_ratio[i] = 0.0
                    moments_c_e_uncert[i] = ( moments_prop_uncert[i] )

        if len(moments_energy) > 0:
            print energy[i], ace_average_values[i], endl_average_values[i], moments_average_values[i], mcnp_average_values[i], ace_c_e_ratio[i], ace_c_e_uncert[i], endl_c_e_ratio[i], endl_c_e_uncert[i], moments_c_e_ratio[i], moments_c_e_uncert[i]
        else:
            print energy[i], ace_average_values[i], endl_average_values[i],mcnp_average_values[i], ace_c_e_ratio[i], ace_c_e_uncert[i], endl_c_e_ratio[i], endl_c_e_uncert[i]

    # calculate % of c/e values within 1,2,3 sigma
    num_ace_c_e_in_one_sigma = 0
    num_ace_c_e_in_two_sigma = 0
    num_ace_c_e_in_three_sigma = 0
    num_endl_c_e_in_one_sigma = 0
    num_endl_c_e_in_two_sigma = 0
    num_endl_c_e_in_three_sigma = 0
    num_moments_c_e_in_one_sigma = 0
    num_moments_c_e_in_two_sigma = 0
    num_moments_c_e_in_three_sigma = 0
    num_below_ace = 0
    num_above_ace = 0
    num_below_endl = 0
    num_above_endl = 0
    num_below_moments = 0
    num_above_moments = 0

    start = 0
    for i in range(start,len(ace_c_e_ratio)):
        if ace_c_e_ratio[i] < 1.0:
            num_below_ace = num_below_ace+1
        else:
            num_above_ace = num_above_ace+1

        if endl_c_e_ratio[i] < 1.0:
            num_below_endl = num_below_endl+1
        else:
            num_above_endl = num_above_endl+1

        if len(moments_energy) > 0:
            if moments_c_e_ratio[i] < 1.0:
                num_below_moments = num_below_moments+1
            else:
                num_above_moments = num_above_moments+1

        diff = abs( 1.0 - ace_c_e_ratio[i] )
        #print ace_c_e_ratio[i], diff, 2*ace_c_e_uncert[i]
        if diff <= ace_c_e_uncert[i]:
            num_ace_c_e_in_one_sigma = num_ace_c_e_in_one_sigma + 1
        if diff <= 2*ace_c_e_uncert[i]:
            num_ace_c_e_in_two_sigma = num_ace_c_e_in_two_sigma + 1
        if diff <= 3*ace_c_e_uncert[i]:
            num_ace_c_e_in_three_sigma = num_ace_c_e_in_three_sigma + 1

        diff = abs( 1.0 - endl_c_e_ratio[i] )
        #print endl_c_e_ratio[i], diff, 2*endl_c_e_uncert[i]
        if diff <= endl_c_e_uncert[i]:
            num_endl_c_e_in_one_sigma = num_endl_c_e_in_one_sigma + 1
        if diff <= 2*endl_c_e_uncert[i]:
            num_endl_c_e_in_two_sigma = num_endl_c_e_in_two_sigma + 1
        if diff <= 3*endl_c_e_uncert[i]:
            num_endl_c_e_in_three_sigma = num_endl_c_e_in_three_sigma + 1

        if len(moments_energy) > 0:
            diff = abs( 1.0 - moments_c_e_ratio[i] )
            #print moments_c_e_ratio[i], diff, 2*moments_c_e_uncert[i]
            if diff <= moments_c_e_uncert[i]:
                num_moments_c_e_in_one_sigma = num_moments_c_e_in_one_sigma + 1
            if diff <= 2*moments_c_e_uncert[i]:
                num_moments_c_e_in_two_sigma = num_moments_c_e_in_two_sigma + 1
            if diff <= 3*moments_c_e_uncert[i]:
                num_moments_c_e_in_three_sigma = num_moments_c_e_in_three_sigma + 1

    print "----------------------------------------------------------------"
    print "% C/E in 1 sigma: "
    print "ACE:",num_ace_c_e_in_one_sigma, len(range(start,len(ace_c_e_ratio))), float(num_ace_c_e_in_one_sigma)/len(range(start,len(ace_c_e_ratio)))
    print "ENDL:",num_endl_c_e_in_one_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_one_sigma)/len(range(start,len(endl_c_e_ratio)))
    if len(moments_energy) > 0:
        print "ENDL:",num_endl_c_e_in_one_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_one_sigma)/len(range(start,len(endl_c_e_ratio)))
    print "----------------------------------------------------------------"
    print "% C/E in 2 sigma: "
    print "ACE:",num_ace_c_e_in_two_sigma, len(range(start,len(ace_c_e_ratio))), float(num_ace_c_e_in_two_sigma)/len(range(start,len(ace_c_e_ratio)))
    print "ENDL:",num_endl_c_e_in_two_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_two_sigma)/len(range(start,len(endl_c_e_ratio)))
    if len(moments_energy) > 0:
        print "ENDL:",num_endl_c_e_in_two_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_two_sigma)/len(range(start,len(endl_c_e_ratio)))
    print "----------------------------------------------------------------"
    print "% C/E in 3 sigma: "
    print "ACE:",num_ace_c_e_in_three_sigma, len(range(start,len(ace_c_e_ratio))), float(num_ace_c_e_in_three_sigma)/len(range(start,len(ace_c_e_ratio)))
    print "ENDL:",num_endl_c_e_in_three_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_three_sigma)/len(range(start,len(endl_c_e_ratio)))
    if len(moments_energy) > 0:
        print "ENDL:",num_endl_c_e_in_three_sigma, len(range(start,len(endl_c_e_ratio))), float(num_endl_c_e_in_three_sigma)/len(range(start,len(endl_c_e_ratio)))
    print "----------------------------------------------------------------"
    print "ACE % below: ", float(num_below_ace)/len(range(start,len(ace_c_e_ratio)))
    print "ACE % above: ", float(num_above_ace)/len(range(start,len(ace_c_e_ratio)))
    print "ENDL % below: ", float(num_below_endl)/len(range(start,len(endl_c_e_ratio)))
    print "ENDL % above: ", float(num_above_endl)/len(range(start,len(endl_c_e_ratio)))
    if len(moments_energy) > 0:
        print "Moments % below: ", float(num_below_moments)/len(range(start,len(moments_c_e_ratio)))
        print "Moments % above: ", float(num_above_moments)/len(range(start,len(moments_c_e_ratio)))

    # Create the output file for C/E data
    f = open(output_c_e, 'w')
    header = ''
    if len(moments_energy) > 0:
        header="#FACEMC/MCNP\tFACEMC-ENDL(LinLin)/MCNP\tFACEMC-Moments/MCNP\n"
    else:
        header="#FACEMC/MCNP\tFACEMC-ENDL(LinLin)/MCNP\n"
    f.write( header )
    f.write( "#Energy    Values    Error\n" )

    for i in range(0,len(energy)):
        line = ' '
        line = line + str(energy[i])
        line = line + ' '
        line = line + str(ace_c_e_ratio[i])
        line = line + ' '
        line = line + str(ace_c_e_uncert[i])
        line = line + ' '
        line = line + str(endl_c_e_ratio[i])
        line = line + ' '
        line = line + str(endl_c_e_uncert[i])
        if len(moments_energy) > 0:
            line = line + ' '
            line = line + str(moments_c_e_ratio[i])
            line = line + ' '
            line = line + str(moments_c_e_uncert[i])
        line = line + '\n'

        f.write( line )

    # Plot data
    fig1 = plt.figure(num=1, figsize=(10,5))
    plt.subplot2grid((2,6),(0, 0), colspan=6)
#    plt.xlabel('Energy (MeV)')
    plt.ylabel('Flux ($\#/cm^2$)', fontsize=15)
    plt.title('Comparison of the Track Length Flux in a Sphere of H', fontsize=16)
#    plt.xlim(0.85,1.0)
    plt.ylim(0.0,800000.0)

    plt.plot( energy, ace_average_values, label='FACEMC-ACE')
    plt.plot( energy, endl_average_values, label='FACEMC-ENDL')
    if len(moments_energy) > 0:
        plt.plot( energy, moments_average_values, label='FACEMC-Hybrid')
    plt.plot( energy, mcnp_average_values, label='MCNP')
    # plot with errorbars
    # plt.errorbar( energy, mcnp_average_values, yerr=mcnp_error, label='MCNP')
    # plt.errorbar( energy, ace_average_values, yerr=ace_error, label='FACEMC-ACE')
    # plt.errorbar( energy, endl_average_values, yerr=endl_error, label='FACEMC-ENDL')

    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    plt.legend(loc=2)
    plt.subplot2grid((2,6),(1, 0), colspan=6)
    plt.xlabel('Energy (MeV)', fontsize=15)
    plt.ylabel('C/E', fontsize=15)
#    plt.ylim(0.9, 1.1)
#    plt.plot( energy, ace_c_e_ratio )
    # plot with errorbars
    plt.errorbar( energy, ace_c_e_ratio, yerr=ace_c_e_uncert, label='FACEMC-ACE')
    plt.errorbar( energy, endl_c_e_ratio, yerr=endl_c_e_uncert, label='FACEMC-ENDL')
    if len(moments_energy) > 0:
        plt.errorbar( energy, moments_c_e_ratio, yerr=moments_c_e_uncert, label='FACEMC-Hybrid')
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)

    plt.show()
    fig1.savefig('./track_flux_comparison.pdf', bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
