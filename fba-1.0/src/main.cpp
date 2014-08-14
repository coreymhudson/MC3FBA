//
//  main.cpp
//  FluxAnalysis
//
//  Created by Andi Dhroso on 10/21/12.
//  Copyright (c) 2012 Andi Dhroso. All rights reserved.
//


#include <unistd.h>
#include <iostream>
#include "Biomass.h"

#define PRINT 1

typedef std::map<std::string, std::string> Parameters;
Parameters loadParameters(const std::string &);
void setParameters(const Parameters &);

void usage();
void displayMenu();


int main(int argc, const char * argv[]) {
    using namespace std;
    
    if (argc < 2) {
        usage();
        return EXIT_FAILURE;
    } else if(strncmp(argv[1], "-h", 2) == 0) {
        displayMenu();
        return EXIT_SUCCESS;
    } else if(strncmp(argv[1], "-", 1) == 0){
        std::cerr << "error: invalid option" << std::endl;
        displayMenu();
    } else {
        Parameters p = loadParameters(argv[1]);
        setParameters(p);
        
        fba::Biomass biomass(fba::MCMC);
        biomass.optimizeBiomass();
    }
    
    return 0;
}

void setParameters(const Parameters &p) {
    void setOption(std::string, std::string);
    Parameters::const_iterator i, end = p.end();
    for (i = p.begin(); i != end; ++i) {
        setOption(i->first, i->second);
    }
}

void setOption(std::string option, std::string value) {
    if (option == "-r") {
        fba::REACTION_FILE = value;
    } else if (option == "-b") {
        fba::BOUNDARY_FILE = value;
    } else if (option == "-f") {
        fba::OPT_CRITERIA = atoi(value.c_str());
    } else if (option == "-e") {
        fba::EPSILON = atoi(value.c_str());
    } else if (option == "-n") {
        fba::NUM_STATES = atoi(value.c_str());
    } else if (option == "-l") {
        fba::LOCAL_SAMPLE_SIZE = atof(value.c_str());
    } else if (option == "-g") {
        fba::GLOBAL_SAMPLE_SIZE = atof(value.c_str());
    } else if (option == "-i") {
        fba::ITERATIONS = atoi(value.c_str());
    } else if (option == "-t") {
        fba::TEMPERATURE = atoi(value.c_str());
    } else if (option == "-o") {
	fba::BIOMASS_THRESHOLD = atof(value.c_str());
    } else if (option == "-s") {
        fba::REACTION_SUBSET_FILE = value;
    } else if(option == "-q") {
        fba::REDUCE_SEARCH_SPACE = atoi(value.c_str());
    }
    else {
        std::cerr << "option \"" << option << " is not recognized ... ignored" << std::endl;
    }
}

Parameters loadParameters(const std::string &filename) {
    Parameters p;
    std::string line;
    std::vector<std::string> tokens;

    std::ifstream in(filename.c_str());
    if (in.is_open()) {
        while (in.good()) {
            std::getline(in, line);
            if (line.empty() || line[0] == '#') continue;
            CTokenizer<CIsSpace>::Tokenize(tokens, line, CIsSpace());
            p[tokens[0]] = tokens[1];
        }
    } else {
        perror("warning: parameter file cannot be opened");
    }
    return p;
}


void usage() {
    std::cerr <<"Usage: ./fba <Parameter_File>" << std::endl;
    //std::cerr <<"Usage: opt_biomass_rxn <matrix file> <boundary file> Biomass_Rxn_# ((constrain_Rxn) (constraint val)) (-skd)  (-sko) (-dkd)" << std::endl;
}

void displayMenu() {
    std::cout << "#########################################################################################" << std::endl;
    std::cout << "# FBA lib description                                                                   #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    std::cout << "# usage ./fba <path to parameter file>                                                  #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    std::cout << "#Paramenters file should contain the following option                                   #" << std::endl;
    std::cout << "# -h show menu - except this option                                                     #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -r reactions input file                                                               #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -b boundary input file                                                                #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
   
    std::cout << "# -f user optimization flag - default 0                                                 #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -e epsilon (print step) - default 0                                                   #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -n number of processes - states default 4                                             #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -i number of iterations - default 1000                                                #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -l local state threshold - default 0.01                                               #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -g global state threshold - default 0.1                                               #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# Temperature is for state switching: higher=> less mixing, lower=>more mixing          #" << std::endl;
    std::cout << "# -t temperature - default 0.01                                                         #" << std::endl;
    
    std::cout << "# -s reactions considered for randomization only - input file                           #" << std::endl;
    std::cout << "#                                                                                       #" << std::endl;
    
    std::cout << "# -q on/off flag for option s - defualt 0                                               #" << std::endl;
    std::cout << "#########################################################################################" << std::endl;
}























