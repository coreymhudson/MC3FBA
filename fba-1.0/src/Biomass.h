 //
//  Biomass2.h
//  xcode
//
//  Created by Andi Dhroso on 11/17/12.
//  Copyright (c) 2012 Andi Dhroso. All rights reserved.
//

#ifndef __xcode__Biomass2__
#define __xcode__Biomass2__

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <cstring>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <stdlib.h>
#include <assert.h>


#include <glpk.h>

#include "stoich_mat.h"
#include "random.h"
#include "Tokenizer.h"

#pragma mark - 

namespace fba {
    //lower - upper boundry
    const extern int      NegInfinity;
    const extern int      PosInfinity;

    extern bool     REDUCE_SEARCH_SPACE;       //randomize entire set of reactions, or just a subset
    extern int      OPT_CRITERIA;               //optimization number provided by user
    extern int      NUM_STATES;                 //number of states
    extern int      EPSILON;                    //skip print step
    extern int      ITERATIONS;                 //total of number of iterations
    extern double   TEMPERATURE;                //temperature for state switching (higher means less mixing lower means more mixing)
        
    extern double   LOCAL_SAMPLE_SIZE;          //local state sample size
    extern double   GLOBAL_SAMPLE_SIZE;         //global state sample size
    
    extern double   MIN_GROWTH; 		//max growth, i.e., biomass
    extern double   BIOMASS_THRESHOLD;     	//biomass threshold

    extern std::string  REACTION_FILE;          //reactions input file
    extern std::string  BOUNDARY_FILE;          //boundary input file

    extern  std::string REACTION_SUBSET_FILE;                   //reaction randomization subset
    typedef int Index, Distribution;
    typedef std::vector<std::pair<Index, Distribution> > ReactionSubset;   //reaction randomization indices
    extern ReactionSubset reactionSubset;                       //reaction randomization indices
    
    void loadReactionSubset(std::vector<std::pair<std::string &,int> > &);
    void convertToLower(std::string &);
    void convertToLower(std::pair<std::string,int> &);
    int getReactionIndex(const std::string &, std::vector<std::pair<std::string, int> > &);
    
    void setReducedSearchSpace(bool = false);
    void setReactionFile(const std::string & );
    void setBoundaryFile(const std::string & );
    void setOptimizationNum(const int);
    
    int getOptimizationValue();
    
    
    
    enum AbstractConstant {
        //Probability distribution
        EXP = 1, NORM, UNIFORM,
        
        //Enable/Disable linear programmer ouput options
        DISABLE_OUTPUT, ENABLE_OUTPUT,
        
        //Ouput solution options
        CONSOLE, FILE,
        
        //Optimization method types
        MCMC, BINARY
    };

    
typedef glp_prob LP;
typedef LP* LPPtr;

//typedef struct state State;
class State {
    public:
    State();
    ~State();
    
    State& operator=(const State &rhs) {
        glp_copy_prob(_lp, rhs.getLP(), true);
        _params     = rhs._params;
        _successes  = rhs._successes;
        _local      = rhs._local;
        _optFlag    = rhs._optFlag;
        _biomass    = rhs._biomass;
        
        return *this;
    }
    
    void operator++() {
        ++_successes;
    }
    
    void clone(const State* state);
    void clone(const State & state);
    LPPtr & getLP();
    LPPtr getLP() const;
    
    glp_smcp & getParams();
    const glp_smcp & getParams() const;
    
    LPPtr       _lp;
    glp_smcp    _params;
    int         _successes;
    int         _optFlag;
    double      _biomass;
    bool        _local;
};
typedef State* StatePtr;
typedef std::vector<StatePtr> States;

    
class Biomass {
public:
    explicit Biomass();
    explicit Biomass(const AbstractConstant type);
    virtual ~Biomass();
    
    #pragma mark Setters
    void setOptimizationType(const AbstractConstant type);
    void setOutputOption(const AbstractConstant option = CONSOLE );
    void setNumGlobalStates(const int globalStates  = 3);
    
    #pragma mark Getters
    ReactionMatrixPtr & getReactionMatrix();
    ReactionMatrixPtr getReactionMatrix() const;

    #pragma mark - Utils
    void optimizeBiomass();
    void optimizeMCMC();
    void optimize(const int stateIndex, const int qid);
    bool updateState(State &t0, State &t1);
    int  shouldSwitch(const int stateId, const int sucsess);
    
    size_t switchState(const StatePtr t0, const size_t from, const int qid);
    
    void findNorm(const State &state, int count, double *mean, double *sd);
    double recipMean(const State &state, int count);
    double metCouple(const int localSuccesses, const int globalSuccesses);
    void updateBoundaries(const State &t0, State &t1, double sampleSize, int samplers);
    
    bool loadData();
//    bool loadData(const std::string &reactions, const std::string &boundaries);
    bool initializeMatrixAttributes();
    void initializeStates();
    bool initializeState(State &state);
    void setStateDefaultBounds(State &state);

    void printRXNS();
    void printSol(const State &state, const int gen, const int i);
    double optimum_objective_value();
    
    private:
    int                 non_zero_coeff;
    int*                row_indices;
    int*                col_indices;
    double*             coefficients;

    ReactionMatrixPtr   m;
    AbstractConstant    optType;                //mcmc, binary ...etc
    AbstractConstant    outputType;             //console, file
    States              states;                 //local, global
    KW_RNG::RNG         myrand;
    
};//end Biomass

template <typename T>
std::string toString(T t) {
    std::stringstream converter;
    converter << t;
    
    return converter.str();
}
    
enum Options {
    BECOME_LOCAL = 1,
    SUCCESS = 2,
    STOP
};
    
struct msg{
    long mtype;               //To:
    struct buff {
        size_t from;          //From:
        int value;
        Options flag;
    }info;
};

void sendMessage(const int qid, const size_t from, const size_t to, const int value, const Options flag);

void optimizeState(fba::State & state);
int queue_init(key_t &key);
void terminate_queue(int qid);
void finish(std::set<pid_t> &threads);
pid_t getTerminatedChildProcess(int childProcessid = -1);
    
}//end namespace fba

#endif /* defined(__xcode__Biomass2__) */
#pragma mark - 





