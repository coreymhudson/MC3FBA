//
//  Biomass2.cpp
//  xcode
//
//  Created by Andi Dhroso on 11/17/12.
//  Copyright (c) 2012 Andi Dhroso. All rights reserved.
//

#include "Biomass.h"

namespace fba {
    //lower - upper boundry
    const int    NegInfinity            = -999999;
    const int    PosInfinity            = 999999;


    int NUM_STATES   = 4;                //number of states
    int EPSILON      = 100;              //skip print step
    int ITERATIONS   = 1000;             //total of number of iterations
    int OPT_CRITERIA = 0;                //optimization number provided by user
    double TEMPERATURE = 0.01;           //temperature for state switching (higher means less mixing lower means more mixing)
    double LOCAL_SAMPLE_SIZE = 0.01;     //local state sample size
    double GLOBAL_SAMPLE_SIZE = 0.1;     //global state sample size

    double MIN_GROWTH;
    double BIOMASS_THRESHOLD = 0.75;

    std::string  REACTION_FILE = "";     //reactions input file
    std::string  BOUNDARY_FILE = "";     //boundary input file

    std::string REACTION_SUBSET_FILE = "";  //reaction randomization subset
    ReactionSubset reactionSubset;          //reaction randomization indices
    int REACTION_SUBSET_SIZE = 0;
    bool REDUCE_SEARCH_SPACE = 0;           //randomize entire set of reactions, or just a subset
    
    void loadReactionSubset(std::vector<std::pair<std::string, int> > &reactionNames) {
        std::string line;
        std::vector<std::string> tokens;
        std::vector<std::pair<std::string, Distribution> > subset;
        
        std::ifstream in(REACTION_SUBSET_FILE.c_str());
        if (in.is_open()) {
            while (in.good()) {
                std::getline(in, line);
                if (line.empty()) continue;
                CTokenizer<CIsSpace>::Tokenize(tokens, line, CIsSpace());
                
                //subset file is expected to have the correct format: reaction_name distribution[0 or 1]
                assert(tokens.size() > 1);
                subset.push_back(std::make_pair(tokens[0], atoi(tokens[1].c_str()))) ;
            }
            in.close();
        } else {
            perror("warning: reaction subset file could not be opened...ignored");
            REDUCE_SEARCH_SPACE = 0;    //if file can't be opened, don't use it
            return ;
        }
        
        //convert reaction names to lower case
        std::for_each(reactionNames.begin(), reactionNames.end(), (void (*)(std::pair<std::string, int> &)) convertToLower);
        std::for_each(subset.begin(), subset.end(), (void (*)(std::pair<std::string, Distribution> &)) convertToLower);
        
        //get reaction index
        int index = 0;
        std::vector<std::pair<std::string, Distribution> >::iterator iter, end = subset.end();
        for (iter = subset.begin(); iter != end; ++iter) {
            index = getReactionIndex(iter->first, reactionNames);
            if ( index != -1) reactionSubset.push_back(std::make_pair(index, iter->second));
        }
        if (reactionSubset.empty()) REDUCE_SEARCH_SPACE = 0;
        else { REACTION_SUBSET_SIZE = static_cast<int>(reactionSubset.size()); }
    }

//    void loadReactionSubset(std::vector<std::pair<std::string, int> > &reactionNames) {
//        std::string line;
//        std::vector<std::string> tokens;
//        std::vector<std::string> subset;
//
//        std::ifstream in(REACTION_SUBSET_FILE.c_str());
//        if (in.is_open()) {
//            while (in.good()) {
//                std::getline(in, line);
//                if (line.empty()) continue;
//                CTokenizer<CIsSpace>::Tokenize(tokens, line, CIsSpace());
//                subset.push_back(tokens[0]);
//            }
//            in.close();
//        } else {
//            perror("warning: reaction subset file could not be opened...ignored");
//            REDUCE_SEARCH_SPACE = 0;    //if file can't be opened, don't use it
//            return ;
//        }
//
//        //convert reaction names to lower case
//        std::for_each(reactionNames.begin(), reactionNames.end(), (void (*)(std::pair<std::string, int> &)) convertToLower);
//        std::for_each(subset.begin(), subset.end(), (void (*)(std::string &)) convertToLower);
//        
//        //get reaction index
//        int index = 0;
//        std::vector<std::string>::iterator iter, end = subset.end();
//        for (iter = subset.begin(); iter != end; ++iter) {
//            index = getReactionIndex(*iter, reactionNames);
//            if ( index != -1) reactionSubset.push_back(index);
//        }
//        if (reactionSubset.empty()) REDUCE_SEARCH_SPACE = 0;
//        else { REACTION_SUBSET_SIZE = static_cast<int>(reactionSubset.size()); }
//    }
    
    void convertToLower(std::pair<std::string, int> &str) {
        std::transform(str.first.begin(), str.first.end(), str.first.begin(), (int(*)(int)) tolower);
    }
    
//    void convertToLower(std::pair<std::string,int> &reaction) {
//        std::transform(reaction.first.begin(), reaction.first.end(), reaction.first.begin(), (int(*)(int)) tolower);
//    }
    
    int getReactionIndex(const std::string &str, std::vector<std::pair<std::string, int> > &reactionNames) {
        std::vector<std::pair<std::string, int> >::iterator it, end = reactionNames.end();
        for (it = reactionNames.begin(); it != end; ++it) {
            if (it->first == str)
                return it->second;
        }
        return -1;
    }
}

#pragma mark - FBA State
fba::State::State() : _successes(0), _optFlag(0), _biomass(0), _local(false) {
    _lp = glp_create_prob();
}

fba::State::~State() {
    glp_delete_prob(_lp);
    _lp = NULL;
}

void fba::State::clone(const State* state) {
    glp_copy_prob(_lp, state->getLP(), true);
    _params     = state->_params;
    _successes  = state->_successes;
    _local      = state->_local;
    _optFlag    = state->_optFlag;
    _biomass    = state->_biomass;
}

void fba::State::clone(const State & state) {
    
    glp_copy_prob(this->_lp, state.getLP(), true);
    this->_params     = state.getParams();
    this->_successes  = state._successes;
    this->_local      = state._local;
    this->_optFlag    = state._optFlag;
    this->_biomass    = state._biomass;
}

fba::LPPtr & fba::State::getLP() {
    return _lp;
}

fba::LPPtr fba::State::getLP() const {
    return _lp;
}

glp_smcp & fba::State::getParams() {
    return _params;
}

const glp_smcp & fba::State::getParams() const {
    return _params;
}

#pragma mark - Biomass class
fba::Biomass::Biomass() : m(NULL),optType(MCMC) {
    myrand.init();
    row_indices = NULL;
    col_indices = NULL;
    coefficients = NULL;

//    for (int i = 0; i < NUM_STATES; ++i)
//        states.push_back(new State);
    
}

fba::Biomass::Biomass(const AbstractConstant type) : m(NULL),  optType(type) {
    myrand.init();
    row_indices = NULL;
    col_indices = NULL;
    coefficients = NULL;
    
//    for (int i = 0; i < NUM_STATES; ++i)
//        states.push_back(new State);
}

fba::Biomass::~Biomass() {
    delete [] row_indices;
    delete [] col_indices;
    delete [] coefficients;
    delete m;

    for (int i = 0; i < states.size(); ++i) {
        delete states[i];
        states[i] = NULL;
    }
    states.clear();
}
#pragma mark -

#pragma mark Biomass Setters
//void fba::Biomass::setReducedSearchSpace(bool reducedSpace) {
//    reducedSearchSpace = reducedSpace;
//}

//void fba::Biomass::setReactionFile(const std::string & file) {
//    reactionFile = file;
//}

//void fba::Biomass::setBoundaryFile(const std::string & file) {
//    boundaryFile = file;
//}

//void fba::Biomass::setOptimizationNum(const int & optNum) {
//    optCriteria = optNum;
//}

//void fba::Biomass::setOptimizationType(const AbstractConstant type) {
//    optCriteria = type;
//}

void fba::Biomass::setOutputOption(const AbstractConstant option ) {
    outputType = option;
}

void fba::Biomass::setNumGlobalStates(const int numStates ) {
    NUM_STATES = numStates;
}

#pragma mark  - Biomass Getters 
ReactionMatrixPtr & fba::Biomass::getReactionMatrix() {
    return m;
}

ReactionMatrixPtr fba::Biomass::getReactionMatrix() const {
    return m;
}

//int fba::Biomass::getOptimizationValue() {
//    return optCriteria;
//}

#pragma mark - Biomass Utils
void fba::Biomass::optimizeBiomass( ) {
    bool ld = loadData();
    bool ma = initializeMatrixAttributes();
    if (ld && ma && NUM_STATES) {

        for (int i = 0; i < NUM_STATES; ++i) states.push_back(new State);

        initializeStates();
        if (REDUCE_SEARCH_SPACE) {
            //get reaction names and their index
            std::vector<std::pair<std::string, int> > reactionNames;
            const int size = m->get_num_reactions();
            for (int i = 0; i < size; ++i){
                reactionNames.push_back(std::make_pair(m->get_reaction(i)->react_name, i));
            }
           
            loadReactionSubset(reactionNames);
            reactionNames.clear();

        }
	const double ov = optimum_objective_value();
	if (ov != -1){
		MIN_GROWTH = ov * BIOMASS_THRESHOLD;
	        optimizeMCMC();
	} else {
		perror("warning: solution not feasible");
	}
    }
}

void fba::Biomass::optimizeMCMC() {
    key_t key = 0;
    int qid = queue_init(key);
    std::set<pid_t> threads;
    pid_t pid;
    
    //create worker processes
    for (int i = 0; i < NUM_STATES; ++i) {
        pid = fork();
        if (pid == 0) {
            optimize(i+1, qid);
            _exit(0);
        } else if(pid > 0) {
            threads.insert(pid);
        }
    }
    
    finish(threads);
    terminate_queue(qid);
}

void fba::Biomass::optimize(const int stateId, const int qid) {
    //Message Q structure
    struct msg buff;
    const size_t size = sizeof(buff) - sizeof(long);
    
    //local data
    int     consecutive_fails = 0;
    bool    sol_found;
    //double step = 0.25, progress = .25;
    
    State &t0 = *states[stateId-1], t1;
    t1 = t0;
    int overallAllSuccess = 1;
    
    //localSuccess = number of successes achieved when this state was the local state
    //globalSuccess = number of successes achieved when this state was a global state
    //localFail = number of fails occured when this state was the local state
    //globalFail = number of fails occured when this state was a global state
    //revertToDefaultValues = frequency of state to revert to default values.
    //selectionFrequency = frequence of state being selected.
    int localSuccess = 0, globalSuccess = 0, localFail = 0,
    globalFail = 0, revertToDefaultValues = 0, selectionFrequency = 0, reselection = 0;
    
    do {
        optimizeState(t1);
        //std::cout << "optimizing state: " << stateId << std::endl;
        sol_found = updateState(t0, t1);
        if (t1._local && !sol_found) {
            if (++consecutive_fails > 50) { //trigger switching states only if num if consecutive failures is greater than 50
                size_t id = stateId;
                while (id == stateId) id = (myrand.rand_int32() % NUM_STATES)+1;        //make sure we don't select the same state
                
                fba::sendMessage(qid, stateId, id, overallAllSuccess, fba::BECOME_LOCAL);
                t0._local = false;
                t1._local = false;
                t0._params.tm_lim = 500;
                t1._params.tm_lim = 500;
            }
            ++localFail;
        } else if (t1._local) {
            overallAllSuccess = t0._successes;    //update number of successes
            if (overallAllSuccess % EPSILON == 0) {
                printSol(t0, overallAllSuccess, stateId);
            }
            int id = shouldSwitch(stateId, overallAllSuccess);
            if (id > 0 && id != stateId) { //switch only if id > 0 and not equal to state this state id, otherwise keep going
                //std::cout << "attempting to switch states from state: " << stateId << " to state id: " << id << std::endl;
                fba::sendMessage(qid, stateId, id, overallAllSuccess, fba::BECOME_LOCAL);
                t0._local = false;
                t1._local = false;
                t0._params.tm_lim = 500;
                t1._params.tm_lim = 500;
            } else if(id == stateId) {
                ++reselection;
            }
          /**  const double percentComplete = (double)overallAllSuccess / ITERATIONS;
            //std::cout << "percent complete: " <<  percentComplete << "threshold: " << std::endl;
            if (percentComplete >= progress) {
                std::cout << "state id: " << stateId << " has processed so far: " << progress*100 << "%" << std::endl;
                progress += step;
            }
           */
            ++localSuccess;       //number of successes this states only -  (testing purposes only) don't pass this to other states
            consecutive_fails = 0; //reset num of consecutive failures if optimization was a success
        } else {
            //global state
            if (sol_found) {
                ++globalSuccess;
                consecutive_fails = 0;
            } else if(++consecutive_fails > 1000){
                ++revertToDefaultValues;
				setStateDefaultBounds(t0);
			} else {
                ++globalFail;
            }
            
            size_t rc = msgrcv(qid, &buff, size, stateId, IPC_NOWAIT);
            if (rc != -1) {
                if (buff.info.flag == fba::BECOME_LOCAL) {
                    t0._local = true;
                    t1._local = true;
                    t1._successes = buff.info.value;
                    t0._successes = buff.info.value;
                    overallAllSuccess = t1._successes;
                    t0._params.tm_lim = 1000;
                    t1._params.tm_lim = 1000;
                    
                    ++selectionFrequency;
                } else if(buff.info.flag == fba::STOP) {
                    break;
                }
            }
        }
        
    } while (!t0._local || (t0._local && overallAllSuccess <= ITERATIONS));   //if state is not local, never stop, unless local state signals for STOP
    
    //tell other processes to finish
    if (t0._local) {
        std::cout << "state id: " << stateId << " issusing STOP command with a total number of successes: " << overallAllSuccess << std::endl;
        for (int i = 1; i <= NUM_STATES; ++i) {
            if (i != stateId)   sendMessage(qid, stateId, i, 0, fba::STOP);
        }
    }
    std::cout << "Worker process with id: " << stateId << " finished with overall num success:" << overallAllSuccess <<
    ", local success: " << localSuccess << ", global success: " << globalSuccess << ", local fail: " << localFail << ", global fail: " << globalFail <<
    ", default values: " << revertToDefaultValues << ", frequency selection: " << selectionFrequency <<
    ", re-selection frequency: " << reselection << std::endl;
}

void fba::optimizeState(fba::State & state) {
    fba::LPPtr & lp = state.getLP();
    glp_simplex(lp, &state.getParams());
    
    state._biomass = glp_get_obj_val(lp);
    state._optFlag = glp_get_prim_stat(lp);
}

bool fba::Biomass::updateState(State &t0, State &t1) {
    bool sol_found;
    if (t1._optFlag == GLP_FEAS && t1._biomass > 9575.775) {
        ++t1;
        sol_found = true;
        t0 = t1;
    } else {
        sol_found = false;
        t1 = t0;
    }
    if (t1._local)  updateBoundaries(t0, t1, LOCAL_SAMPLE_SIZE, 1);
    else            updateBoundaries(t0, t1, GLOBAL_SAMPLE_SIZE, 2);
    
    return sol_found;
}

void fba::sendMessage(const int qid, const size_t from, const size_t to, const int value, const fba::Options flag) {
    struct fba::msg buff;
    const size_t size = sizeof(struct fba::msg) - sizeof(long);
    
    buff.mtype = to ;
    buff.info.from = from;
    buff.info.value = value;
    buff.info.flag = flag;
    
    size_t rc = msgsnd(qid, &buff, size, 0);
    assert(rc != -1);
}

//return a state index
int fba::Biomass::shouldSwitch(const int stateId, const int success) {
	double max=0.0;
	int max_index = 0;
	double switch_val;
	
    for(int i = 0; i < NUM_STATES; i++){
        switch_val = 1.0 / ( ((double)i+1.0) * (myrand.rand_open01()+TEMPERATURE) );
        if(switch_val > max){
            max = switch_val;
            max_index = i;
        }
    }
    
    if(max_index==0)
        return stateId;     //does this mean we don't switch
    else
        return (myrand.rand_int32() % NUM_STATES) +1;
    
    //return -1;
}

double fba::Biomass::recipMean(const State & state, int count) {
	double sum = 0.0;
	double prim;
	double j = 0.0;
	for(int i = 0; i < count; i++){
		prim = glp_get_col_prim(state._lp, i+1);
		if(fabs(prim) != PosInfinity){
			sum += fabs(prim);
			j += 1.0;
		}
	}
	return (1.0 / (sum/j));
}

double fba::Biomass::metCouple(const int localSuccesses, const int globalSuccesses) {
	return(static_cast<double>(localSuccesses * LOCAL_SAMPLE_SIZE) / static_cast<double>(globalSuccesses * GLOBAL_SAMPLE_SIZE));
}

void fba::Biomass::findNorm(const State & state, int count, double * mean, double * sd){
	double sum = 0.0;
	double sumsq = 0.0;
	double j = 0.0;
	double var;
	int i;
	double prim;
	for(i = 0; i < count; i++){
		prim = glp_get_col_prim(state._lp, i+1);
		if(fabs(prim) != static_cast<double>(PosInfinity)){
			sum += prim;
			sumsq += prim*prim;
			j += 1.0;
		}
	}
	*mean = sum / j;
	var = (((sum*sum)-sumsq)/(j-1.0))/j;
	*sd = sqrt(var);
}

void fba::Biomass::updateBoundaries(const State &t0, State &t1, double sampleSize, int samplers) {
    double prim;
    int i, index, distribution, numReactions = m->get_num_reactions();
	int bit;
	double mean;
	double sd;
	double exp_dist;
	double norm_dist;
	size_t change = (REDUCE_SEARCH_SPACE) ? myrand.rand_int32() % REACTION_SUBSET_SIZE : myrand.rand_int32() % numReactions;     //random reaction index
	while(change == OPT_CRITERIA || m->get_reaction(change-1)->exchange == true) {                                               //loop until exchange is false or random != optNum
		change = static_cast<int>(myrand.rand_int32() % numReactions);
	}
	if(!(REDUCE_SEARCH_SPACE)){	
	    size_t change2 = change;
		if(samplers != 1){
			change2 = myrand.rand_int32() % numReactions;
			while(change2 == OPT_CRITERIA || m->get_reaction(change2-1)->exchange == true || change2 == change){
				change2 = static_cast<int> (myrand.rand_int32() % numReactions);
			}
		}
	    const LPPtr &lp0 = t0.getLP();
    	LPPtr &lp1 = t1.getLP();
		double lambda = recipMean(t1, numReactions);
    	for(i = 1; i <= numReactions; i++) {
			prim = glp_get_col_prim( lp0, i);
			if (i == change || i == change2) {
				bit = random_bit();
				exp_dist = myrand.exponential(lambda);
				while(exp_dist > PosInfinity) {
					exp_dist = myrand.exponential(lambda);
				}
				if(bit == 0)
        	        glp_set_col_bnds(lp1, i, GLP_DB, -exp_dist, exp_dist);    //lower and upper boundaries on current state
				
        	    else if(bit == 1){
					int neg_bit = random_bit();
					if(neg_bit == 0)
						glp_set_col_bnds(lp1, i, GLP_DB, NegInfinity, -exp_dist);
        	        
        	        else if(neg_bit == 1)
        	            glp_set_col_bnds(lp1, i, GLP_DB, exp_dist, PosInfinity);
				}
			} else if (i != OPT_CRITERIA && m->get_reaction(i-1)->exchange == false) {
        	    if(fabs(prim) > 0.001 && fabs(prim) < PosInfinity && myrand.rand_open01() < sampleSize)
					glp_set_col_bnds(lp1, i, GLP_FX, prim, prim);
        	
        	} else {
        	    switch(m->get_reaction(i-1)->type) {
        	        case BD_FREE:
        	            glp_set_col_bnds(lp1, i, GLP_FR, 0.0, 0.0);
        	            break;
        	        case BD_LOWER:
        	            glp_set_col_bnds(lp1, i, GLP_LO, m->get_reaction(i-1)->lower_bound, 0.0);
        	            break;
        	        case BD_UPPER:
        	            glp_set_col_bnds(lp1, i, GLP_UP, 0.0, m->get_reaction(i-1)->upper_bound);
        	            break;
        	        case BD_BOTH:
            	        glp_set_col_bnds(lp1, i, GLP_DB, m->get_reaction(i-1)->lower_bound, m->get_reaction(i-1)->upper_bound);
                	    break;
            	}
			}
		}
	}
	else{
		const LPPtr &lp0 = t0.getLP();
    	LPPtr &lp1 = t1.getLP();
		findNorm(t1, numReactions, &mean, &sd);
		double lambda = 1.0/mean;
		for(std::vector<int>::size_type i = 0; i != reactionSubset.size(); i++){
			index = reactionSubset.at(i).first;     //first: index, second: distribution
			distribution = reactionSubset.at(i).second;
			prim = glp_get_col_prim(lp0, index);
			if(index == change){
				if(distribution == 0){
					exp_dist = myrand.exponential(lambda);
					if(random_bit() == 0){
						glp_set_col_bnds(lp1, index, GLP_DB, -exp_dist, exp_dist);
					} else{
						if(random_bit() == 0){
							glp_set_col_bnds(lp1, index, GLP_DB, NegInfinity, -exp_dist);
						} else {
							glp_set_col_bnds(lp1, index, GLP_DB, exp_dist, PosInfinity);
						}
					}
				} else if(distribution == 1){
					norm_dist = myrand.normal(mean, sd);
					norm_dist = fabs(norm_dist);
					if(random_bit() == 0){
						glp_set_col_bnds(lp1, index, GLP_DB, -norm_dist, norm_dist);
					} else{
						if(random_bit() == 0){
							glp_set_col_bnds(lp1, index, GLP_DB, NegInfinity, -norm_dist);
						} else{
							glp_set_col_bnds(lp1, index, GLP_DB, norm_dist, PosInfinity);
						}
					}
				}
			} else{
				if(fabs(prim) > 0.001 && fabs(prim) < PosInfinity && myrand.rand_open01() < sampleSize){
					glp_set_col_bnds(lp1, index, GLP_FX, prim, prim);
				} else {
					switch(m->get_reaction(index-1)->type) {
						case BD_FREE:
							glp_set_col_bnds(lp1, index, GLP_FR, 0.0, 0.0);
							break;
						case BD_LOWER:
							glp_set_col_bnds(lp1, index, GLP_LO, m->get_reaction(index-1)->lower_bound, 0.0);
							break;
						case BD_UPPER:
							glp_set_col_bnds(lp1, index, GLP_UP, 0.0, m->get_reaction(index-1)->upper_bound);
							break;
						case BD_BOTH:
							glp_set_col_bnds(lp1, index, GLP_DB, m->get_reaction(index-1)->lower_bound, m->get_reaction(index-1)->upper_bound);
							break;
					}
					
				}
			}
		}
	}
}

bool fba::Biomass::loadData() {
    if (REACTION_FILE.empty() || BOUNDARY_FILE.empty()) return false;
    if (m) {
        delete m;
        m = NULL;
    }

    std::ifstream in1(REACTION_FILE.c_str()), in2(BOUNDARY_FILE.c_str());
    if (!in1.is_open() || !in2.is_open()) {
        perror("error: Reaction or Boundary file not found.");
        return false;
    }
        
    size_t i;
    m = new Reaction_matrix(REACTION_FILE.c_str() );
    std::cout << "Read " << m->get_num_reactions() << " reactions and " << m->get_num_metabolites() << " metabolites" << std::endl;
    
    double max = 0;
    for(i = 0; i < m->get_reaction(OPT_CRITERIA)->num_reactants; i++) {
        if(fabs(m->get_reaction(OPT_CRITERIA)->react_stoich[i]) > max)
            max = fabs(m->get_reaction(OPT_CRITERIA)->react_stoich[i]);
    }
    
    if(max > 10000) {
        for(i = 0; i < m->get_reaction(OPT_CRITERIA)->num_reactants; i++)
            m->get_reaction(OPT_CRITERIA)->react_stoich[i] = m->get_reaction(OPT_CRITERIA)->react_stoich[i] / 1000;
        
        for(i = 0; i < m->get_reaction(OPT_CRITERIA)->num_products; i++)
            m->get_reaction(OPT_CRITERIA)->prod_stoich[i] = m->get_reaction(OPT_CRITERIA)->prod_stoich[i] / 1000;
    }
    
    //load bounds file
    get_bounds(m, BOUNDARY_FILE.c_str());
    
    return true;
}

bool fba::Biomass::initializeMatrixAttributes() {
    if (!m) return false;
    
    size_t i, j;
    non_zero_coeff = 0;
    
    size_t numReactions = m->get_num_reactions(), numReactants, numProducts;
    for(j = 0; j < numReactions; j++) {
        numReactants = m->get_reaction(j)->num_reactants;
        for(i = 0; i < numReactants; i++)
            non_zero_coeff++;
        
        numProducts = m->get_reaction(j)->num_products;
        for(i = 0; i < numProducts; i++)
            non_zero_coeff++;
    }
    
    row_indices = new int [non_zero_coeff+1];
    col_indices = new int [non_zero_coeff+1];
    coefficients = new double [non_zero_coeff+1];
    int cnt = 0;
    numReactions = m->get_num_reactions();
    for(j = 0; j < numReactions; j++) {
        numReactants = m->get_reaction(j)->num_reactants;
        for(i = 0; i < numReactants; i++) {
            cnt++;
            row_indices[cnt] = m->get_reaction(j)->reactants[i]->met_num+1;
            col_indices[cnt] = static_cast<int>(j+1);
            coefficients[cnt] = (double)( m->get_reaction(j)->react_stoich[i]);
        }
        
        numProducts = m->get_reaction(j)->num_products;
        for(i = 0; i < numProducts; i++) {
            cnt++;
            row_indices[cnt] = m->get_reaction(j)->products[i]->met_num+1;
            col_indices[cnt] = static_cast<int>(j+1);
            coefficients[cnt] = (double)( m->get_reaction(j)->prod_stoich[i]);
        }
    }
    
    return true;
}

void fba::Biomass::initializeStates() {
    for (int i = 0; i < NUM_STATES; ++i) 
        initializeState(*states[i]);
    
    (*states[0])._local = true;           //first state is always local
    (*states[0])._params.tm_lim = 1000;   //local state time out before returning with an answer
}

bool fba::Biomass::initializeState(State &state) {
    
    if (!m) return false;
    
    std::string row, col, dummy;
    int i;
    
    LPPtr &lp = state.getLP();
    //create problem object
    glp_set_prob_name(lp, "FBA");
    
    //set optimizaton direction
    glp_set_obj_dir(lp, GLP_MAX);
    
    //Create rows
    glp_add_rows(lp, m->get_num_metabolites());
    
    //Set the constraints for those rows to fixed 0
    for(i = 1; i <= m->get_num_metabolites(); i++) {
        row = "M";
        dummy = toString(i);
        row.append(dummy);
        
        glp_set_row_name(lp, i, row.c_str());
        glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
    }
    
    int optValue = OPT_CRITERIA;
    //Create a column for each reaction
    glp_add_cols(lp, m->get_num_reactions() );
    for(i = 1; i <= m->get_num_reactions(); i++) {
        col = "R";
        dummy = toString(i);
        col.append(dummy);
        
        glp_set_col_name(lp, i, col.c_str());
        if (i == (optValue + 1))
            glp_set_obj_coef(lp, i, 1.0);
        else
            glp_set_obj_coef(lp, i, 0.0);
    }
    setStateDefaultBounds(state);
    
    glp_load_matrix(lp, non_zero_coeff, row_indices, col_indices, coefficients);
    
    glp_smcp &params = state.getParams();
    glp_init_smcp(&params);
    params.msg_lev = GLP_MSG_OFF;
    params.tm_lim = 500;
    //params.presolve = GLP_ON;
    //glp_term_out(DISABLE_OUTPUT);
    
    return true;
}

void fba::Biomass::setStateDefaultBounds(State &state) {
    int i, numReactions = m->get_num_reactions();
    LPPtr &lp = state.getLP();
    for(i = 1; i <= numReactions; i++) {
        switch(m->get_reaction(i-1)->type) {
            case BD_FREE:
                glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
                break;
            case BD_LOWER:
                glp_set_col_bnds(lp, i, GLP_LO, m->get_reaction(i-1)->lower_bound, 0.0);
                break;
            case BD_UPPER:
                glp_set_col_bnds(lp, i, GLP_UP, 0.0, m->get_reaction(i-1)->upper_bound);
                break;
            case BD_BOTH:
                glp_set_col_bnds(lp, i, GLP_DB, m->get_reaction(i-1)->lower_bound, m->get_reaction(i-1)->upper_bound);
                break;
        }
    }
}

void fba::Biomass::printRXNS() {
	std::cout << "GEN\t";
	int i;
	for(i = 1; i<= m->get_num_reactions(); i++){
		std::cout << m->get_reaction(i-1)->react_name << "\t";
	}
	std::cout << std::endl;
}

void fba::Biomass::printSol(const State &state, const int gen, const int thread) {
    if (!state._lp) return ;
    
    int i;
    double val;
    //    printf("%d\t%d\t", gen, thread);
    //printf("gen: %d\tthread: %d\t", gen, thread);
    std::cout << "gen: " << gen << "\tthread: " << thread <<"\t";
    for(i=1; i <= m->get_num_reactions(); i++) {
        m->get_reaction(i-1)->flux = 0.0;
        val = glp_get_col_prim(state.getLP(), i);
        //        printf("%f\t", val);
        //printf("%.2f\t", val);
		if(fabs(val) < 0.001){
			val = 0.0;
		}
		else if((int)val == 999999 || (int)val == -999999){
        	val = 0.0;
        }
		std::cout << val << "\t";
		m->get_reaction(i-1)->flux = val;
    }
    std::cout << std::endl;
}

double fba::Biomass::optimum_objective_value(){
	glp_prob* lp = glp_create_prob();
	glp_copy_prob(lp, states[0]->getLP(), true);
	glp_smcp params;
	glp_init_smcp(&params);
	params.msg_lev = GLP_MSG_OFF;
	params.tm_lim = 1000;
	glp_simplex(lp, &params);
	double objective_value = glp_get_obj_val(lp);
	int sol_status = glp_get_prim_stat(lp);
	glp_delete_prob(lp);
	return(sol_status == GLP_FEAS) ? objective_value : -1;
}

pid_t fba::getTerminatedChildProcess(int childProcessid) {
    int status;	// catch the status of the child
    pid_t processPid = 0;
    do  {
        pid_t w = waitpid(childProcessid, &status, WUNTRACED | WCONTINUED);
        if (w == -1) {
            perror("waitpid");
            break;
        } else {
            processPid = w;
        }
        
        if (WIFEXITED(status)) {
            if (status > 0) {
                std::cerr << "Child process ("<< processPid <<") exited with non-zero status of " << WEXITSTATUS(status) << std::endl;
                continue;
            } else {
                //std::cout << "Child process ("<< processPid <<") exited with status of " << WEXITSTATUS(status) << std::endl;
                continue;
            }
        } else if (WIFSIGNALED(status)) {
            std::cout << "Child process ("<< processPid <<") killed by signal (" << WTERMSIG(status) << ")" << std::endl;
            continue;
        } else if (WIFSTOPPED(status)) {
            std::cout << "Child process ("<< processPid <<") stopped by signal (" << WSTOPSIG(status) << ")" << std::endl;
            continue;
        } else if (WIFCONTINUED(status)){
            std::cout << "Child process ("<< processPid <<") continued" << std::endl;
            continue;
        }
    }
    while (!WIFEXITED(status) && !WIFSIGNALED(status));
    return processPid;
}

void fba::finish(std::set<pid_t> &threads) {
    while (!threads.empty()) {
        const pid_t pid = getTerminatedChildProcess();
        threads.erase(pid);
    }
}

int fba::queue_init(key_t &key) {
    int qid;
    if( (key = ftok( REACTION_FILE.c_str(), 'a' )) == -1 ) {
        perror("ftok: queue_init" );
        exit(1);
    }
    
    if( (qid = msgget( key, 0666 | IPC_CREAT)) == -1 ) {
        perror("msgget: queue_init" );
        exit(1);
    }
    return qid;
}

void fba::terminate_queue(int qid) {
    if ( msgctl(qid, IPC_RMID, NULL) == -1 ) {
        perror("mcgct: terminate_queue. remove msg queue");
        assert(false);
    }
}

#pragma mark - 












