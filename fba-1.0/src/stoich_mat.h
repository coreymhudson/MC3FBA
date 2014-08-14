


//extern "C" {
//#include "cblas.h"
//#include "clapack.h"
//}

#ifndef ___STOICH_MAT_H___
#define ___STOICH_MAT_H___
#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <vector>

//#include <vecLib/cblas.h>
//#include <vecLib/clapack.h>

#include <math.h>


//#include "gen_dna_funcs.h"

#define NAME_LEN 300
enum BOUND_TYPE {BD_FREE, BD_UPPER, BD_LOWER, BD_BOTH};

enum RXN_RELATE {PLANE, LINE, EMPTY};

enum BOOL {FALSE, TRUE};
typedef int __CLPK_integer;

class Metabolite {
public:
	int met_num;
	char met_name[NAME_LEN];
	BOOL visit;
	int operator==(Metabolite test);
	Metabolite& operator= (Metabolite &assign_from);
};
typedef Metabolite* MetabolitePtr;
typedef std::vector<MetabolitePtr> Metabolites;

class Reaction {
public:
	size_t react_num, num_reactants, num_products;
	double *react_stoich, *prod_stoich, upper_bound, lower_bound, flux, min, max;
	char react_name[NAME_LEN];
	Metabolite **reactants, **products; 
	BOUND_TYPE type;
	BOOL reversible, exchange;
	
	Reaction();
	Reaction& operator= (Reaction& assign_from);
	void set_products_reactants(Metabolites & reactant_list, Metabolites & product_list, std::vector<double> &r_stoich, std::vector<double> & p_stoich);
	~Reaction();
};
typedef Reaction* ReactionPtr;
typedef std::vector<ReactionPtr> Reactions;

class Reaction_matrix {
public:
	int dim_null_space;
	__CLPK_integer *nrhs, *fortran_cols, *fortran_rows, *info, *lwork, *pivots, *iwork;
	char *trans;
	double *fortran_matrix, *sing_vals, *u_mat, *v_mat,  *work,
	**null_space_span;
	
	Reaction_matrix();
    
	Reaction_matrix(const char* infile);

    ~Reaction_matrix();
	
    Reaction_matrix& operator= (Reaction_matrix& assign_from);
    
	int get_num_reactions() {
        return(num_reactions);
    }
    
	int get_num_metabolites() {
        return(num_metabolites);
    }
	
    void create_fortran_stoic_matrix();
	
    void compute_ortho_basis();
	
    Metabolite* get_metabolite(size_t i)  {
        return(metabolites[i]);
    }
    
	ReactionPtr & get_reaction(size_t i) {
        return (reactions[i]);
    }
	
    //get a reference of reactions
    Reaction*** getReactions( ) {
        return &reactions;
    }
    
    size_t get_omit_reaction() {
        return(omit_reaction);
    }
	
    void set_omit_reaction(int s);

	void real_omit_reaction(int s);
	
protected:
	int num_reactions, num_metabolites, omit_reaction, omit_size;
	BOOL internal_data;
	Metabolite **metabolites;
	Reaction **reactions;       
};
typedef Reaction_matrix* ReactionMatrixPtr;


class Space_compare {
public:
	int cnt_preclude;
	long int *fortran_cols, *fortran_rows, *nrhs, *info, *lwork;
	char *trans;
	double *fortran_matrix, *fortran_vecs, *work;
	
	Space_compare();
	double do_comparison (Reaction_matrix *full_matrix, Reaction_matrix *reduced_matrix);
	~Space_compare();
};

RXN_RELATE paired_rxn_line_soln(Reaction_matrix * the_matrix, int rxn1, int rnx2);

void get_bounds (Reaction_matrix *curr_matrix, const char *bound_file);
double string_to_float(const char *instring);
#endif


