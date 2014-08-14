

#include "stoich_mat.h"

#define NULL_TOL 1e-12

using namespace std;

int Metabolite::operator== (Metabolite test) {
	int retval = 0;
	if(strcmp(test.met_name, met_name) == 0) retval = 1;
	return(retval);
}

Metabolite& Metabolite::operator= (Metabolite &assign_from) {
	strcpy(met_name, assign_from.met_name);
	met_num=assign_from.met_num;
	return(*this);
}

Reaction::Reaction() {
	reactants = products = NULL;
	react_stoich = prod_stoich = NULL;
	num_reactants = num_products = 0;
	reversible = TRUE;
	type = BD_FREE;
	lower_bound = upper_bound = flux = 0;
}

//TODO: Before uncommenting code below
//make sure Reaction_matrix destructor is
//fixed, otherwise will crash due over
//releasing objects.
template <typename T, typename SIZE >
void release2D(T &t, const SIZE & size) {
//    for (SIZE i = 0; i < size; ++i) {
//        delete t[i];
//    }
    delete [] t;
}

template <typename T >
void release1D(T & t) {
    delete [] t;
}

Reaction& Reaction::operator= (Reaction& assign_from)
{
	int i;
	strcpy(react_name, assign_from.react_name);
	react_num = assign_from.react_num;
	
	num_reactants = assign_from.num_reactants;
	num_products = assign_from.num_products;
	
	if (reactants != NULL)
        release2D(reactants, num_reactants);
    
	if (products != NULL)
        release2D(products, num_products);
    
	if (react_stoich != NULL)
        release1D(react_stoich);
    
	if (prod_stoich != NULL) release1D(prod_stoich);
	
	reversible = assign_from.reversible;
    type = assign_from.type;
	upper_bound = assign_from.upper_bound;
	lower_bound = assign_from.lower_bound;
	flux = assign_from.flux;
	exchange = assign_from.exchange;
	
    reactants = NULL;
	if (num_reactants > 0) {
		reactants = new Metabolite*[num_reactants];
		react_stoich = new double[num_reactants];
	}
	
    products = NULL;
	if (num_products > 0) {
		products = new Metabolite*[num_products];
		prod_stoich = new double[num_products];
	}
	
	for(i = 0; i < assign_from.num_reactants; i++) {
        //	reactants[i]=assign_from.reactants[i];
		react_stoich[i] = assign_from.react_stoich[i];
	}
	
	for(i = 0; i < assign_from.num_products; i++) {
        //	products[i]=assign_from.products[i];
		prod_stoich[i] = assign_from.prod_stoich[i];
	}
	return(*this);
}

void Reaction::set_products_reactants(Metabolites & reactant_list, Metabolites & product_list, std::vector<double> & r_stoich, std::vector<double> & p_stoich)
{
	int i;
	num_products = product_list.size();
	num_reactants = reactant_list.size();
	
	if (reactants != NULL)
        release2D(reactants, num_reactants);
	if (products != NULL)
        release2D(products, num_products);
	if (react_stoich != NULL)
        release1D(react_stoich);
	if (prod_stoich != NULL)
        release1D(prod_stoich);
	
	reactants = NULL;
	if (num_reactants > 0) {
		reactants = new Metabolite*[num_reactants];
		react_stoich = new double[num_reactants];
	}
	
    products = NULL;
	if (num_products > 0) {
		products = new Metabolite*[num_products];
		prod_stoich = new double[num_products];
	}
	
    for(i = 0; i < num_reactants; i++) {
		reactants[i] = reactant_list.at(i);
		react_stoich[i] = r_stoich.at(i);
	}
	for(i = 0; i < num_products; i++) {
		products[i] = product_list.at(i);
		prod_stoich[i] = p_stoich.at(i);
	}
}


Reaction::~Reaction() { 
    if (reactants != NULL) {
		release2D(reactants, num_reactants);
		release1D(react_stoich);
	}
    
	if(products != 0) {
        release2D(products, num_products);
        release1D(prod_stoich);
	}
}

Reaction_matrix::Reaction_matrix() {
	num_reactions = num_metabolites = 0;
	omit_reaction = -1;
	metabolites = NULL;
	reactions = NULL;
	fortran_matrix = NULL;
	u_mat = NULL;
	v_mat = NULL;
	sing_vals = NULL;
	null_space_span = NULL; 
	
	fortran_cols = new __CLPK_integer;
	fortran_rows = new __CLPK_integer;
	nrhs = new __CLPK_integer;
	trans = new char;
	lwork = new __CLPK_integer;
	info = new __CLPK_integer;
	
	internal_data = FALSE;
	omit_size = 0;
}


Reaction_matrix::Reaction_matrix(const char* infile)
{
    //Goal rid this class of incessant iterating through arrays to find metabolite values
	//To be done through a c++ map
	//Metabolite_hash<string metabolite_name, int metabolite_num>
	//At the very end this structure will be iterated once to put all the names
	//and numbers into a single Metabolite class
	map<string, int> Metabolite_hash;
	
	int i, num_react, num_prod, reaction_cnt;
	double stoich;
	char dump, line[1600], new_name[NAME_LEN];

	ifstream fin;
	Metabolite a_metabolite;
    MetabolitePtr my_metabolite, new_metabolite;
	Reaction* a_reaction;
    
    std::vector<double> react_stoich, prod_stoich;
	Metabolites reactant_list, product_list;
    
    Reactions reaction_list;
    Metabolites metabolite_list;
	
	fortran_matrix = NULL;
	sing_vals = NULL;
	u_mat = NULL;
	v_mat = NULL;
	null_space_span = NULL;
	
	fortran_cols = new __CLPK_integer;
	fortran_rows = new __CLPK_integer;
	nrhs = new __CLPK_integer;
	trans = new char;
	lwork = new __CLPK_integer;
	info = new __CLPK_integer;
	internal_data = TRUE;
	omit_reaction = -1;
	omit_size = 0;
	fin.open(infile);
	
	num_metabolites = 0;
	reaction_cnt = 0;
	
	fin.getline(line, 1599);
	fin >> new_name >> dump;
	while(strcmp(new_name, "(Exchange") != 0) {
        reactant_list.clear();
		product_list.clear();
        react_stoich.clear();
		prod_stoich.clear();
        
		a_reaction = new Reaction;
		
		num_react = num_prod = 0;
		a_reaction->react_num = reaction_cnt;
		a_reaction->exchange = FALSE;
		strcpy(a_reaction->react_name, new_name);
		reaction_cnt++;
		
		if ((dump == 'i') || (dump == 'I')) a_reaction->reversible = FALSE;
		else a_reaction->reversible = TRUE;
		
        int metNum = 0;
        std::map<std::string, int>::iterator it;
		while (dump != '\n') {
			fin >> stoich >> a_metabolite.met_name;
			
            it = Metabolite_hash.find(a_metabolite.met_name);
			if(it == Metabolite_hash.end() ) {
				Metabolite_hash[a_metabolite.met_name] = num_metabolites;
				
                new_metabolite = new Metabolite;
				a_metabolite.met_num = num_metabolites;
				
                (*new_metabolite) = a_metabolite;
				metabolite_list.push_back(new_metabolite);
                num_metabolites++;
				my_metabolite = new_metabolite;
			}
			else {
				my_metabolite = metabolite_list.at(it->second);
			}
            
            if (stoich < 0) {
				reactant_list.push_back(my_metabolite);
				react_stoich.push_back(stoich);
				num_react++;
			} else {
                product_list.push_back(my_metabolite);
				prod_stoich.push_back(stoich);
				num_prod++;
			}
			
			metNum++;
			fin.get(dump);
		}
        a_reaction->set_products_reactants(reactant_list, product_list, react_stoich, prod_stoich);
		reaction_list.push_back(a_reaction);
		
		fin >> new_name >> dump;
	}
	
	fin.getline(line, 1599);
	
	reactant_list.clear();
	react_stoich.clear();
	num_react = 0;
	
	while(! fin.eof()) {
		fin >> a_metabolite.met_name >> line;
		
		if (! fin.eof()) {
			product_list.clear();
			prod_stoich.clear();
			num_prod=1;
			
			//Get rid of this search and access hash
            std::map<std::string, int>::iterator it = Metabolite_hash.find(a_metabolite.met_name);
			if(it == Metabolite_hash.end()){
				new_metabolite = new Metabolite;
				a_metabolite.met_num = num_metabolites;
				
                (*new_metabolite) = a_metabolite;
				metabolite_list.push_back(new_metabolite);
				cout << "Found new: |" << new_metabolite->met_name << "| " << it->second << endl;
				my_metabolite = new_metabolite;
                num_metabolites++;
			}
			else{
				my_metabolite = metabolite_list.at(it->second);
			}
			
			product_list.push_back(my_metabolite);
			stoich = 1.0;
			prod_stoich.push_back(stoich);
			
			a_reaction = new Reaction;
			
			strcpy(a_reaction->react_name, "Exchng_");
			strcat(a_reaction->react_name, my_metabolite->met_name);
			a_reaction->exchange = TRUE;
			a_reaction->react_num = reaction_cnt;
			reaction_cnt++;
			
			a_reaction->set_products_reactants(reactant_list, product_list, react_stoich, prod_stoich);
			reaction_list.push_back(a_reaction);
		}
	}
	num_reactions = static_cast<int>(reaction_list.size());
	num_metabolites = static_cast<int>(metabolite_list.size());
	
	metabolites = new Metabolite*[num_metabolites];
	for(i = 0; i < num_metabolites; i++)
		metabolites[i] = metabolite_list.at(i);
	
    reactions = new Reaction *[num_reactions];
	for(i = 0; i < num_reactions; i++) {
		reactions[i] = reaction_list.at(i);
		if (reactions[i]->reversible == TRUE)
			reactions[i]->type = BD_FREE;
		else
			reactions[i]->type = BD_LOWER;
		
		reactions[i]->lower_bound = 0.0;
		reactions[i]->upper_bound = 0.0;
		
	}
#if 0
	for(i=0; i<num_metabolites; i++)
		cout<<metabolites[i]->met_name<<"\t"<<metabolites[i]->met_num<<"\t"<<metabolites[i]<<endl;
	
	for(i=0; i<num_reactions; i++) {
		cout<<reactions[i]->num_react<<"\t"<<reactions[i]->num_prod<<"\t"<<reactions[i]->react_name<<"\t"<<reactions[i]<<endl;
		for(j=0; j<reactions[i]->num_react; j++) {
			cout<<"R: "<<reactions[i]->reactants[j]->met_name<<"\t"<<reactions[i]->reactants[j]->met_num<<"\t"<<reactions[i]->reactants[j]<<endl;
		}
		for(j=0; j<reactions[i]->num_prod; j++) {
			cout<<"P: "<<reactions[i]->products[j]->met_name<<"\t"<<reactions[i]->products[j]->met_num<<"\t"<<reactions[i]->products[j]<<endl;
		}
		
	}
	
#endif
}


Reaction_matrix& Reaction_matrix::operator= (Reaction_matrix& assign_from)
{
	int i, j, k, target;
	
	if (fortran_matrix != 0) delete[] fortran_matrix;
	if (u_mat != 0) delete[] u_mat;
	if (v_mat != 0) delete[] v_mat;
	if(sing_vals != 0) delete[] sing_vals;
	
	for(i=0; i < num_metabolites; i++) delete metabolites[i];
	for(i=0; i < num_reactions; i++) delete reactions[i];
	
	delete[] metabolites;
	delete[] reactions;
	
	num_metabolites = assign_from.get_num_metabolites();
	num_reactions = assign_from.get_num_reactions();
	
	metabolites = new Metabolite* [num_metabolites];
	reactions = new Reaction* [num_reactions];
	
	
	for(i = 0; i < num_metabolites; i++) {
		metabolites[i] = new Metabolite();
		(*metabolites[i]) = (*assign_from.get_metabolite(i));
	}
	
	for(i = 0; i < num_reactions; i++) {
		reactions[i] = new Reaction();
		(*reactions[i]) = (*assign_from.get_reaction(i));
	}
	
	for(i = 0; i < num_reactions; i++) {
		for(j = 0; j < reactions[i]->num_reactants; j++) {
			target=assign_from.get_reaction(i)->reactants[j]->met_num;
			k = 0;
			while(metabolites[k]->met_num != target) k++;
			reactions[i]->reactants[j] = metabolites[k];
		}
		for(j = 0; j<reactions[i]->num_products; j++) {
			target = assign_from.get_reaction(i)->products[j]->met_num;
			k = 0;
			while(metabolites[k]->met_num != target) k++;
			reactions[i]->products[j] = metabolites[k];
		}
	}
	return(*this);
}


void Reaction_matrix::set_omit_reaction(int s)
{
	int i, j, new_met_num;
	
	omit_reaction=s;
	omit_size=1;
	
	for(i = 0; i < num_metabolites; i++)
		metabolites[i]->visit = FALSE;
	
	new_met_num = 0;
	
	for(i = 0; i < num_reactions; i++) {
		if (i != omit_reaction) {
			for(j = 0; j < reactions[i]->num_reactants; j++) {
				if (reactions[i]->reactants[j]->visit == FALSE) {
					reactions[i]->reactants[j]->met_num = new_met_num;
					reactions[i]->reactants[j]->visit = TRUE;
					new_met_num++;
				}
			}
			for(j = 0; j < reactions[i]->num_products; j++) {
				if (reactions[i]->products[j]->visit == FALSE) {
					reactions[i]->products[j]->met_num = new_met_num;
					reactions[i]->products[j]->visit = TRUE;
					new_met_num++;
				}
			}
		}
	}
	
	for(j = 0; j<reactions[omit_reaction]->num_reactants; j++) {
		if (reactions[omit_reaction]->reactants[j]->visit == FALSE) {
			reactions[omit_reaction]->reactants[j]->met_num = new_met_num;
			reactions[omit_reaction]->reactants[j]->visit = TRUE;
			new_met_num++;
		}
	}
	for(j = 0; j<reactions[omit_reaction]->num_products; j++) {
		if (reactions[omit_reaction]->products[j]->visit == FALSE) {
			reactions[omit_reaction]->products[j]->met_num = new_met_num;
			reactions[omit_reaction]->products[j]->visit = TRUE;
			new_met_num++;
		}
	}
}


void Reaction_matrix::real_omit_reaction(int s)
{
	int i,j, cnt, new_met_num;
	Reaction **new_reactions;
	
	omit_reaction=s;
	
	cout<<"Reaction name is "<<reactions[s]->react_name<<endl;
	cout<<reactions[s]->reactants[0]->met_name<<": "<<reactions[s]->reactants[0]->met_num<<endl;
	cout<<reactions[s]->products[0]->met_name<<": "<<reactions[s]->products[0]->met_num<<endl;
	
	new_reactions = new Reaction * [num_reactions-1];
	
	cnt = 0;
	for(i = 0; i<num_reactions; i++) {
		if (i != omit_reaction) {
			new_reactions[cnt] = reactions[i];
			cnt++;
		}
	}
	
	delete[] reactions;
	reactions = new_reactions;
    
	omit_reaction = -1;
	num_reactions = num_reactions-1;
	omit_size = 0;
	
	new_met_num = 0;
	
	for(i = 0; i<num_metabolites; i++)
		metabolites[i]->visit = FALSE;
	
    //	for(i=0; i<num_metabolites; i++)
    //		cout<<metabolites[i]->met_num<<": "<<metabolites[i]->met_num<<" V: "<<metabolite[i]
	
	for(i = 0; i < num_reactions; i++) {
		for(j = 0; j < reactions[i]->num_reactants; j++) {
			if (reactions[i]->reactants[j]->visit == FALSE) {
				reactions[i]->reactants[j]->met_num = new_met_num;
				reactions[i]->reactants[j]->visit = TRUE;
				new_met_num++;
			}
		}
		for(j = 0; j<reactions[i]->num_products; j++) {
			if (reactions[i]->products[j]->visit == FALSE) {
				reactions[i]->products[j]->met_num = new_met_num;
				reactions[i]->products[j]->visit = TRUE;
				new_met_num++;
			}
		}
	}
	
	cout << "Now have " << new_met_num << " reactants\n";
}


void Reaction_matrix::create_fortran_stoic_matrix()
{
	int i, j, offset;
	
	*fortran_cols = num_reactions - omit_size;
	*fortran_rows = num_metabolites;
	
	//cout<<"Omitsize: "<<omit_size<<" Omit: "<<omit_reaction<<endl;
	//cout<<"M: "<<*fortran_rows<<"\tN:"<<*fortran_cols<<endl;
	
	fortran_matrix = new double [ (*fortran_cols)*(*fortran_rows)];
	
	u_mat = new double [(*fortran_rows)*(*fortran_rows)];
	v_mat = new double [(*fortran_cols)*(*fortran_cols)];
	sing_vals = new double [(*fortran_rows)];
	
	for(i = 0; i < *fortran_cols; i++) {
		for(j = 0; j < num_metabolites; j++) {
			fortran_matrix[(i*num_metabolites)+j] = 0.0;
		}
	}
	
	offset = 0;
	for(i = 0; i < num_reactions; i++) {
		//cout<<"Reaction is "<<reactions[i]->react_name<<endl;
		if (i != omit_reaction) {
			for(j = 0; j < reactions[i]->num_reactants; j++) {
				//	dummy=reactions[i];
				//	cout<<"Updating "<<reactions[i]->reactants[j]->met_name<<": "<<reactions[i]->reactants[j]->met_num<<endl;
				fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->reactants[j]->met_num]=reactions[i]->react_stoich[j];
			}
			for(j = 0; j < reactions[i]->num_products; j++) {
				fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->products[j]->met_num]=reactions[i]->prod_stoich[j];
			}
		}
		else {
            offset = 1;
        }
	}
}


void Reaction_matrix::compute_ortho_basis()
{
//	int i,j, k, extra_zero_vecs=0, null_cnt, offset;
//	double prod;
//	BOOL has_non_null, has_error;
//	
//	*trans = 'A';
//	
//	*lwork = 10.0*(*fortran_rows)*(*fortran_rows) + 4*(*fortran_rows);
//	
//	iwork = new __CLPK_integer [8*(*fortran_rows)];
//	work = new double[*lwork];
//	
//	
//	if (fortran_matrix == 0)
//		create_fortran_stoic_matrix();
//	
//	//	for(i=0; i<*fortran_cols; i++) {
//	//	for(j=0; j<*fortran_rows; j++)
//	//		cout<<fortran_matrix[(i*(*fortran_rows))+j]<<"\t";
//	//	cout<<"\n";
//	//}
//	
//	
//	dgesdd_(trans, fortran_rows, fortran_cols, fortran_matrix, fortran_rows, sing_vals, u_mat, fortran_rows, v_mat, fortran_cols, work, lwork, iwork, info );
//	delete[] work;
//	
//	/*cout<<"Vmatrix\n";
//	 for(i=0; i<*fortran_cols; i++) {
//	 for(j=0; j<*fortran_cols; j++)
//	 cout<<v_mat[(i*(*fortran_cols))+j]<<"\t";
//	 cout<<endl;
//	 }
//	 
//	 cout<<"Umatrix\n";
//	 for(i=0; i<*fortran_rows; i++) {
//	 for(j=0; j<*fortran_rows; j++)
//	 cout<<u_mat[(i*(*fortran_rows))+j]<<"\t";
//	 cout<<endl;
//	 }*/
//	
//	has_error = FALSE;
//	
//	for(i = 0; i < *fortran_rows; i++) {
//		prod=0;
//		for(j = 0; j < (num_reactions-omit_size); j++)
//			prod += v_mat[i+(*fortran_cols)*(j)]*v_mat[i+(*fortran_cols)*(j)];
//		if(fabs(prod-1.0) > 1e-7) {
//			has_error = TRUE;
//			// cerr<<"Error: vector with sv: "<<sing_vals[i]<<" has length: "<<prod<<endl;
//		}
//	}
//	
//	if (has_error == TRUE) {
//		cout<<"Re-runnning SVD of reaction " << omit_reaction<<endl;
//		if (3*(*fortran_rows)+*fortran_cols > 5*(*fortran_rows)) {*lwork =3*(*fortran_rows)+*fortran_cols;}
//		else {
//            *lwork = 5*(*fortran_rows);
//        }
//		work = new double [*lwork];
//		
//		for(i = 0; i < *fortran_cols; i++) {
//			for(j = 0; j < num_metabolites; j++) {
//				fortran_matrix[(i*num_metabolites)+j] = 0.0;
//			}
//		}
//		
//		offset = 0;
//		for(i = 0; i < num_reactions; i++) {
//			//cout<<"Reaction is "<<reactions[i]->react_name<<endl;
//			if (i != omit_reaction) {
//				for(j = 0; j < reactions[i]->num_reactants; j++) {
//					//	dummy=reactions[i];
//					//	cout<<"Updating "<<reactions[i]->reactants[j]->met_name<<": "<<reactions[i]->reactants[j]->met_num<<endl;
//					fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->reactants[j]->met_num]=reactions[i]->react_stoich[j];
//				}
//				for(j = 0; j < reactions[i]->num_products; j++) {
//					fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->products[j]->met_num]=reactions[i]->prod_stoich[j];
//				}
//			}
//			else {
//                offset=1;
//            }
//		}
//		dgesvd_(trans, trans, fortran_rows, fortran_cols, fortran_matrix, fortran_rows, sing_vals, u_mat, fortran_rows, v_mat, fortran_cols,
//				work, lwork, info );
//		delete[] work;
//	}
//	
//	dim_null_space = 0;
//	has_error = FALSE;
//	//cout<<"Singular vals:\n";
//	for(i = 0; i<*fortran_rows; i++) {
//		if (fabs(sing_vals[i]) <= NULL_TOL) {
//			//		cout<<sing_vals[i]<<endl;
//			dim_null_space++;
//			extra_zero_vecs++;
//			
//			prod = 0;
//			for(j = 0; j<(num_reactions-omit_size); j++)
//				prod += v_mat[i+(*fortran_cols)*(j)]*v_mat[i+(*fortran_cols)*(j)];
//			if(fabs(prod-1.0) > 1e-7) {
//				has_error = TRUE;
//				cerr<<"Error After correction: vector with sv: "<<sing_vals[i]<<" has length: "<<prod<<endl;
//			}
//		}
//		
//	}
//	
//	for(i = num_metabolites; i < (num_reactions-omit_size); i++) {
//		has_non_null=FALSE;
//		for(j = 0; j<(num_reactions-omit_size); j++) {
//			if (v_mat[i+(*fortran_cols)*(j)] > NULL_TOL) has_non_null = TRUE;
//		}
//		
//		if (has_non_null == TRUE) dim_null_space++;
//	}
//	
//	cout<<"The null space of this matrix is of dimension: "<<dim_null_space<<endl;
//	cout<<"There were "<<extra_zero_vecs<<" extra vectors of singular value 0\n";
//	//cout<<"Omit size: "<<omit_size<<endl;
//	
//	null_space_span=new double*[num_reactions];
//	for(i = 0; i<num_reactions; i++)
//		null_space_span[i]=new double[dim_null_space];
//	
//	null_cnt = 0;
//	
//	for(i= 0; i<num_metabolites; i++) {
//		offset = 0;
//		if  (fabs(sing_vals[i]) <=NULL_TOL) {
//			for(j = 0; j<num_reactions; j++) {
//				if (j != omit_reaction)
//					null_space_span[j][null_cnt] = v_mat[i+(*fortran_cols)*(j-offset)];
//				else {
//					offset=1;
//					null_space_span[j][null_cnt]=0;
//				}
//			}
//			null_cnt++;
//		}
//	}
//	
//	//OFFSET FOR OMIT????
//	
//	for(i = num_metabolites; i < (num_reactions-omit_size); i++) {
//		
//		offset = 0;
//		has_non_null = FALSE;
//		for(j = 0; j<(num_reactions-omit_size); j++) {
//			if (v_mat[i+(*fortran_cols)*(j)] > NULL_TOL) has_non_null = TRUE;
//		}
//		if (has_non_null == TRUE) {
//			for(j = 0; j < num_reactions; j++){
//				if (j != omit_reaction)
//					null_space_span[j][null_cnt] = v_mat[i+(*fortran_cols)*(j-offset)];
//				else {
//					offset = 1;
//					null_space_span[j][null_cnt] = 0;
//				}
//			}
//			
//			null_cnt++;
//		}
//	}
//	
//	for (i = 0; i < dim_null_space; i++) {
//		prod = 0;
//		for(k = 0; k < num_reactions; k++)
//			prod += null_space_span[k][i]*null_space_span[k][i];
//		
//		if (fabs(1.0-prod) >1e-7) cerr<<"ERROR: vector "<<i<<"not of unit length: "<<prod<<"\n";
//		
//		for(j = i+1; j < dim_null_space; j++) {
//			prod = 0;
//			for(k = 0; k < num_reactions; k++)
//				prod += null_space_span[k][i]*null_space_span[k][j];
//			
//			if (fabs(prod) >1e-7) cerr << "ERROR: vectors " << i << " and " << j << " are not orthogonal: "<< prod << "\n";
//		}
//	}
//	
//	//for(i=0; i<num_reactions; i++) {
//	//	for(j=0; j<dim_null_space; j++)
//	//		cout<<null_space_span[i][j]<<"\t";
//	//	cout<<"\n";
//	//}
//	
//	//cout<<"Info: "<<*info<<endl;
//	delete[] iwork;
//	work = 0;
}


Reaction_matrix::~Reaction_matrix()
{
	int i;
	if (fortran_matrix !=0)
		delete[] fortran_matrix;
	
	if (u_mat != 0) delete[] u_mat;
	if (v_mat != 0) delete[] v_mat;
	if (sing_vals != 0) delete[] sing_vals;
	if (null_space_span !=0) {
		for(i = 0; i<num_reactions; i++)
			delete[] null_space_span[i];
		delete[] null_space_span;
		
	}
	
	delete fortran_cols;
	delete fortran_rows;
	delete nrhs;
	delete trans;
	delete lwork;
	delete info;
	
    if (internal_data == TRUE) {
		for(i = 0; i < num_reactions; i++)
			delete reactions[i];
	}
    delete[] reactions;
    
    if (internal_data == TRUE) {
		for(i = 0; i < num_metabolites; i++)
			delete metabolites[i];
	}
	delete[] metabolites;
}

Space_compare::Space_compare()
{
	trans = new char;
	fortran_cols = new long int;
	fortran_rows = new long int;
	nrhs = new long int;
	info = new long int;
	lwork = new long int;
	
	fortran_matrix=fortran_vecs=work=0;
    ;
}

double Space_compare::do_comparison (Reaction_matrix *full_matrix, Reaction_matrix *reduced_matrix)
{
//	int i, j, start;
//	long int k, rows, cols;
//	double *soln_vec, *space_vector, alpha, beta, *proj_vec, sqr, retval = 0;
//	BOOL all_zero;
//    
//    start = 0;
//	while((start < full_matrix->dim_null_space) &&
//		  (fabs(full_matrix->null_space_span[reduced_matrix->get_omit_reaction()][start]) < NULL_TOL)) {start++;}
//	
//	if (start < full_matrix->dim_null_space) {
//		rows = reduced_matrix->dim_null_space;
//		cols = 1;
//		k = reduced_matrix->get_num_reactions();
//		fortran_matrix = new double[reduced_matrix->dim_null_space*reduced_matrix->get_num_reactions()];
//		
//		
//		for(i = 0; i < reduced_matrix->get_num_reactions(); i++) {
//			for(j = 0; j < reduced_matrix->dim_null_space; j++) {
//				//fortran_matrix[(i*reduced_matrix->dim_null_space)+j]=reduced_matrix->null_space_span[i][j];
//				fortran_matrix[(j*reduced_matrix->get_num_reactions())+i] = reduced_matrix->null_space_span[i][j];
//				//	cout<<	fortran_matrix[(j*reduced_matrix->get_num_reactions())+i]<<"\t";
//			}
//			//cout<<endl;
//		}
//		
//		space_vector = new double[full_matrix->get_num_reactions()];
//		for(i = 0; i < full_matrix->get_num_reactions(); i++)
//            space_vector[i] = full_matrix->null_space_span[i][start];
//		//for(i=0; i<full_matrix->get_num_reactions(); i++) cout<<space_vector[i]<<endl;
//		soln_vec = new double[reduced_matrix->dim_null_space];
//		
//		//Multply one of the full null-space vectors by the reduced space
//		
//		alpha = 1.0;
//		beta = 0.0;
//		
//		//cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, rows, cols, k, alpha, fortran_matrix, k, space_vector,
//		//			k, beta, soln_vec, rows);
//		cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, rows, cols, k, alpha, fortran_matrix, k, space_vector,
//					k, beta, soln_vec, rows);
//		
//		//cout<<"Combination to project vector to reduced space\n";
//		//for(i=0; i<reduced_matrix->dim_null_space; i++)
//		//	cout<<soln_vec[i]<<endl;
//		
//		proj_vec = new double[reduced_matrix->get_num_reactions()];
//		
//		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, reduced_matrix->get_num_reactions(), 1, reduced_matrix->dim_null_space,
//					alpha, fortran_matrix, reduced_matrix->get_num_reactions(), soln_vec, reduced_matrix->dim_null_space, beta, proj_vec, reduced_matrix->get_num_reactions());
//		
//		//cout<<"Vector projected to reduced space\n";
//		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
//		//	cout<<proj_vec[i]<<endl;
//		
//		sqr = 0;
//		for(i = 0; i < reduced_matrix->get_num_reactions(); i++) {
//			space_vector[i] = full_matrix->null_space_span[i][start]-proj_vec[i];
//			sqr += space_vector[i] * space_vector[i];
//		}
//		
//		//cout<<"Orthologonal vector\n";
//		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
//		//	cout<<space_vector[i]<<endl;
//		
//		sqr = sqrt(sqr);
//		//cout<<"Normalizing from length: "<<sqr<<endl;
//		for(i = 0; i < reduced_matrix->get_num_reactions(); i++)
//			space_vector[i] = space_vector[i]/sqr;
//		
//		cnt_preclude=0;
//		
//		for(i = 0; i < reduced_matrix->get_num_reactions(); i++) {
//			if (space_vector[i] > NULL_TOL) {
//				all_zero = TRUE;
//				j = 0;
//				while((j < reduced_matrix->dim_null_space) && (all_zero == TRUE)) {
//					if (reduced_matrix->null_space_span[i][j] > NULL_TOL)
//						all_zero = FALSE;
//					else j++;
//				}
//				
//				if(all_zero == TRUE) cnt_preclude++;
//			}
//			
//			if (space_vector[i] != 0.0) {
//				sqr = space_vector[i]*space_vector[i];
//				retval += sqr*(log(sqr)/log(2));
//			}
//		}
//		
//		//cout<<"Orthologonal vector\n";
//		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
//		//	cout<<space_vector[i]<<endl;
//		
//		delete[] fortran_matrix;
//		fortran_matrix = 0;
//		delete[] soln_vec;
//		delete[] space_vector;
//		//dgemm_(transa,transb,fortran_rows,fortran_cols,k,alpha,fortran_matrix,k,space_vector,k,beta,soln_vec,fortran_rows);
//		
//        //return retval;
//	}
//    return retval;
    return 0;
}


Space_compare::~Space_compare()
{
	delete trans;
	delete fortran_cols;
	delete fortran_rows;
	delete nrhs;
	delete info;
	delete lwork;
	
	if (fortran_matrix != 0) delete[] fortran_matrix;
	if(fortran_vecs != 0) delete[] fortran_vecs;
	if(work != 0) delete[] work;
}



RXN_RELATE paired_rxn_line_soln(Reaction_matrix * the_matrix, int rxn1, int rxn2)
{
	int i,start;
	double first_ratio, curr_ratio, len;
	RXN_RELATE ret_val=EMPTY;
	
	i = 0;
	
	while ((i < the_matrix->dim_null_space) && ((fabs(the_matrix->null_space_span[rxn1][i]) <= NULL_TOL) || (fabs(the_matrix->null_space_span[rxn2][i]) <= NULL_TOL))) i++;
	
	if (i < the_matrix->dim_null_space) {
		//If we have tranverse the whole null space, then rxn1 and rxn2 are never non-zero together.  Otherwise, we enter this if statement
		ret_val = LINE;
		len = sqrt( (the_matrix->null_space_span[rxn1][i]*the_matrix->null_space_span[rxn1][i]) +
				 (the_matrix->null_space_span[rxn2][i]*the_matrix->null_space_span[rxn2][i]) );
		first_ratio = fabs(the_matrix->null_space_span[rxn1][i]/len);
		
		start = i;
		for(i = start; i < the_matrix->dim_null_space; i++) {
			if ((fabs(the_matrix->null_space_span[rxn1][i]) > NULL_TOL) && (fabs(the_matrix->null_space_span[rxn2][i]) > NULL_TOL)) {
				len=sqrt( (the_matrix->null_space_span[rxn1][i]*the_matrix->null_space_span[rxn1][i]) +
						 (the_matrix->null_space_span[rxn2][i]*the_matrix->null_space_span[rxn2][i]) );
				curr_ratio = fabs(the_matrix->null_space_span[rxn1][i]/len);
				
				if (fabs(curr_ratio - first_ratio) > 1e-5) ret_val = PLANE;
			}
		}
	}
	return(ret_val);
}


void get_bounds (Reaction_matrix *curr_matrix, const char *bound_file)
{
	int i;
	double up, low;
	BOOL low_inf, high_inf;
	char name[200], low_string[100], high_string[100];
	ifstream fin;
	
	fin.open(bound_file);
	
	if (fin.fail()) {
		cerr<<"Error: could not open bound file " << bound_file << endl;
		return;
	}
	
	while(! fin.eof()) {
#if 0
		fin >> name >> low >> up;
#endif
#if 1
		fin >> name >> low_string >> high_string;
		//cout<<name<<" Chck by "<<low_string[1]<<"\t"<<high_string[1]<<endl;
		if ((strlen(low_string) == 1) || ((low_string[1] >= 48) && (low_string[1] <= 57))) {
			low_inf = FALSE;
			low = string_to_float(low_string);
		} else {
			low_inf = TRUE;
			low = -1e10;
		}
        
		if ((strlen(high_string) == 1) || ((high_string[1] >= 48) && (high_string[1] <= 57))) {
			high_inf = FALSE;
			up = string_to_float(high_string);
		}else {
			high_inf = TRUE;
			up = 1e10;
		}
#endif
		//cout<<name<<": "<<low<<"\t"<<up<<endl;
		i = 0;
		
		while((i < curr_matrix->get_num_reactions()) && (strcmp(curr_matrix->get_reaction(i)->react_name, name) != 0)) i++;
		
		if(i >= curr_matrix->get_num_reactions()) cerr << "Error: could not find reaction for bound name " << name<<endl;
		else {
#if 1
			if (low_inf == FALSE) {
				if (high_inf == FALSE) {
					if (curr_matrix->get_reaction(i)->reversible == TRUE)
						curr_matrix->get_reaction(i)->lower_bound = low;
					else
						curr_matrix->get_reaction(i)->lower_bound = fabs(low);
					
					curr_matrix->get_reaction(i)->upper_bound = up;
					curr_matrix->get_reaction(i)->type = BD_BOTH;
				}
				else {
					curr_matrix->get_reaction(i)->lower_bound = low;
					curr_matrix->get_reaction(i)->type = BD_LOWER;
				}
			}
			else {
				if (high_inf == FALSE) {
					curr_matrix->get_reaction(i)->upper_bound = up;
					curr_matrix->get_reaction(i)->type = BD_UPPER;
				}
				else {
					curr_matrix->get_reaction(i)->type = BD_FREE;
				}
				
			}
#endif
#if 0
			if ((low == -999999) && (up == 0)) {
				cout << "UPPER bound " << name << "\t" << low <<"\t" << up << endl;
				curr_matrix->get_reaction(i)->lower_bound = low;
				curr_matrix->get_reaction(i)->upper_bound = up;
				curr_matrix->get_reaction(i)->type = BD_UPPER;
			}
			else {
				if ((low == 0) && (up == 999999)) {
					cout << "LOWER " << name << "\t" << low << "\t" << up << endl;
					curr_matrix->get_reaction(i)->lower_bound = low;
					curr_matrix->get_reaction(i)->upper_bound = up;
					curr_matrix->get_reaction(i)->type = BD_LOWER;
					
				}
				else {
					if ((low == -999999) && (up == 999999)) {
                        cout << "FREE " << name << "\t" << low << "\t" << up << endl;
						curr_matrix->get_reaction(i)->lower_bound = low;
						curr_matrix->get_reaction(i)->upper_bound = up;
						curr_matrix->get_reaction(i)->type = BD_FREE;
					}
					else {
                        cout << "STRICT " << name << "\t" << low << "\t" << up << endl;
						if (curr_matrix->get_reaction(i)->reversible = TRUE) 
							curr_matrix->get_reaction(i)->lower_bound=low;
						else 
							curr_matrix->get_reaction(i)->lower_bound = fabs(low);
						
						curr_matrix->get_reaction(i)->upper_bound = up;
						curr_matrix->get_reaction(i)->type = BD_BOTH;
					}
				}
			}
#endif
		}
	}
	fin.close();
}

double string_to_float(const char *instring)
{
    int i=0;
    BOOL neg=FALSE;
    double value, tenths=0.1;
    
    while(instring[i]<48 || instring[i]>57)
        i++;
    
    if (instring[i-1]=='-')
        neg=TRUE;
    value=instring[i]-48;
    i++;
    
    while (instring[i]>=48 && instring[i]<=57)
    {
        value=value*10+(instring[i]-48);
        i++;
    }
    if (instring[i]=='.')
    {
        i++;
        while (instring[i]>=48 && instring[i]<=57)
        {
            value=value+((instring[i]-48)*tenths);
            tenths*=0.1;
            i++;
        }
    }
    
    if (neg==TRUE)
        value=-1.0*value;
    
    return(value);
    
}




