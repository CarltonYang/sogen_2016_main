/*
 Simulation for zebrafish segmentation
 Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "sim.hpp" // Function declarations
#include "array2D.hpp"
#include "debug.hpp"
#include "feats.hpp"
#include "init.hpp"
#include "io.hpp"
#include "structs.hpp"
#include <vector>
using namespace std;

extern terminal* term; // Declared in init.cpp

/* simulate_all_params simulates every parameter set after initialization is finished
	parameters:
 ip: the program's input parameters
 rs: the current simulation's rates
 sd: the current simulation's data
 sets: the array of parameter sets
 mds: the array of all mutant data
 file_passed: a pointer to the output file stream of the passed file
 file_scores: a pointer to the output file stream of the scores file
 dirnames_cons: the array of mutant directory paths
 file_features: a pointer to the output file stream of the features file
 file_conditions: a pointer to the output file stream of the conditions file
	returns: nothing
	notes:
	todo:
 TODO consolidate ofstream parameters.
 */
void simulate_all_params (input_params& ip, rates_static& rs, rates_dynamic& rs_d, sim_data& sd, double** sets, mutant_data mds[], ofstream* file_passed, ofstream* file_scores, char** dirnames_cons, ofstream* file_features, ofstream* file_conditions) {
    // Initialize score data
    int sets_passed = 0;
    double score[ip.num_sets];
    
    // Initialize the concentration levels structs
    int max_cl_size = MAX(sd.steps_til_growth, sd.max_delay_size + sd.steps_total - sd.steps_til_growth) / sd.big_gran + 1;

	con_levels cls[ip.num_active_mutants];
	con_levels baby_cls[ip.num_active_mutants];
	for (int i = 0; i < ip.num_active_mutants; i++) {
    	cls[i]= con_levels(NUM_CON_STORE, max_cl_size, sd.cells_total, sd.active_start); // Concentration levels for analysis and storage
		baby_cls[i]= con_levels(NUM_CON_LEVELS, sd.max_delay_size, sd.cells_total, sd.active_start); // Concentration levels for simulating (time in this cl is treated cyclically)	
	}
    
    
    // Simulate every parameter set
    for (int i = 0; i < ip.num_sets; i++) {
        memcpy(rs_d.rates_base, sets[i], sizeof(double) * NUM_RATES); // Copy the set's rates to the current simulation's rates
        score[i] = simulate_param_set(i, ip, sd, rs, rs_d, cls, baby_cls, mds, file_passed, file_scores, dirnames_cons, file_features, file_conditions);
        sets_passed += determine_set_passed(sd, i, score[i]); // Calculate the maximum score and whether the set passed
    }
    
    // Pipe the scores if piping specified by the user
    if (ip.piping) {
        write_pipe(score, ip, sd);
    }
    
    cout << endl << term->blue << "Done: " << term->reset << sets_passed << "/" << ip.num_sets << " parameter sets passed all conditions" << endl;
}

/* simulate_param_set simulates the given parameter set with every specified mutant
	parameters:
 sd: the current simulation's data
 score: the set's received score
	returns: whether or not the set passed
	notes:
	todo:
 */
bool determine_set_passed (sim_data& sd, int set_num, double score) {
    cout << term->blue << "Done: " << term->reset << "set " << set_num << " scored ";
    double max_score = sd.no_growth ? sd.max_scores[SEC_POST] : sd.max_score_all;
    bool passed = score == max_score;
    if (passed) {
        cout << term->blue;
    } else {
        cout << term->red;
    }
    cout << score << " / " << max_score << term->reset << endl;
    return passed;
}

/* simulate_param_set simulates the given parameter set with every specified mutant
	parameters:
 set_num: the index of the parameter set to simulate
 ip: the program's input parameters
 sd: the current simulation's data
 rs: the current simulation's rates
 cl: the concentration levels used for analysis and storage
 baby_cl: the concentration levels used for simulating
 mds: the array of all mutant data
 file_passed: a pointer to the output file stream of the passed file
 file_scores: a pointer to the output file stream of the scores file
 dirnames_cons: the array of mutant directory paths
 file_features: a pointer to the output file stream of the features file
 file_conditions: a pointer to the output file stream of the conditions file
	returns: the cumulative score of every mutant
	notes:
	todo:
 */
double simulate_param_set (int set_num, input_params& ip, sim_data& sd, rates_static& rs, rates_dynamic& rs_d, con_levels cls[], con_levels baby_cls[], mutant_data mds[], ofstream* file_passed, ofstream* file_scores, char** dirnames_cons, ofstream* file_features, ofstream* file_conditions) {
    // Prepare for the simulations
    cout << term->blue << "Simulating set " << term->reset << set_num << " . . ." << endl;
    int num_passed = 0;
    double scores[NUM_SECTIONS * NUM_MUTANTS] = {0};
    if (!ip.reset_seed) { // Reset the seed for each set if specified by the user
        init_seeds(ip, set_num, set_num > 0, true);
    }
    for (int i = 0; i < ip.num_active_mutants; i++) {
    	cls[i].reset(); // Reset the concentration levels for each set
	}
    (*mds).feat.reset();
    
    // Simulate every mutant in the posterior before moving on to the anterior
    int end_section = SEC_ANT * !(sd.no_growth);
    for (int i = SEC_POST; i <= end_section; i++) {
        sd.section = i;
        num_passed += simulate_section(set_num, ip, sd, rs, rs_d, cls, baby_cls, mds, dirnames_cons, scores);
    }
    
    // Calculate the total score
    double total_score = 0;
    for (int i = 0; i < NUM_SECTIONS * ip.num_active_mutants; i++) {
        total_score += scores[i];
    }
    
    // Print the mutant's results (if not short circuiting)
    //if (total_score == sd.max_scores[SEC_POST] + sd.max_scores[SEC_WAVE] + sd.max_scores[SEC_ANT]) {
    if (total_score == sd.max_scores[SEC_POST]  + sd.max_scores[SEC_ANT]) {
        print_passed(ip, file_passed, rs_d);
    }
    print_osc_features(ip, file_features, mds, set_num, num_passed);
    print_conditions(ip, file_conditions, mds, set_num);
    print_scores(ip, file_scores, set_num, scores, total_score);
    
    return total_score;
}

/* simulate_section simulates the given section with every specified mutant
	parameters:
 set_num: the index of the parameter set to simulate
 ip: the program's input parameters
 sd: the current simulation's data
 rs: the current simulation's rates
 cl: the concentration levels used for analysis and storage
 baby_cl: the concentration levels used for simulating
 mds: the array of all mutant data
 dirnames_cons: the array of mutant directory paths
 scores: the array of scores to populate with each mutant's results
	returns: the number of mutants that passed
	notes:
	todo:
 */
int simulate_section (int set_num, input_params& ip, sim_data& sd, rates_static& rs,rates_dynamic& rs_d, con_levels cls[], con_levels baby_cls[], mutant_data mds[], char** dirnames_cons, double scores[]) {
    // Prepare for this section's simulations
    int num_passed = 0;
	//changed array to a 2d array to store all the data
    double temp_rates[ip.num_active_mutants][2]; // Array of knockout rates so knockouts can be quickly applied and reverted
    determine_start_end(sd);
    reset_mutant_scores(ip, mds);

    rates_dynamic rs_ds[ip.num_active_mutants];
    for (int i = 0; i < ip.num_active_mutants; i++) {
        rs_ds[i].cells=sd.cells_total;
        rs_ds[i].rates_active.initialize(NUM_RATES,sd.cells_total);
		rs_ds[i].rates_cell.initialize(NUM_RATES,sd.cells_total);
		memcpy(rs_ds[i].rates_base,rs_d.rates_base, sizeof(double) * NUM_RATES);
    }

    // Simulate each mutant
    for (int i = 0; i < ip.num_active_mutants; i++) {
        mutant_sim_message(mds[i], i);
        store_original_rates(rs_ds[i], mds[i], temp_rates[i]); // will be used to revert original rates after current mutant
        knockout (rs_ds[i], mds[i], 0);
	}
		//cout<<"in simulate_section"<<endl;
        std::vector<double> current_score(ip.num_active_mutants);
		current_score = simulate_mutant(set_num, ip, sd, rs, rs_ds, cls, baby_cls, mds, mds[MUTANT_WILDTYPE].feat, dirnames_cons, temp_rates);


	 for (int i = 0; i < ip.num_active_mutants; i++) {
        scores[sd.section * ip.num_active_mutants + i] = current_score.at(i);
        baby_cls[i].reset();

		//maybe no longer needed?
        revert_knockout(rs_ds[i], mds[i], temp_rates[i]); // this should still happen at the end
        
        if (current_score[i] == mds[i].max_cond_scores[sd.section]) { // If the mutant passed, increment the passed counter
            ++num_passed;
        } else if (ip.short_circuit) { // Exit both loops if the mutant failed and short circuiting is active
            return num_passed;
        }
    }
    
    return num_passed;
}

/* determine_start_end determines the start and end points for the current simulation based on the current section
	parameters:
 sd: the current simulation's data
	returns: nothing
	notes:
	todo:
 */
void determine_start_end (sim_data& sd) {
    if (sd.section == SEC_POST) {
        sd.time_start = 1;
        sd.time_end = MIN(sd.steps_til_growth + 1, sd.steps_total);
    } else {
        sd.time_start = sd.max_delay_size;
        sd.time_end = sd.max_delay_size + sd.steps_total - sd.steps_til_growth;
        //cout<<"MAX DELAY: "<< sd.max_delay_size<<" "<< sd.steps_til_growth<<endl;
    }
}

/* reset_mutant_scores resets each mutant's condition scores
	parameters:
 ip: the program's input parameters
 mds: the array of all mutant data
	returns: nothing
	notes:
	todo:
 */
void reset_mutant_scores (input_params& ip, mutant_data mds[]) {
    for (int i = 0; i < ip.num_active_mutants; i++) {
        for (int j = 0; j < NUM_SECTIONS; j++) {
            for (int k = 0; k < 1 + MAX_CONDS_ANY; k++) {
                mds[i].conds_passed[j][k] = 0;
            }
            mds[i].secs_passed[j] = false;
            
        }
    }
    //mds[MUTANT_WILDTYPE].conds_passed[SEC_WAVE][0]=1;
    //mds[MUTANT_WILDTYPE].conds_passed[SEC_WAVE][1]=1;
    //mds[MUTANT_WILDTYPE].conds_passed[SEC_WAVE][2]=1;
    //mds[MUTANT_WILDTYPE].conds_passed[SEC_WAVE][3]=1;
    mds[MUTANT_WILDTYPE].conds_passed[SEC_ANT][6]=1;
    mds[MUTANT_WILDTYPE].conds_passed[SEC_ANT][7]=1;
}

/* mutant_sim_message prints a message indicating that the given mutant is being simulated in the given section
	parameters:
 md: the mutant being simulated
 section: the section being simulated
	returns: nothing
	notes:
	todo:
 */
inline void mutant_sim_message (mutant_data& md, int section) {
    term->verbose() << term->blue << "  Simulating " << term->reset << md.print_name << " ";
    if (section == SEC_POST) {
        term->verbose() << "before";
    } else {
        term->verbose() << "with";
    }
    term->verbose() << " growth . . ." << endl;
}

/* store_original_rates stores the rates which might be knocked out later in the simulation
	parameters:
 rs: the current simulation's rates
 md: the mutant being simulated
 orig_rates: the array to store the original rates
	returns: nothing
	notes:
	todo:
 */
inline void store_original_rates (rates_dynamic& rs, mutant_data& md, double orig_rates[]) {
    for (int i = 0; i < md.num_knockouts; i++) {
        orig_rates[i] = rs.rates_base[md.knockouts[i]];
    }
}

/* knockout knocks out the given mutant's rates and stores their originals in the given array
	parameters:
 rs: the current simulation's rates
 md: the mutant being simulated
 orig_rates: the array to store the original rates
	returns: nothing
	notes:
	todo:
 
	151221: For DAPT Mutant, knockout after induction
 */
inline void knockout (rates_dynamic& rs, mutant_data& md, bool induction) {
    for (int i = 0; i < md.num_knockouts; i++) {
        if (!(md.index==MUTANT_DAPT && induction == 0)){
            rs.rates_base[md.knockouts[i]] = 0;
        }
    }
}

/* revert_knockout reverts a knockout for the given mutant using the given original rates
	parameters:
 rs: the current simulation's rates
 md: the mutant being simulated
 orig_rates: the array of original rates
	returns: nothing
	notes:
	todo:
 */
inline void revert_knockout (rates_dynamic& rs, mutant_data& md, double orig_rates[]) {
    for (int i = 0; i < md.num_knockouts; i++) {
        rs.rates_base[md.knockouts[i]] = orig_rates[i];
    }
}

/* simulate_mutant simulates the given parameter set with the given mutant
	parameters:
 set_num: the index of the parameter set to simulate
 ip: the program's input parameters
 sd: the current simulation's data
 rs: the current simulation's rates
 cl: the concentration levels used for analysis and storage
 baby_cl: the concentration levels used for simulating
 md: the mutant to simulate
 wtfeat: the oscillation features the wild type produced (or will if the mutant is the wild type)
 dirname_cons: the directory path of the mutant
	returns: the score of the mutant
	notes:
	todo:
 TODO Break up this enormous function.
 */
vector<double> simulate_mutant (int set_num, input_params& ip, sim_data& sd, rates_static& rs, rates_dynamic rs_ds[], con_levels cls[], con_levels baby_cls[], mutant_data mds[], features& wtfeat, char** dirname_cons, double temp_rates[][2]) {

	//initiate the return array
	std::vector<double> score_array(ip.num_active_mutants);

	for (int i = 0; i < ip.num_active_mutants; i++) {
    	reset_seed(ip, sd); // Reset the seed for each mutant
    	baby_cls[i].reset(); // Reset the concentrations levels used for simulating
    	perturb_rates_all(rs,rs_ds[i]); // Perturb the rates of all starting cells
	    //cout<<"in simulate_mutant"<<endl;
	    // Initialize active record data and neighbor calculations
	    sd.initialize_active_data();
	    if (sd.height > 1) {
	        calc_neighbors_2d(sd);
	    }
	    cls[i].active_start_record[0] = sd.active_start;
	    baby_cls[i].active_start_record[0] = sd.active_start;
	    
	    // Copy the posterior results if this is an anterior simulation
	    if (sd.section == SEC_ANT) {
	        copy_mutant_to_cl(sd, baby_cls[i], mds[i]);
	    }
	}    
    // Simulate the mutant
    bool passed = model(ip, sd, rs, rs_ds, cls, baby_cls, mds, temp_rates);
    

	for (int i = 0; i < ip.num_active_mutants; i++) {
    // Analyze the simulation's oscillation features
    	term->verbose() << term->blue << "    Analyzing " << term->reset << "oscillation features . . . ";
    	double score = 0;
    	if (sd.section == SEC_POST) { // Posterior analysis
        	osc_features_post(sd, ip, cls[i], mds[i].feat, wtfeat, mds[i], dirname_cons[i], sd.time_start / sd.big_gran, sd.time_end / sd.big_gran, set_num);
        	term->verbose() << term->blue << "Done" << endl;
    	} else { // Anterior analysis
        	if (ip.ant_features) {
        	    term->verbose() << endl;
        	}
        	osc_features_ant(sd, ip, wtfeat, dirname_cons[i], cls[i], mds[i], 0, sd.height, 0, 5, set_num);
        	if (ip.ant_features) {
        	    term->verbose() << term->blue << "    Done " << term->reset << "analyzing oscillation features" << endl;
        	} else {
        	    term->verbose() << term->blue << "Done" << endl;
        	}
    	}
    
   	 // Copy and print the appropriate data
    	print_concentrations(ip, sd, cls[i], mds[i], dirname_cons[i], set_num);
    	if (sd.section == SEC_ANT) { // Print concentrations of columns of cells from posterior to anterior to a file if the user specified it
    	    print_cell_columns(ip, sd, cls[i], dirname_cons[i], set_num);
    	}
    	if (!sd.no_growth && sd.section == SEC_POST && !ip.short_circuit) { // Copy the concentration levels to the mutant data (if not 	short circuiting)
    	    copy_cl_to_mutant(sd, baby_cls[i], mds[i]);
    	}
    	
    	// Print how the mutant performed and finish book-keeping
    	term->verbose() << "  " << term->blue << "Done: " << term->reset << mds[i].print_name << " scored ";
    	if (passed) {
    	    mds[i].secs_passed[sd.section] = true; // Mark that this mutant has passed this simulation
    	    score += mds[i].tests[sd.section](mds[i], wtfeat);
    	    double max_score = mds[i].max_cond_scores[sd.section];
    	    if (score == max_score) {
    	        // Copy the concentration levels to the mutant data if this is a posterior simulation (if short circuiting)
    	        if (!sd.no_growth && sd.section == SEC_POST && ip.short_circuit) {
    	            copy_cl_to_mutant(sd, baby_cls[i], mds[i]);
    	        }
    	        term->verbose() << term->blue;
    	    } else {
    	        term->verbose() << term->red;
    	    }
    	    term->verbose() << score << " / " << max_score;
    	} else {
    	    term->verbose() << term->red << "failed to complete the simulation";
    	}
    	term->verbose() << term->reset << endl;
		score_array[i]= score;
	}
    return score_array;
}

/* model performs the biological functions of a simulation
	parameters:
 sd: the current simulation's data
 rs: the current simulation's rates
 cl: the concentration levels used for analysis and storage
 baby_cl: the concentration levels used for simulating
 md: the mutant to simulate
	returns: the score of the mutant
	notes:
	todo:
 */

__global__ void firstCUDA(sim_data sd, params p[], array2D<double> rates_active_arr[], con_levels baby_cls[], int baby_j, int j, int time_prev, bool past_induction[],bool past_recovery[]){

	unsigned int k = threadIdx.x;
	unsigned int i = blockIdx.x;
    baby_cls[i].cons[BIRTH][baby_j][k] = baby_cls[i].cons[BIRTH][time_prev][k];
    baby_cls[i].cons[PARENT][baby_j][k] = baby_cls[i].cons[PARENT][time_prev][k];
    

	if (sd.width_current == sd.width_total || k % sd.width_total <= sd.active_start) { // Compute only existing (i.e. already grown) cells
    			// Calculate the cell indices at the start of each mRNA and protein's delay
                int old_cells_mrna[NUM_INDICES];
                int old_cells_protein[NUM_INDICES];
                calculate_delay_indices(sd, baby_cls[i], baby_j, j, k, rates_active_arr[i], old_cells_mrna, old_cells_protein);
                
                // Perform biological calculations
                st_context stc(time_prev, baby_j, k);
                protein_synthesis(sd, rates_active_arr[i], baby_cls[i], stc, old_cells_protein);
                dimer_proteins(sd,rates_active_arr[i], baby_cls[i], stc);
                mRNA_synthesis(sd, rates_active_arr[i], baby_cls[i], stc, old_cells_mrna, p[i], past_induction[i], past_recovery[i]);
    }
	


}

void launchkernel(int cells_total,sim_data& sd, params p_darray[], array2D<double> rates_active_darray[],con_levels baby_cls_darray[], int baby_j,int j,int time_prev,bool past_induction[], bool past_recovery[]){
	//Start timer
	cudaEvent_t start,stop;
	
	//Set dimensions
	dim3 dimBlock(cells_total,1,1); //each cell had own thread
	dim3 dimGrid(6,1,1); //simulation done on single block
	cudaDeviceSetLimit(cudaLimitStackSize, 65536);
	//Run kernel
	//findtauCUDA<<<dimGrid,dimBlock>>>(st,critical,cx,cHOR,cells);
	firstCUDA<<<dimGrid, dimBlock>>>(sd, p_darray, rates_active_darray, baby_cls_darray, baby_j, j, time_prev,past_induction, past_recovery);

	CUDA_ERRCHK(cudaDeviceSynchronize());	
	//convert back to CPU
	if (cudaPeekAtLastError() != cudaSuccess) {
        cout << "Kernel launch error: " << cudaPeekAtLastError() << "\n";
    }

	CUDA_ERRCHK(cudaDeviceSynchronize());
}


bool model (input_params& ip, sim_data& sd, rates_static& rs, rates_dynamic rs_ds[], con_levels cls[], con_levels baby_cls[], mutant_data mds[], double temp_rates[][2]) {
    int steps_elapsed = sd.steps_split; // Used to determine when to split a column of cells
    for (int i = 0; i < ip.num_active_mutants; i++) {
        update_rates(rs, rs_ds[i], sd.active_start); // Update the active rates based on the base rates, perturbations, and gradients
    }
    // Iterate through each time step
    int j; // Absolute time used by cl
    int baby_j; // Cyclical time used by baby_cl
    bool past_induction[ip.num_active_mutants]; // Whether we've passed the point of induction of knockouts or overexpression
    bool past_recovery[ip.num_active_mutants]; // Whether we've recovered from the knockouts or overexpression


    for (int i = 0; i < ip.num_active_mutants; i++) {
	//Copy arrays to GPU
		baby_cls[i].cons.allocateGPU();
		rs_ds[i].rates_active.allocateGPU();
		sd.neighbors.allocateGPU();
		baby_cls[i].allocateGPU();
	
		baby_cls[i].cons.swapToGPU();
		rs_ds[i].rates_active.swapToGPU();
		sd.neighbors.swapToGPU();
		baby_cls[i].swapToGPU();
	}

	array2D<double> rates_active_arr[ip.num_active_mutants];
	for (int i = 0; i < ip.num_active_mutants; i++) {
		rates_active_arr[i]=rs_ds[i].rates_active;
	}

	params p[ip.num_active_mutants];
    for (int i = 0; i < ip.num_active_mutants; i++) {
		p[i].initialize(mds[i]);
	}

	
    for (int i = 0; i < ip.num_active_mutants; i++) {
         past_induction[i]=false;
         past_recovery[i]=false;
    }

	bool *past_induction_darray;
	int size_past_array=ip.num_active_mutants* sizeof(bool);
	CUDA_ERRCHK(cudaMalloc((void**)&past_induction_darray, size_past_array));
	CUDA_ERRCHK(cudaMemcpy(past_induction_darray,past_induction,size_past_array,cudaMemcpyHostToDevice));

    bool *past_recovery_darray;
	CUDA_ERRCHK(cudaMalloc((void**)&past_recovery_darray, size_past_array));
	CUDA_ERRCHK(cudaMemcpy(past_recovery_darray,past_recovery,size_past_array,cudaMemcpyHostToDevice));

	array2D<double> *rates_active_darray;
	int size_array2D=ip.num_active_mutants* sizeof(array2D<double>);
	CUDA_ERRCHK(cudaMalloc((void**)&rates_active_darray, size_array2D));
	CUDA_ERRCHK(cudaMemcpy(rates_active_darray,rates_active_arr,size_array2D,cudaMemcpyHostToDevice));

    con_levels *baby_cls_darray;
	int size_array_babycl=ip.num_active_mutants* sizeof(con_levels);
	CUDA_ERRCHK(cudaMalloc((void**)&baby_cls_darray, size_array_babycl));
	CUDA_ERRCHK(cudaMemcpy(baby_cls_darray,baby_cls,size_array_babycl,cudaMemcpyHostToDevice));

    params *p_darray;
	int size_array_p=ip.num_active_mutants* sizeof(params);
	CUDA_ERRCHK(cudaMalloc((void**)&p_darray, size_array_p));
	CUDA_ERRCHK(cudaMemcpy(p_darray,p,size_array_p,cudaMemcpyHostToDevice));


	



    int time_prev=0;
    for (j = sd.time_start, baby_j = 0; j < sd.time_end; j++, baby_j = WRAP(baby_j + 1, sd.max_delay_size)) {
        if (j % 100 == 0) {
		    cout << "iter " << j << '\n';
        }
		
		for (int i = 0; i < ip.num_active_mutants; i++) {
			// only updates base_rates and cell_rates, which are not directly involved in simulation
			//happeing once each simulation
    	    if (!past_induction[i] && !past_recovery[i] && (j  > anterior_time(sd,mds[i].induction))) {
    	        //cout<<anterior_time(sd,md.induction)<<endl;
    	        knockout(rs_ds[i], mds[i], 1); //knock down rates after the induction point
    	        perturb_rates_all(rs,rs_ds[i]); //This is used for knockout the rate in the existing cells, may need modification
    	        past_induction[i] = true;
    	    }
			
    	    //if (past_induction && (j + sd.steps_til_growth > md.recovery +sd.max_delay_size)) {
    	    if (past_induction[i] && (j > anterior_time(sd,mds[i].recovery))) {
    	        revert_knockout(rs_ds[i], mds[i], temp_rates[i]);
    	        past_recovery[i] = true;
    	    }
    	    
    	    
    	    
    	    time_prev = WRAP(baby_j - 1, sd.max_delay_size); // Time is cyclical, so time_prev may not be baby_j - 1
			// try to put this part into GPU calculation
    	    
    	    // Iterate through each extant cell
			//generate params struct with values
		}

		
	    
        

		launchkernel(sd.cells_total,sd, p_darray, rates_active_darray, baby_cls_darray, baby_j,j, time_prev, past_induction_darray, past_recovery_darray);
       
		for (int i = 0; i < ip.num_active_mutants; i++) {
        // Split cells periodically in anterior simulations
	        if (sd.section == SEC_ANT && (steps_elapsed % sd.steps_split) == 0) {
				baby_cls[i].cons.swapToCPU();
				rs_ds[i].rates_active.swapPointerToCPU();
    	        split(sd, rs, rs_ds[i], baby_cls[i], baby_j, j);
    	        update_rates(rs, rs_ds[i], sd.active_start);
    	        steps_elapsed = 0;
				rs_ds[i].rates_active.swapToGPU();
				baby_cls[i].cons.swapToGPU();
    	    }
    	    
			// all operations can be executed on CPU, data already on CPU
    	    // Update the active record data and split counter
    	    steps_elapsed++;
			// Record of the start of the active PSM at each time step
			baby_cls[i].swapPointerToCPU();
    	    baby_cls[i].active_start_record[baby_j] = sd.active_start;
    	    baby_cls[i].active_end_record[baby_j] = sd.active_end;
			baby_cls[i].swapToGPU();
			// unilateral copy from baby_cl to cl, which means only form GPU to CPU
    	    // Copy from the simulating cl to the analysis cl
    	    if (baby_j % (sd.max_delay_size)  == 0) {
				baby_cls[i].cons.swapToCPU();
				baby_cls[i].swapPointerToCPU();
				for (int l =0, m=j-sd.max_delay_size; l<sd.max_delay_size;l++,m++ ){
    		       		baby_to_cl(baby_cls[i], cls[i], l, m / sd.big_gran);
				}
				baby_cls[i].cons.swapPointerToGPU();
				baby_cls[i].swapPointerToGPU();
    	    }
			if (j==sd.time_end-2) {
				baby_cls[i].cons.swapToCPU();
				baby_cls[i].swapPointerToCPU();
				for (int l =0, m=j-baby_j; l<baby_j;l++,m++ ){
    	       		baby_to_cl(baby_cls[i], cls[i], l, m / sd.big_gran);
				}
				baby_cls[i].cons.swapPointerToGPU();
				baby_cls[i].swapPointerToGPU();
    	    }
		}	
    }
	for (int i = 0; i < ip.num_active_mutants; i++) {
		baby_cls[i].cons.swapToCPU();
		rs_ds[i].rates_active.swapToCPU();
		sd.neighbors.swapToCPU();
		baby_cls[i].swapToCPU();

		baby_cls[i].cons.deallocateGPU();
		rs_ds[i].rates_active.deallocateGPU();
		sd.neighbors.deallocateGPU();
		baby_cls[i].deallocateGPU();
    	// Copy the last time step from the simulating cl to the analysis cl and mark where the simulating cl left off time-wise
	    baby_to_cl(baby_cls[i], cls[i], WRAP(baby_j - 1, sd.max_delay_size), (j - 1) / sd.big_gran);
	    sd.time_baby = baby_j;
	}
	return true;
}

/* calculate_delay_indices calculates where the given cell was at the start of all mRNA and protein delays
	parameters:
 sd: the current simulation's data
 cl: the concentration levels used for analysis and storage
 baby_time: the cyclical time used by baby_cl, the cl for simulating
 time: the absolute time used by cl, the cl for analysis
 cell_index: the current cell index
 active_rates: the active rates
 old_cells_mrna: the indices of the cell's old positions for mRNA delays
 old_cells_protein: the indices of the cell's old positions for protein delays
	returns: nothing
	notes:
	todo:
 */
__device__ void calculate_delay_indices (sim_data& sd, con_levels& cl, int baby_time, int time, int cell_index, array2D<double>& active_rates, int old_cells_mrna[], int old_cells_protein[]) {
    if (sd.section == SEC_POST) { // Cells in posterior simulations do not split so the indices never change
        for (int l = 0; l < NUM_INDICES; l++) {
            old_cells_mrna[IMH1 + l] = cell_index;
            old_cells_protein[IPH1 + l] = cell_index;
        }
    } else { // Cells in anterior simulations split so with long enough delays the cell must look to its parent for values, causing its effective index to change over time
        for (int l = 0; l < NUM_INDICES; l++) {
            old_cells_mrna[IMH1 + l] = index_with_splits(sd, cl, baby_time, time, cell_index, active_rates[RDELAYMH1 + l][cell_index]);
            old_cells_protein[IPH1 + l] = index_with_splits(sd, cl, baby_time, time, cell_index, active_rates[RDELAYPH1 + l][cell_index]);
        }
    }
}

/* index_with_splits calculates where the given cell was at the start of the given delay
	parameters:
 sd: the current simulation's data
 cl: the concentration levels used for analysis and storage
 baby_time: the cyclical time used by baby_cl, the cl for simulating
 time: the absolute time used by cl, the cl for analysis
 cell_index: the current cell index
 delay: the amount of time (in minutes) the delay takes
	returns: the index of the cell at the start of the delay
	notes:
	todo:
 */
__device__ inline int index_with_splits (sim_data& sd, con_levels& cl, int baby_time, int time, int cell_index, double delay) {
    int delay_steps = delay / sd.step_size;
    if (time - delay_steps < 0 || time - delay_steps >= cl.cons[BIRTH][baby_time][cell_index]) { // If the delay is longer than the simulation has run or shorter than the cell's age then return the cell's index
        return cell_index;
    } else { // If the delay is not longer than the simulation has run but longer than the cell's age then move to the cell's parent and check again
        return index_with_splits(sd, cl, baby_time, time, cl.cons[PARENT][baby_time][cell_index], delay);
    }
}

/* any_less_than_0 checks if any of the concentrations in the given range at the given time are less than 0
	parameters:
 cl: the concentration levels
 time: the time step
	returns: true if any concentrations are less than 0, false otherwise
	notes:
 This function checks only the first cell's concentrations for the sake of speed
	todo:
 */
inline bool any_less_than_0 (con_levels& cl, int time) {
    for (int i = MIN_CON_LEVEL; i <= MAX_CON_LEVEL; i++) {
        if (cl.cons[i][time][0] < 0) { // This checks only the first cell
            return true;
        }
    }
    return false;
}

/* concentrations_too_high checks if any concentrations are higher than the given maximum threshold
	parameters:
 cl: the concentration levels used for analysis and storage
 time: the time step
 max_con_thresh: the maximum concentration threshold
	returns: true if any concentrations are too high, false otherwise
	notes:
 This function checks only the first cell's concentrations for the sake of speed
	todo:
 */
inline bool concentrations_too_high (con_levels& cl, int time, double max_con_thresh) {
    if (max_con_thresh != INFINITY) {
        for (int i = MIN_CON_LEVEL; i <= MAX_CON_LEVEL; i++) {
            if (cl.cons[i][time][0] > max_con_thresh) { // This checks only the first cell
                return true;
            }
        }
    }
    return false;
}

/* split splits the posterior-most column of cells
	parameters:
 sd: the current simulation's data
 rs: the current simulation's rates
 cl: the concentration levels used for simulating
 baby_time: the cyclical time used by baby_cl, the cl for simulating
 time: the absolute time used by cl, the cl for analysis
	returns: nothing
	notes:
	todo:
 */
void split (sim_data& sd, rates_static& rs, rates_dynamic& rs_d, con_levels& cl, int baby_time, int time) {
    // Calculate the next active start and current width
    int next_active_start = (sd.active_start + 1) % sd.width_total;
    sd.width_current = MIN(sd.width_current + 1, sd.width_total);
    
    // Determine which child corresponds to which parent
    int parents[sd.height];
    for (int i = 0; i < sd.height; i++) {
        parents[i] = -1;
    }
    for (int i = 0; i < sd.height; i++) {
        int index;
        bool dup;
        do { // Ensure each parent produces exactly one child
            dup = false;
            index = random_int(pair<int, int>(0, sd.height - 1));
            for (int j = 0; j < i; j++) {
                if (index == parents[j]) {
                    dup = true;
                    break;
                }
            }
        } while (dup);
        parents[i] = index;
    }
    
    // Transfer each new cell's parent concentration levels to the new cell
    for (int i = 0; i < NUM_CON_LEVELS; i++) {
        for (int k = 0; k < sd.height; k++) {
            cl.cons[i][baby_time][next_active_start + k * sd.width_total] = cl.cons[i][baby_time][sd.active_start + parents[k] * sd.width_total];
        }
    }
    
    // Set each new cell's birth to the current time and store its assigned parent
    for (int k = 0; k < sd.height; k++) {
        cl.cons[BIRTH][baby_time][next_active_start + k * sd.width_total] = time;
        cl.cons[PARENT][baby_time][next_active_start + k * sd.width_total] = sd.active_start + parents[k] * sd.width_total;
    }
    
    // Perturb the new cells and update the active record data
    perturb_rates_column(sd, rs, rs_d, next_active_start);
    sd.active_start = next_active_start;
    sd.active_end = WRAP(sd.active_start - sd.width_current + 1, sd.width_total);
}

/* copy_records copies each cell's birth and parent to the given time step
	parameters:
 sd: the current simulation's data
 cl: the concentration levels used for analysis and storage
 time: the time step
 time_prev: the previous time step
	returns: nothing
	notes:
	todo:
 */
inline void copy_records (sim_data& sd, con_levels& cl, int time, int time_prev) {
    for (int k = 0; k < sd.cells_total; k++) {
        cl.cons[BIRTH][time][k] = cl.cons[BIRTH][time_prev][k];
        cl.cons[PARENT][time][k] = cl.cons[PARENT][time_prev][k];
    }
}

/* update_rates updates the rates_active array in the given rates struct to account for perturbations and gradients
	parameters:
 rs: the current simulation's rates
 active_start: the column at the start of the posterior
	returns: nothing
	notes:
	todo:
 */
void update_rates (rates_static& rs, rates_dynamic& rs_d, int active_start) {
    if (rs.using_gradients) { // If at least one rate has a gradient
        for (int i = 0; i < NUM_RATES; i++) {
            if (rs.has_gradient[i]) { // If this rate has a gradient
                // Iterate through every cell
                for (int k = 0; k < rs.cells; k++) {
                    // Calculate the cell's index relative to the active start
                    int col = k % rs.width;
                    int gradient_index;
                    if (col <= active_start) {
                        gradient_index = active_start - col;
                    } else {
                        gradient_index = active_start + rs.width - col;
                    }
                    
                    // Set the cell's active rate to its perturbed rate modified by its position's gradient factor
                    rs_d.rates_active[i][k] = rs_d.rates_cell[i][k] * rs.factors_gradient[i][gradient_index];
                }
            } else { // If this rate does not have a gradient then set every cell's active rate to its perturbed rate
                for (int k = 0; k < rs.cells; k++) {
                    rs_d.rates_active[i][k] = rs_d.rates_cell[i][k];
                }
            }
        }
    } else { // If no rates have gradients then set every cell's active rate to its perturbed rate
        for (int i = 0; i < NUM_RATES; i++) {
            for (int k = 0; k < rs.cells; k++) {
                rs_d.rates_active[i][k] = rs_d.rates_cell[i][k];
            }
        }
    }
}

/* protein_synthesis calculates the concentrations of every protein for a given cell
	parameters:
 sd: the current simulation's data
 rs: the active rates
 cl: the concentration levels for simulating
 stc: the spatiotemporal context, i.e. cell and time steps
 old_cells_protein: an array of the cell's indices at the start of each protein's delay
	returns: nothing
	notes:
 Every dimerizing gene must be added to every other dimerizing gene's section in this function.
	todo:
 
	151221: Added prtein synthesis for mespa and mespb
 */
//const cph_indices MESPB_INDICES(CMMESPB, CPMESPB, CPMESPBMESPB, RPSMESPB, RPDMESPB, RDAMESPBMESPB, RDDIMESPBMESPB, RDELAYPMESPB, IMESPB, IPMESPB);
__device__ void protein_synthesis (sim_data& sd,array2D<double>& rs, con_levels& cl, st_context& stc, int old_cells_protein[]) {
    double dimer_effects[NUM_HER_INDICES] = {0}; // Heterodimer calculations
    di_args dia(rs, cl, stc, dimer_effects); // WRAPper for repeatedly used structs
    cp_args cpa(sd, rs, cl, stc, old_cells_protein, dimer_effects); // WRAPper for repeatedly used indices
    
    /// Dimerizing genes
    
    // Her1
    //dim_int(dia, di_indices(CPH1, CPH7,  CPH1H7,  RDAH1H7,  RDDIH1H7,  IH1));
    //dim_int(dia, di_indices(CPH1, CPH13, CPH1H13, RDAH1H13, RDDIH1H13, IH1));
    con_protein_her(cpa, cph_indices(CMH1, CPH1, CPH1H1, RPSH1, RPDH1, RDAH1H1, RDDIH1H1, RDELAYPH1, IH1, IPH1));
    
    
    
    if (sd.section == SEC_ANT) {
        // MespA
        dim_int(dia, di_indices(CPMESPA, CPMESPB, CPMESPAMESPB, RDAMESPAMESPB, RDDIMESPAMESPB, IMESPA));
        con_protein_her(cpa, cph_indices(CMMESPA, CPMESPA, CPMESPAMESPA, RPSMESPA, RPDMESPA, RDAMESPAMESPA, RDDIMESPAMESPA, RDELAYPMESPA, IMESPA, IPMESPA));
        
        // MespB
        dim_int(dia, di_indices(CPMESPB, CPMESPA, CPMESPAMESPB, RDAMESPAMESPB, RDDIMESPAMESPB, IMESPB));
        con_protein_her(cpa, cph_indices(CMMESPB, CPMESPB, CPMESPBMESPB, RPSMESPB, RPDMESPB, RDAMESPBMESPB, RDDIMESPBMESPB, RDELAYPMESPB, IMESPB, IPMESPB));
    }
    
    
    // Delta
    con_protein_delta(cpa, cpd_indices(CMDELTA, CPDELTA, RPSDELTA, RPDDELTA, RDELAYPDELTA, IPDELTA));
}

/* dim_int calculates the dimer interactions for a given protein (a step in protein_synthesis)
	parameters:
 a: a struct containing the arguments needed
 dii: a struct containing the indices needed
	returns: nothing
	notes:
	todo:
 */
__device__ inline void dim_int (di_args& a, di_indices dii) {
    array2D<double>& r = a.rs;
    concentration_level<double>& c = a.cl.cons;
    int tp = a.stc.time_prev;
    int cell = a.stc.cell;
    
    // The part of the given Her protein concentration's differential equation that accounts for heterodimers
    a.dimer_effects[dii.dimer_effect] =
    a.dimer_effects[dii.dimer_effect]
    - r[dii.rate_association][cell] * c[dii.con_protein_self][tp][cell] * c[dii.con_protein_other][tp][cell]
    + r[dii.rate_dissociation][cell] * c[dii.con_dimer][tp][cell];
}

/* con_protein_her calculates the protein concentration of the given Her gene
	parameters:
 a: a struct containing the arguments needed
 dii: a struct containing the indices needed
	returns: nothing
	notes:
	todo:
 */
__device__ inline void con_protein_her (cp_args& a, cph_indices i) {
    array2D<double>& r = a.rs;
    concentration_level<double>& c = a.cl.cons;
    int cell = a.stc.cell;
    int delay_steps = r[i.delay_protein][cell] / a.sd.step_size;
    int tc = a.stc.time_cur;
    int tp = a.stc.time_prev;
    int td = WRAP(tc - delay_steps, a.sd.max_delay_size);
    
    // The part of the given Her protein concentration's differential equation that accounts for everything but heterodimers, whose influence is calculated in dim_int
    c[i.con_protein][tc][cell] =
    c[i.con_protein][tp][cell]
    + a.sd.step_size *
    (r[i.rate_synthesis][cell] * c[i.con_mrna][td][a.old_cells[i.old_cell]]
     - r[i.rate_degradation][cell] * c[i.con_protein][tp][cell]
     - 2 * r[i.rate_association][cell] * SQUARE(c[i.con_protein][tp][cell])
     + 2 * r[i.rate_dissociation][cell] * c[i.con_dimer][tp][cell]
     + a.dimer_effects[i.dimer_effect]);
}

/* con_protein_delta calculates the protein concentration of the given Delta gene
	parameters:
 a: a struct containing the arguments needed
 dii: a struct containing the indices needed
	returns: nothing
	notes:
	todo:
 */
__device__ inline void con_protein_delta (cp_args& a, cpd_indices i) {
    array2D<double>& r = a.rs;
    concentration_level<double>& c = a.cl.cons;
    int cell = a.stc.cell;
    int delay_steps = r[i.delay_protein][cell] / a.sd.step_size;
    int tc = a.stc.time_cur;
    int tp = a.stc.time_prev;
    int td = WRAP(tc - delay_steps, a.sd.max_delay_size);
    
    // The Delta protein concentration's differential equation (no dimerization occurs)
    c[i.con_protein][tc][cell] =
    c[i.con_protein][tp][cell]
    + a.sd.step_size * (r[i.rate_synthesis][cell] * c[i.con_mrna][td][a.old_cells[i.old_cell]]
                        - r[i.rate_degradation][cell] * c[i.con_protein][tp][cell]);
}

/* dimer_proteins calculates the concentrations of every dimer for a given cell
	parameters:
 sd: the current simulation's data
 rs: the active rates
 cl: the concentration levels for simulating
 stc: the spatiotemporal context, i.e. cell and time steps
	returns: nothing
	notes:
	todo:
 
	151221: added dimerization for mespamespa, mespamespb, mespbmespb
 */
__device__ void dimer_proteins (sim_data& sd, array2D<double>& rs, con_levels& cl, st_context& stc) {
    cd_args cda(sd, rs, cl, stc); // WRAPper for repeatedly used structs
    
    
    con_dimer(cda, CPH1H1, 0, cd_indices(CPH1, RDAH1H1, RDDIH1H1, RDDGH1H1));
    
    if (sd.section == SEC_ANT) {
        for (int i = CPMESPAMESPA, j = 0;   i <= CPMESPAMESPB; i++, j++) {
            con_dimer(cda, i, j, cd_indices(CPMESPA, RDAMESPAMESPA, RDDIMESPAMESPA, RDDGMESPAMESPA));
        }
        con_dimer(cda, CPMESPBMESPB, 0, cd_indices(CPMESPB, RDAMESPBMESPB, RDDIMESPBMESPB, RDDGMESPBMESPB));
    }
    
}

/* con_dimer calculates the given dimer's concentration
	parameters:
 a: a struct containing the arguments needed
 con: the index of the dimer
 offset: an index offset for reusability purposes
 i: a struct containing the indices needed
	returns: nothing
	notes:
	todo:
 TODO consolidate these parameters
 
	151221: pay attention to the index for mesp genes
 */
__device__ inline void con_dimer (cd_args& a, int con, int offset, cd_indices i) {
    array2D<double>& r = a.rs;
    concentration_level<double>& c = a.cl.cons;
    int tc = a.stc.time_cur;
    int tp = a.stc.time_prev;
    int cell = a.stc.cell;
    int con_offset=offset;
    if (i.con_protein == CPH1 && offset == 2){
        con_offset=4;
    }
    
    // The given dimer concentration's differential equation
    c[con][tc][cell] =
    c[con][tp][cell]
    + a.sd.step_size * (r[i.rate_association + offset][cell] * c[i.con_protein][tp][cell] * c[i.con_protein + con_offset][tp][cell]
                        - r[i.rate_dissociation + offset][cell] * c[con][tp][cell]
                        - r[i.rate_degradation + offset][cell] * c[con][tp][cell]);
}

/* mRNA_synthesis calculates the concentrations of every mRNA for a given cell
	parameters:
 sd: the current simulation's data
 rs: the active rates
 cl: the concentration levels for simulating
 stc: the spatiotemporal context, i.e. cell and time steps
 old_cells_mrna: an array of the cell's indices at the start of each mRNA's delay
 md: the currently simulating mutant's data
	returns: nothing
	notes:
	todo:
 
	151221: Added mRNA transcription for meps genes, pay attention to index of mesp genes
 */
__device__ void mRNA_synthesis (sim_data& sd, array2D<double>& rs, con_levels& cl, st_context& stc, int old_cells_mrna[], params& p, bool past_induction, bool past_recovery ) {
    // Translate delays from minutes to time steps
    int delays[NUM_INDICES];
    for (int j = 0; j < NUM_INDICES; j++) {
        delays[j] = rs[RDELAYMH1 + j][stc.cell] / sd.step_size;
    }
    
    // Calculate the influence of the given cell's neighbors (via Delta-Notch signaling)
    double avg_delays[NUM_DD_INDICES]; // Averaged delays for each mRNA concentration caused by the given cell's neighbors' Delta protein concentrations
    if (sd.height == 1) { // For 2-cell and 1D simulations
        if (sd.width_current > 2) { // For 1D simulations
            // Each cell has 2 neighbors so calculate where they and the active start and end were at the start of each mRNA concentration's delay
            int neighbors[NUM_DD_INDICES][NEIGHBORS_1D];
            for (int j = 0; j < NUM_INDICES; j++) {
                int old_time = WRAP(stc.time_cur - delays[j], sd.max_delay_size);
                int old_active_start = cl.active_start_record[old_time];
                int old_active_end = cl.active_end_record[old_time];
                calc_neighbors_1d(sd, neighbors[IMH1 + j], old_cells_mrna[IMH1 + j], old_active_start, old_active_end);
            }
            
            // For each mRNA concentration, average the give cell's neighbors' Delta protein concentrations
            for (int j = 0; j < NUM_DD_INDICES; j++) {
                int* cells = neighbors[IMH1 + j];
                int cell = old_cells_mrna[IMH1 + j];
                int time = WRAP(stc.time_cur - delays[j], sd.max_delay_size);
                
                if (cell % sd.width_total == cl.active_start_record[time]) {
                    avg_delays[IMH1 + j] = cl.cons[CPDELTA][time][cells[0]];
                } else if (cell % sd.width_total == cl.active_end_record[time]) {
                    avg_delays[IMH1 + j] = cl.cons[CPDELTA][time][cells[1]];
                } else {
                    avg_delays[IMH1 + j] = (cl.cons[CPDELTA][time][cells[0]] + cl.cons[CPDELTA][time][cells[1]]) / 2;
                }
            }
        } else { // For 2-cell simulations
            // Both cells have one neighbor each so no averaging is required
            for (int j = 0; j < NUM_DD_INDICES; j++) {
                avg_delays[IMH1 + j] = cl.cons[CPDELTA][WRAP(stc.time_cur - delays[j], sd.max_delay_size)][1 - old_cells_mrna[IMH1 + j]];
            }
        }
    } else { // For 2D simulations
        // Each cell has 6 neighbors so calculate where they and the active start and end were at the start of each mRNA concentration's delay
        //int neighbors[NUM_DD_INDICES][NEIGHBORS_2D];
        //for (int j = 0; j < NUM_DD_INDICES; j++) {
            // 2D neighbors are precalculated and simply copied from the structure as needed
            //int cell = old_cells_mrna[IMH1 + j];
            //memcpy(neighbors[IMH1 + j], sd.neighbors[cell], sizeof(int) * NEIGHBORS_2D);
			
        //}
        
        // For each mRNA concentration, average the given cell's neighbors' Delta protein concentrations
        for (int j = 0; j < NUM_DD_INDICES; j++) {


            int cell = old_cells_mrna[IMH1 + j];
			array2D<int>::cell cells = sd.neighbors[cell];
            //int* cells = neighbors[IMH1 + j];
            int time = WRAP(stc.time_cur - delays[j], sd.max_delay_size);
            concentration_level<double>::cell cur_cons = cl.cons[CPDELTA][time];
            double sum=0;
            if (cell % sd.width_total == cl.active_start_record[time]) {
                sum = (cur_cons[cells[0]] + cur_cons[cells[3]] + cur_cons[cells[4]] + cur_cons[cells[5]]) / 4;
            } else if (cell % sd.width_total == cl.active_start_record[time]) {
                sum = (cur_cons[cells[0]] + cur_cons[cells[1]] + cur_cons[cells[2]] + cur_cons[cells[3]]) / 4;
            } else {
                sum = (cur_cons[cells[0]] + cur_cons[cells[1]] + cur_cons[cells[2]] + cur_cons[cells[3]] + cur_cons[cells[4]] + cur_cons[cells[5]]) / 6;
            }
            avg_delays[IMH1 + j] = sum;
        }
    }
    
    // Calculate every mRNA concentration
    for (int j = 0; j < NUM_INDICES; j++) {
        double mtrans;
        //if (j == IMH13) { // her13 mRNA is not affected by dimers' repression
        //	mtrans = rs[RMSH13][stc.cell];
        //} else {
        double avgpd;
        if (j >= IMH1 && j <= IMMESPB) {
            avgpd = avg_delays[IMH1 + j];
        }
        else { // delta mRNA is not affected by Delta-Notch signaling
            avgpd = 0;
        }
        
        
        double oe = 0;
        /*
         if (past_induction && !past_recovery && (j == md.overexpression_rate)) {
         oe = md.overexpression_factor;
         //cout<<md.overexpression_factor<<endl;
         }
         */
        
        if (j==p.overexpression_rate){
            //cout<< md.overexpression_rate<<endl;
            if (past_induction) {
                //cout<<"PAST"<<endl;
                if (!past_recovery) {
                    oe = p.overexpression_factor;
                    //cout<< oe<<endl;
                }
            }
            
        }
        
        if (j == IMMESPA && sd.section == SEC_ANT) {
            mtrans = transcription_mespa(rs, cl, WRAP(stc.time_cur - delays[j], sd.max_delay_size), old_cells_mrna[IMH1 + j], avgpd, rs[RMSH1 + j][stc.cell], oe, sd.section);
            //cout<<"mespa"<<mtrans<<endl;
        } else if (j == IMMESPB && sd.section == SEC_ANT) {
            mtrans = transcription_mespb(rs, cl, WRAP(stc.time_cur - delays[j], sd.max_delay_size), old_cells_mrna[IMH1 + j], avgpd, rs[RMSH1 + j][stc.cell], oe);
            
        } else {
            mtrans = transcription(rs, cl, WRAP(stc.time_cur - delays[j], sd.max_delay_size), old_cells_mrna[IMH1 + j], avgpd, rs[RMSH1 + j][stc.cell], oe);
        }
        
        //}
        
        // The current mRNA concentration's differential equation
        cl.cons[CMH1 + j][stc.time_cur][stc.cell] =
        cl.cons[CMH1 + j][stc.time_prev][stc.cell]
        + sd.step_size * (mtrans - rs[RMDH1 + j][stc.cell] * cl.cons[CMH1 + j][stc.time_prev][stc.cell]);
    }
}

/* transcription calculates mRNA transcription, taking into account the effects of dimer repression
	parameters:
 rs: the active rates
 cl: the concentration levels for simulating
 time: the current time
 cell: the current cell
 avgpd: the neighbor-averaged Delta protein concentration
 ms: the rate of mRNA synthesis
 oe: the rate of mRNA overexpression
	returns: the transcription factor
	notes:
	todo:
 TODO clean up these parameters
 */
__device__ inline double transcription (array2D<double>& rs, con_levels& cl, int time, int cell, double avgpd, double ms, double oe) {
    double th1h1= 0, tdelta;
    th1h1 = rs[RCRITPH1H1][cell] == 0 ? 0 : cl.cons[CPH1H1][time][cell] / rs[RCRITPH1H1][cell];
    tdelta = rs[RCRITPDELTA][cell] == 0 ? 0 : avgpd / rs[RCRITPDELTA][cell];
    return ms * (oe + (1 + tdelta) / (1 + tdelta + SQUARE(th1h1)  ));
}

/* 151221: transcription_mespa calculates mRNA transcription for mespa, taking into account the effects of dimer repression
	parameters:
 rs: the active rates
 cl: the concentration levels for simulating
 time: the current time
 cell: the current cell
 avgpd: the neighbor-averaged Delta protein concentration
 ms: the rate of mRNA synthesis
 oe: the rate of mRNA overexpression
	returns: the transcription factor
	notes:
	todo:
 TODO clean up these parameters
 */
__device__ inline double transcription_mespa (array2D<double>& rs, con_levels& cl, int time, int cell, double avgpd, double ms, double oe, int section) {
    double th1h1,  tmespbmespb = 0, tdelta;
    th1h1 = rs[RCRITPH1H1][cell] == 0 ? 0 : cl.cons[CPH1H1][time][cell] / rs[RCRITPH1H1][cell];
    tmespbmespb = rs[RCRITPMESPBMESPB][cell] == 0 ? 0 : cl.cons[CPMESPBMESPB][time][cell] / (rs[RCRITPMESPBMESPB][cell]);
    tdelta = rs[RCRITPDELTA][cell] == 0 ? 0 : avgpd / (rs[RCRITPDELTA][cell]);
    
    return ms * (oe + (tdelta) / (tdelta + rs[NS1][cell] * SQUARE(th1h1)  + SQUARE(tmespbmespb)));
}

/* 151221: transcription_mespb calculates mRNA transcription for mespb, taking into account the effects of dimer repression
	parameters:
 rs: the active rates
 cl: the concentration levels for simulating
 time: the current time
 cell: the current cell
 avgpd: the neighbor-averaged Delta protein concentration
 ms: the rate of mRNA synthesis
 oe: the rate of mRNA overexpression
	returns: the transcription factor
	notes:
	todo:
 TODO clean up these parameters
 */
__device__ inline double transcription_mespb (array2D<double>& rs, con_levels& cl, int time, int cell, double avgpd, double ms, double oe) {
    double tmespamespa = 0, tmespamespb = 0, tmespbmespb = 0, tdelta;
    //th1h1 = rs[RCRITPH1H1][cell] == 0 ? 0 : cl.cons[CPH1H1][time][cell] / rs[RCRITPH1H1][cell];
    //th7h13 = rs[RCRITPH7H13][cell] == 0 ? 0 : cl.cons[CPH7H13][time][cell] / rs[RCRITPH7H13][cell];
    //if (section == SEC_ANT) {
    tmespamespa = rs[RCRITPMESPAMESPA][cell] == 0 ? 0 : cl.cons[CPMESPAMESPA][time][cell] / (rs[RCRITPMESPAMESPA][cell]);
    tmespamespb = rs[RCRITPMESPAMESPB][cell] == 0 ? 0 : cl.cons[CPMESPAMESPB][time][cell] / (rs[RCRITPMESPAMESPB][cell]);
    tmespbmespb = rs[RCRITPMESPBMESPB][cell] == 0 ? 0 : cl.cons[CPMESPBMESPB][time][cell] / (rs[RCRITPMESPBMESPB][cell]);
    //}
    tdelta = rs[RCRITPDELTA][cell] == 0 ? 0 : rs[NS2][cell] * avgpd / (rs[RCRITPDELTA][cell]);
    //avg_delta+= tdelta;
    //avg_rest+=(SQUARE(tmespamespa) + SQUARE(tmespbmespb));
    //cout<<"Tdelta: "<< tdelta<<" REST: "<<(SQUARE(tmespamespa) + SQUARE(tmespamespb) + SQUARE(tmespbmespb))<<endl;
    //cout<<"OE in mespbover: "<<oe<<endl;
    
    //+ SQUARE(tmespamespb)
    return ms * (oe + (1 + tdelta) / (1 + tdelta + SQUARE(tmespamespa) + SQUARE(tmespamespb)+SQUARE(tmespbmespb)));
}

/* calc_neighbors_1d calculates a given cell's neighbors in a 1D simulation
	parameters:
 sd: the current simulation's data
 neighbors: the 2-element array of neighbor indices to fill
 cell: the index of the cell whose neighbors should be calculated
 active_start: the start of the active PSM
 active_end: the end of the active PSM
	returns: nothing
	notes:
	todo:
 */
__device__ void calc_neighbors_1d (sim_data& sd, int neighbors[], int cell, int active_start, int active_end) {
    if (cell == active_start) {
        neighbors[0] = cell - 1;
        neighbors[1] = cell - 1;
    } else if (cell == active_end) {
        neighbors[0] = cell + 1;
        neighbors[1] = cell + 1;
    } else {
        neighbors[0] = cell - 1;
        neighbors[1] = cell + 1;
    }
    
    neighbors[0] = WRAP(neighbors[0], sd.width_total);
    neighbors[1] = WRAP(neighbors[1], sd.width_total);
}

/* calc_neighbors_2d calculates a given cell's neighbors in a 2D simulation
	parameters:
 sd: the current simulation's data
	returns: nothing
	notes:
 2D simulations use a hexagonal grid of cells indexed like this:
 ____  ____  ____  ____
 /    \/    \/    \/    \
 | 0  || 1  || 2  || 3  |
 \____/\____/\____/\____/__
 /    \/    \/    \/    \
 | 4  || 5  || 6  || 7  |
 __\____/\____/\____/\____/
 /    \/    \/    \/    \
 | 8  || 9  || 10 || 11 |
 \____/\____/\____/\____/__
 /    \/    \/    \/    \
 | 12 || 13 || 14 || 15 |
 \____/\____/\____/\____/
 
 This function should be called only when necessary due to the time cost; the populated neighbors array should be reused until invalid.
	todo:
 */
void calc_neighbors_2d (sim_data& sd) {
    for (int i = 0; i < sd.cells_total; i++) {
        if (i % 2 == 0) {																		// All even column cells
            sd.neighbors[i][0] = (i - sd.width_total + sd.cells_total) % sd.cells_total;			// Top
            sd.neighbors[i][1] = (i - sd.width_total + 1 + sd.cells_total) % sd.cells_total;		// Top-right
            sd.neighbors[i][2] = (i + 1) % sd.cells_total;											// Bottom-right
            sd.neighbors[i][3] = (i + sd.width_total) % sd.cells_total;								// Bottom
            if (i % sd.width_total == 0) {														// Left edge
                sd.neighbors[i][4] = i + sd.width_total - 1;										// Bottom-left
                sd.neighbors[i][5] = (i - 1 + sd.cells_total) % sd.cells_total;						// Top-left
            } else {																			// Not a left edge
                sd.neighbors[i][4] = i - 1;															// Bottom-left
                sd.neighbors[i][5] = (i - sd.width_total - 1 + sd.cells_total) % sd.cells_total;	// Top-left
            }
        } else {																				// All odd column cells
            sd.neighbors[i][0] = (i - sd.width_total + sd.cells_total) % sd.cells_total;			// Top
            if (i % sd.width_total == sd.width_total - 1) {											// Right edge
                sd.neighbors[i][1] = i - sd.width_total + 1;										// Top-right
                sd.neighbors[i][2] = (i + 1) % sd.cells_total;										// Bottom-right
            } else {																			// Not a right edge
                sd.neighbors[i][1] = i + 1;															// Top-right
                sd.neighbors[i][2] = (i + sd.width_total + 1 + sd.cells_total) % sd.cells_total;	// Nottom-right
            }																					// All odd column cells
            sd.neighbors[i][3] = (i + sd.width_total) % sd.cells_total;								// Bottom
            sd.neighbors[i][4] = (i + sd.width_total - 1) % sd.cells_total;							// Bottom-left
            sd.neighbors[i][5] = (i - 1 + sd.cells_total) % sd.cells_total;							// Top-left
        }
    }
}

/* perturb_rates_all perturbs every rate for every cell
	parameters:
 rs: the current simulation's rates
	returns: nothing
	notes:
	todo:
 */
void perturb_rates_all (rates_static& rs, rates_dynamic& rs_d) {
    for (int i = 0; i < NUM_RATES; i++) {
        if (rs.factors_perturb[i] == 0) { // If the current rate has no perturbation factor then set every cell's rate to the base rate
            for (int j = 0; j < rs.cells; j++) {
              
                rs_d.rates_cell[i][j] = rs_d.rates_base[i];
            }
        } else { // If the current rate has a perturbation factor then set every cell's rate to a randomly perturbed positive or negative variation of the base with a maximum perturbation up to the rate's perturbation factor
            for (int j = 0; j < rs.cells; j++) {
                rs_d.rates_cell[i][j] = rs_d.rates_base[i] * random_perturbation(rs.factors_perturb[i]);
                
            }
        }
    }
}

/* perturb_rates_column perturbs every rate for every cell in the given column
	parameters:
 sd: the current simulation's data
 rs: the current simulation's rates
 column: the column of cells to perturb
	returns: nothing
	notes:
	todo:
 */
void perturb_rates_column (sim_data& sd, rates_static& rs, rates_dynamic& rs_d, int column) {
    for (int i = 0; i < NUM_RATES; i++) {
        if (rs.factors_perturb[i] != 0) { // Alter only rates with a perturbation factor
            for (int j = 0; j < sd.height; j++) {
                rs_d.rates_cell[i][j * sd.width_total + column] = rs_d.rates_base[i] * random_perturbation(rs.factors_perturb[i]);
            }
        }
    }
}

/* random_perturbation calculates a random perturbation with a maximum absolute change up to the given perturbation factor
	parameters:
 perturb: the perturbation factor
	returns: the random perturbation
	notes:
	todo:
 */
inline double random_perturbation (double perturb) {
    return random_double(pair<double, double>(1 - perturb, 1 + perturb));
}

/* baby_to_cl copies the data from the given time step in baby_cl to cl
	parameters:
 baby_cl: the concentration levels for simulating
 cl: the concentration levels for analysis and storage
 baby_time: the time step to access baby_cl with
 time: the time step to access cl with
	returns: nothing
	notes:
	todo:
 */
void baby_to_cl (con_levels& baby_cl, con_levels& cl, int baby_time, int time) {
    for (int i = 0; i < cl.num_con_levels; i++) {
        for (int k = 0; k < cl.cells; k++) {
            cl.cons[i][time][k] = baby_cl.cons[i][baby_time][k];
        }
    }
    cl.active_start_record[time] = baby_cl.active_start_record[baby_time];
    cl.active_end_record[time] = baby_cl.active_end_record[baby_time];
}

/* anterior_time converts the given time step to its equivalent in anterior time
	parameters:
 sd: the current simulation's data
 time: the time step
	returns: the converted time step
	notes:
	todo:
 */
int anterior_time (sim_data& sd, int time) {
    return time / sd.big_gran - sd.steps_til_growth / sd.big_gran + sd.max_delay_size;
}

