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

/*
feats.cpp contains functions to analyze and test the oscillation features of simulations.
*/

#include <cfloat> // Needed for DBL_MAX
#include <algorithm>
#include "feats.hpp" // Function declarations

#include "io.hpp"
#include "sim.hpp"

using namespace std;

extern terminal* term; // Declared in init.cpp

/*

In the future need to find a way to unify the functions in this file into some more consistent ones.
For the time being, I'm focusing on obtaining the needed data.

*/
int get_peaks_and_troughs1 (sim_data& sd, con_levels& cl, int actual_cell, int time_start, growin_array& crit_points, growin_array& type, growin_array& position, int mr) {
	/*
	Calculates all the peaks and troughs for the her1 mRNA oscillations in a cell within a time range. The cell number is given as a
	number relative to the entire PSM.
	*/
	int num_points = 0;
	int col = actual_cell % sd.width_total;
    //cout<< "time start :  " <<time_start<<endl;
	concentration_level<double>::timespan conc = cl.cons[mr];
	for (int j = time_start + 1; j < sd.time_end - 1 && cl.cons[BIRTH][j][actual_cell] == cl.cons[BIRTH][j - 1][actual_cell] && cl.cons[BIRTH][j][actual_cell] == cl.cons[BIRTH][j + 1][actual_cell]; j++) {
       		int pos = 0;
		if (cl.active_start_record[j] >= col) {
			pos = cl.active_start_record[j] - col;
		} else {
			pos = cl.active_start_record[j] + sd.width_total - col;
		}
	
		// check if the current point is a peak
		bool is_peak = true;
		for (int k = MAX(j - (2 / sd.step_size / sd.big_gran), time_start); k <= MIN(j + 2 / sd.step_size / sd.big_gran, sd.time_end - 1); k++) {
			if (conc[j][actual_cell] <= conc[k][actual_cell] && j!=k) {
				is_peak = false;
			}
		}
		if (is_peak) {
			crit_points[num_points] = j;
			type[num_points] = 1;
			position[num_points] = pos;
			num_points++;
		}
		
		// check if the current point is a trough
		bool is_trough = true;
		for (int k = MAX(j - (2 / sd.step_size / sd.big_gran), time_start); k <= MIN(j + 2 / sd.step_size / sd.big_gran, sd.time_end - 1); k++) {
			if (conc[j][actual_cell] >= conc[k][actual_cell] && j!=k) {
				is_trough = false;
			}
		}
		if (is_trough) {
			crit_points[num_points] = j;
			type[num_points] = -1;
			position[num_points] = pos;
			num_points++;
		}
	}
	//cout<<max-min<<endl;
	return num_points;
}


int get_peaks_and_troughs2 (sim_data& sd, con_levels& cl, int actual_cell, int time_start, growin_array& crit_points, growin_array& type, growin_array& position, int mr, double* mh1_comp, double* mespa_comp, double* mespb_comp) {
	/*
	151221: record the concentration value of mh1, mespa, mespb in the time specified in order to calculate the complementary expression score of mespa and mespb. 
	*/
	int num_points = 0;
	int col = actual_cell % sd.width_total;
	
	concentration_level<double>::timespan conc = cl.cons[mr];
	int compl_count=0;
	for (int j = time_start + 1; j < sd.time_end - 1 && cl.cons[BIRTH][j][actual_cell] == cl.cons[BIRTH][j - 1][actual_cell] && cl.cons[BIRTH][j][actual_cell] == cl.cons[BIRTH][j + 1][actual_cell]; j++) {
		// calculate position in the PSM of the cell
        
        
        mh1_comp[compl_count]=cl.cons[CMH1][j][actual_cell];                 //record concentration value of mh1 151221
        mespa_comp[compl_count]=cl.cons[CMMESPA][j][actual_cell];           //record concentration value of mespa 151221
        mespb_comp[compl_count]=cl.cons[CMMESPB][j][actual_cell];             //record concentration value of mespb 151221
        compl_count++;
        
		int pos = 0;
		if (cl.active_start_record[j] >= col) {
			pos = cl.active_start_record[j] - col;
		} else {
			pos = cl.active_start_record[j] + sd.width_total - col;
		}
	
		// check if the current point is a peak
		bool is_peak = true;
		for (int k = MAX(j - (2 / sd.step_size / sd.big_gran), time_start); k <= MIN(j + 2 / sd.step_size / sd.big_gran, sd.time_end - 1); k++) {
			if (conc[j][actual_cell] <= conc[k][actual_cell] && j!=k) {
				is_peak = false;
				
			}
		}
		if (is_peak) {
			crit_points[num_points] = j;
			type[num_points] = 1;
			position[num_points] = pos;
			num_points++;
		}
		
		// check if the current point is a trough
		bool is_trough = true;
		for (int k = MAX(j - (2 / sd.step_size / sd.big_gran), time_start); k <= MIN(j + 2 / sd.step_size / sd.big_gran, sd.time_end - 1); k++) {
			if (conc[j][actual_cell] >= conc[k][actual_cell] && j!=k) {
				is_trough = false;
			}
		}
		if (is_trough) {
			crit_points[num_points] = j;
			type[num_points] = -1;
			position[num_points] = pos;
			num_points++;
		}
	}
	//cout<<max-min<<endl;
	return num_points;
}

double test_complementary (sim_data& sd, con_levels& cl, int time, int con1, int con2) {   // calculate the complementary expression score of mespa and mespb. not used
	double avg_row_con1[sd.width_total];
    double avg_row_con2[sd.width_total];
	memset(avg_row_con1, 0, sizeof(double) * sd.width_total);
    memset(avg_row_con2, 0, sizeof(double) * sd.width_total);

	for (int y = 0; y < sd.width_total; y++) {
		for (int x = 0; x < sd.height; x++) {
			int cell = x * sd.width_total + y;
			avg_row_con1[y] += cl.cons[con1][time][cell];
            avg_row_con2[y] += cl.cons[con2][time][cell];
		}
		
		avg_row_con1[y] /= sd.height;
        avg_row_con2[y] /= sd.height;
	}

	return pearson_correlation(avg_row_con1, avg_row_con2, 0.6*sd.width_total, sd.width_total);  //JY WT.10 rounding?
}

double test_compl(sim_data& sd, double* con1, double* con2) {  // 151221: calculate the complementary expression score of mespa and mespb
	return pearson_correlation(con1, con2, (int)(0.6*(sd.width_total*sd.steps_split - 2)),sd.width_total*sd.steps_split - 2);
}



void osc_features_ant (sim_data& sd, input_params& ip, features& wtfeat, char* filename_feats, con_levels& cl, mutant_data& md, int start_line, int end_line, int start_col, int end_col, int set_num) {
	static int con[4] = {CMH1,  CMDELTA, CMMESPA, CMMESPB};
	static int ind[4] = {IMH1,  IMDELTA, IMMESPA, IMMESPB};
	static const char* concs[4] = {"mh1",  "mdelta", "mespa", "mespb"};
	static const char* feat_names[NUM_FEATURES] = {"period", "amplitude", "sync"};
	
	
	growin_array crit_points(sd.steps_total / (20/sd.step_size)); // Array that will hold all the critical points in the graph
	growin_array type(sd.steps_total / (20/sd.step_size)); // Array that will specify whether each critical point is a peak or a trough (-1 for trough, 1 for peak)
	growin_array position(sd.steps_total / (20/sd.step_size)); // Array that will hold the position in the PSM of the critical points
	
	int strlen_set_num = INT_STRLEN(set_num); // How many bytes the ASCII representation of set_num takes
	char* str_set_num = (char*)mallocate(sizeof(char) * (strlen_set_num + 1));
	sprintf(str_set_num, "%d", set_num);
	
	int num_cell = (end_line - start_line) * (end_col - start_col);
	double mh1_comp[sd.width_total*sd.steps_split - 2];   //151221: structure to store concentration value for mh1, used in get_peaks_and_troughs2 
	double mespa_comp[sd.width_total*sd.steps_split - 2]; //151221: structure to store concentration value for mespa, used in get_peaks_and_troughs2
	double mespb_comp[sd.width_total*sd.steps_split - 2]; //151221: structure to store concentration value for mespb, used in get_peaks_and_troughs2
	double comp_score_a = 0; //151221: complementary score for mespa
	double comp_score_b = 0; //151221: complementary score for mespa
	memset(mh1_comp, 0, sizeof(double) * (sd.width_total*sd.steps_split - 2));
	memset(mespa_comp, 0, sizeof(double) * (sd.width_total*sd.steps_split - 2));
	memset(mespb_comp, 0, sizeof(double) * (sd.width_total*sd.steps_split - 2));
	
	
	for (int i = 0; i < 4; i++) {
		ofstream features_files[NUM_FEATURES]; // Array that will hold the files in which to output the period and amplitude
	
		int mr = con[i];
		int index = ind[i];
		if (ip.ant_features) {
			for (int j = 0; j < NUM_FEATURES; j++) {
				char* filename = (char*)mallocate(sizeof(char) * strlen(filename_feats) + strlen("set_") + strlen_set_num + 1 + strlen(feat_names[j]) + 1 + strlen(concs[i]) + strlen("_ant.feats") + 1);
				sprintf(filename, "%sset_%s_%s_%s_ant.feats", filename_feats, str_set_num, feat_names[j], concs[i]);
				cout << "      ";
				open_file(&(features_files[j]), filename, false);
				mfree(filename);
			}
			
			features_files[PERIOD] << sd.height << "," << sd.width_total << endl;
			features_files[AMPLITUDE] << sd.height << "," << sd.width_total << endl;
		}
		if (md.index == MUTANT_WILDTYPE && (mr == CMH1 || mr == CMMESPA)) { // following only done for CMH1 and CMMESPA in wildtype to prevent unnecessary repeating
            double period_avg = 0;
            int time_start= anterior_time(sd, sd.steps_til_growth + (sd.width_total - sd.width_initial) * sd.steps_split); // time after which the PSM is full of cells
            int num_cells_passed = 0;
            //cout<<"time: " << time_start<< " "<<sd.steps_til_growth + (sd.width_total - sd.width_initial) * sd.steps_split<<" "<< sd.time_end<<endl;
            //cout<<"time:  "<<anterior_time(sd,93000)<<endl;
            for (int col = start_col; col < end_col; col++) {						
                for (int line = start_line; line < end_line; line++) {
                    int pos = cl.active_start_record[time_start]; // always looking at cell at position active_start because that is the newest cell
                    
                    int cell = line * sd.width_total + pos;
                    //cout<<"pos: "<<pos<< " cell: "<< cell<<endl;
                    int num_points = 0;
                    if ( mr != CMMESPA) {
                        num_points = get_peaks_and_troughs1(sd, cl, cell, time_start, crit_points, type, position, mr);
                        //cout<<"number of points: "<<num_points<<endl;
                    } else {
                        
                        num_points = get_peaks_and_troughs2(sd, cl, cell, time_start, crit_points, type, position, mr, mh1_comp, mespa_comp, mespb_comp);  // 151221: get_peaks_and_troughs2 records the concentration value of mh1, mespa and mespb and store them in mh1_comp, mespa_comp, mespb_comp
                        comp_score_a+=test_compl(sd, mh1_comp, mespa_comp);
                        comp_score_b+=test_compl(sd, mh1_comp, mespb_comp);
                        //cout<< "num_points: "<< num_points<<endl;

                    } 
                    double periods[num_points];
                    double per_pos[num_points];
                    double per_time[num_points];
                    double amplitudes[num_points];
                    double amp_pos[num_points];
                    
                    memset(periods, 0, sizeof(double) * num_points);
                    memset(amplitudes, 0, sizeof(double) * num_points);
                    memset(per_pos, 0, sizeof(double) * num_points);
                    memset(amp_pos, 0, sizeof(double) * num_points);
                    
                    int pers = 0;
                    int amps = 0;
                    if (num_points >= 3) { 
                        // Calculate all the periods and amplitudes
                        
                        // Find the first peak and the first trough occurring in the graph
                        int cur_point = 0;
                        //cout<<"num_points: "<< num_points<<endl;
                        //cout<<"start: "<<cell << " ";
                        for (; cur_point < num_points; cur_point++) {		
                            // Check for period
                            if (type[cur_point] == 1 && cur_point >= 2) {
                                periods[pers] = (crit_points[cur_point] - crit_points[cur_point - 2]) * sd.step_size * sd.big_gran;
                                per_pos[pers] = position[cur_point - 2] + (position[cur_point] - position[cur_point - 2]) / 2;
                                //if (index==IMMESPB){cout<< periods[pers]<<" "<<per_pos[pers]<<" ";}
                                per_time[pers] = (crit_points[cur_point] - crit_points[cur_point - 2]) / 2 * sd.step_size * sd.big_gran;
                                pers++;
                            }
                        }
                        
                        //cout<<"pers: "<<pers<< "  amps: "<<amps<<endl;
                        bool passed = true;
                        
                        //last period
                        double last_period=0;
                        if (type[num_points-1]==1) {
                           last_period = (crit_points[num_points-1]-crit_points[num_points-3])*sd.step_size*sd.big_gran;
                        }
                        else{
                           last_period = (crit_points[num_points-2]-crit_points[num_points-4])*sd.step_size*sd.big_gran;
                        }
                        
                        //if (index==IMH1){cout<<"last period: "<< last_period<< endl;}
                        
                        double giuratio= last_period/  md.feat.period_post;
                        if (giuratio< 1.4 || giuratio> 2.2) {
                            //cout<< giuratio<<endl;
                            passed=false;
                        }
                        if (passed && mr == CMH1) {
                            num_cells_passed++;
                        }
                        period_avg += (periods[0]);
                    } else {
                        
                        period_avg+= (period_avg < INFINITY ? INFINITY : 0);;
                    }
                    // Printing to files for the anterior features
                    if (ip.ant_features) {
                        for (int j = 0; j < pers; j++) {
                            features_files[PERIOD] << per_pos[j] << ",";
                        }
                        features_files[PERIOD] << endl;
                        for (int j = 0; j < pers; j++) {
                            features_files[PERIOD] << periods[j] << ",";
                        }
                        features_files[PERIOD] << endl;
                        for (int j = 0; j < amps; j++) {
                            features_files[AMPLITUDE] << amp_pos[j] << ",";
                        }
                        features_files[AMPLITUDE] << endl;
                        for (int j = 0; j < amps; j++) {
                            features_files[AMPLITUDE] << amplitudes[j] << ",";
                        }
                        features_files[AMPLITUDE] << endl;
                    }

                    type.reset(sd.steps_total / (20/sd.step_size));
                    crit_points.reset(sd.steps_total / (20/sd.step_size));
                }		
                time_start += sd.steps_split / sd.big_gran; // skip in time until a new column of cells has been formed
            }
            if (ip.ant_features) {
                features_files[PERIOD].close();
                features_files[AMPLITUDE].close();
            }

            
            period_avg /= (end_line - start_line) * (end_col - start_col);
            
            if (index == IMH1) {
                int threshold = 0.7 * (end_line - start_line) * (end_col - start_col); // 151221: originally 80%, now changed to 70%
                //cout<<num_cells_passed<<endl;
                md.conds_passed[SEC_ANT][0] = (num_cells_passed >= threshold);
            }
            if (index == IMMESPA) {
                md.feat.comp_score_ant_mespa = comp_score_a / num_cell;
                md.feat.comp_score_ant_mespb = comp_score_b / num_cell;
            }
		}
		
		if (md.index == MUTANT_WILDTYPE) {
            
            //condition wildtype sync score
            if (index == 0) { //only measure her1
                int time_600 = anterior_time(sd,(600/sd.step_size));
                int time_600_end = anterior_time(sd,(600+30)/sd.step_size);
                //cout<<"WT mh1 sync: ";
                cout<<"WT delta2 mh1 amp: ";
                for (; time_600 < time_600_end; time_600+= 6/sd.step_size){
                    md.feat.sync_score_ant[IMH1]+=ant_sync(sd,cl,CMH1, time_600);                  //WT 4
                    //cout<<ant_sync(sd,cl,CMH1, time_600)<< " ";
                    md.feat.amplitude_ant_time[IMH1][0]+=amp_10(sd, cl, CMH1, time_600);           //delta 2 comparison
                    cout<<amp_10(sd, cl, CMH1, time_600)<< ", ";
                }
                cout<<endl;
                
                int time_630 = anterior_time(sd,(630/sd.step_size));

                int time_three = anterior_time(sd, (840)/sd.step_size);
                int time_three_end = anterior_time(sd, (870)/sd.step_size);
                cout<<"WT dapt2 mh1 amp: ";
                for (;time_three<time_three_end; time_three+=(6/sd.step_size)){
                    md.feat.amplitude_ant_time[IMH1][2]+=amp_10(sd, cl, CMH1, time_three);         //DAPT 2 comparison
                    cout<<amp_10(sd, cl, CMH1, time_three)<< ", ";
                }
                cout<<endl;
                md.feat.sync_score_ant[IMH1]/=5;
                md.feat.amplitude_ant_time[IMH1][0]/=5;
                md.feat.amplitude_ant_time[IMH1][1]=amp_10(sd, cl, CMH1, time_630);
                md.feat.amplitude_ant_time[IMDELTA][1]=amp_10(sd, cl, CMDELTA, time_630);
                md.feat.amplitude_ant_time[IMMESPA][1]=amp_10(sd, cl, CMMESPA, time_630);
                cout<<"WT her1:" <<md.feat.amplitude_ant_time[IMH1][1]<<" "<<md.feat.amplitude_ant_time[IMDELTA][1]<<" "<<md.feat.amplitude_ant_time[IMMESPA][1]<<endl;
                //cout<<"HER1OVER "<<md.feat.amplitude_ant_time[index][1]<< " "<<md.feat.amplitude_ant_time[IMDELTA][1]<< " "<< md.feat.amplitude_ant_time[IMMESPA][1]<<endl;
                md.feat.amplitude_ant_time[IMH1][2]/=5;
            }
            
             
			if (index == 2 || index ==1){
                //complementray score of mespa and mespb
                //for WT calculations
                int time_600 = anterior_time(sd,(600/sd.step_size));
                int time_600_end = anterior_time(sd,(600+30)/sd.step_size);
                for (; time_600 < time_600_end; time_600+= 6/sd.step_size){
                    md.feat.sync_score_ant[index]+=mesp_sync(sd,cl,index+1, time_600);
                }
                md.feat.sync_score_ant[index]/=5;                   //mesp a/b complementary score in sync_score_ant[1/2]
                //cout<<"COMP SCORE: "<< index<< " "<<md.feat.sync_score_ant[index]<<endl;
                
                
			}
            if (index ==1){//mespa
                int time_600 = anterior_time(sd,(600/sd.step_size));
                int time_600_end = anterior_time(sd,(600+30)/sd.step_size);
                cout<<"WT dapt4 mespa amp: ";
                for (; time_600 < time_600_end; time_600+= 6/sd.step_size){
                    md.feat.amplitude_ant_time[IMDELTA][0]+=amp_10(sd, cl, CMMESPA, time_600);            //delta 4 comparison, using amplitude_ant_time[delta][0] to save space
                    md.feat.sync_score_ant[IMDELTA]+= ant_sync(sd, cl, CMMESPB, time_600);          //put wt mespb sync in delta since sync mespa/b are used for 2nd comp scores
                    cout<<amp_10(sd, cl, CMMESPA, time_600)<< ", ";
                }
                cout<<endl;
                
                
                int time_two = anterior_time(sd, (720)/sd.step_size);
                int time_two_end = anterior_time(sd, (750)/sd.step_size);
                cout<<"WT dapt3 mespa amp: ";
                for (;time_two<time_two_end; time_two+=(6/sd.step_size)){
                    md.feat.amplitude_ant_time[IMMESPA][2]+=amp_10(sd, cl, CMMESPA, time_two);    //DAPT 3 comparison
                    cout<<amp_10(sd, cl, CMMESPA, time_two)<< ", ";
                }
                cout<<endl;
                int time_660 = anterior_time(sd,(660/sd.step_size));
                /*
                int time_660_end = anterior_time(sd,(690)/sd.step_size);
                cout<< "WT MESP AMP: ";
                for (; time_660 < time_660_end; time_660+= 6/sd.step_size){
                    md.feat.amplitude_ant_time[IMMESPB][0]+=amp_10(sd, cl, CMMESPB, time_660);            //MESP  comparison
                    md.feat.amplitude_ant_time[IMMESPA][0]+=amp_10(sd, cl, CMMESPA, time_660);            //MESP  comparison
                    cout <<md.feat.amplitude_ant_time[IMMESPB][0]<< " " ;
                }
                cout <<endl;
                */
                md.feat.amplitude_ant_time[IMMESPB][0]=amp_10(sd, cl, CMMESPB, time_660);
                //cout<<md.feat.amplitude_ant_time[IMMESPB][0]<<endl;
                
                md.feat.amplitude_ant_time[IMMESPA][0]=amp_10(sd, cl, CMMESPA, time_660);
                md.feat.amplitude_ant_time[IMMESPA][2]/=5;
                md.feat.sync_score_ant[IMDELTA]/=5;
                cout<<"WT MESPB SNYC: "<<md.feat.sync_score_ant[IMDELTA]<<endl;

                md.feat.amplitude_ant_time[IMDELTA][0]/=5;
                
            }
      
		}

		if (md.index==MUTANT_DELTA){
            if (index==0) {
				int time_600 = anterior_time(sd,(600/sd.step_size));
				int time_600_end = anterior_time(sd,(630/sd.step_size));
                //cout<<"Delta mh1 sync: ";
                //cout<<"Delta mespb sync: ";
                for (;time_600<time_600_end; time_600+=(6/sd.step_size)){
                    md.feat.sync_score_ant[IMH1]+=ant_sync(sd, cl, CMH1, time_600);             //delta 3
                    md.feat.sync_score_ant[IMMESPB]+= ant_sync(sd, cl, CMMESPB, time_600);        //delta 5
                    md.feat.amplitude_ant_time[IMH1][0]+=amp_10(sd, cl, CMH1, time_600);                //delta 2
                    md.feat.amplitude_ant_time[IMDELTA][0]+= amp_10(sd, cl, CMMESPA, time_600);         //delta 4   **using amplitude_ant_time[delta][0] to save space
                    //cout<<ant_sync(sd, cl, CMH1, time_600)<< " ";
                    //cout<<amp_10(sd, cl, CMMESPA, time_600)<< ", ";
                    //cout<<ant_sync(sd, cl, CMMESPB, time_600)<<",";
				}
                //cout<<endl;
                md.feat.amplitude_ant_time[IMDELTA][0]/=5;
				md.feat.amplitude_ant_time[IMH1][0]/=5;
                md.feat.sync_score_ant[IMH1]/=5;
				md.feat.sync_score_ant[IMMESPB]/=5;
            }
			
		}

		if (md.index==MUTANT_HER1OVER){
			if (index == 0 ){
                int time_630 = anterior_time(sd,(630/sd.step_size));
                //take only one time step since concentration may recover
                md.feat.amplitude_ant_time[IMH1][1]=amp_10(sd, cl, CMH1, time_630);
                md.feat.amplitude_ant_time[IMDELTA][1]=amp_10(sd, cl, CMDELTA, time_630);
                md.feat.amplitude_ant_time[IMMESPA][1]=amp_10(sd, cl, CMMESPA, time_630);
                //cout<<"Her1over mh1,delta,mespa amp: "<<amp_10(sd, cl, CMH1, time_630)<<" "<<amp_10(sd, cl, CMDELTA, time_630)<<" "<<amp_10(sd, cl, CMMESPA, time_630)<<endl;
            }
            
            if (index ==2){
                int time_690 = anterior_time(sd,(690/sd.step_size));
                //take only one time step since concentration may recover
                 md.feat.sync_score_ant[IMMESPB]=ant_sync(sd, cl, IMMESPB, time_690);
            }
		}

		if (md.index==MUTANT_DAPT){
			if (index == 0){
				int time_three = anterior_time(sd, (840)/sd.step_size);
                int time_three_end = anterior_time(sd, (870)/sd.step_size);
                //cout<<"DAPT mh1 sync: ";
                //cout<<"DAPT mespb sync: ";
				for (;time_three<time_three_end; time_three+=(6/sd.step_size)){
					md.feat.amplitude_ant_time[IMH1][2]+=amp_10(sd, cl, CMH1, time_three);          //DAPT 2
					md.feat.sync_score_ant[IMH1]+=ant_sync(sd, cl, CMH1, time_three);               //DAPT 1
                    md.feat.sync_score_ant[IMMESPB]+= ant_sync(sd, cl, CMMESPB, time_three);        //DAPT 4
                    //cout<<ant_sync(sd, cl, CMMESPB, time_three)<< " ";
				}
                //cout<<endl;
                md.feat.amplitude_ant_time[IMH1][2]/=5;
				md.feat.sync_score_ant[IMH1]/=5;
                md.feat.sync_score_ant[IMMESPB]/=5;
			}

			if (index==1) {
				int time_two = anterior_time(sd, (720)/sd.step_size);
				int time_two_end = anterior_time(sd, (750)/sd.step_size);
                //cout<<"DAPT immespa amp: ";
				for (;time_two<time_two_end; time_two+=(6/sd.step_size)){
					md.feat.amplitude_ant_time[IMMESPA][2]+=amp_10(sd, cl, CMMESPA, time_two);      //DAPT 3
                    //cout<<amp_10(sd, cl, CMMESPA, time_two)<< ", ";
				}
                //cout<<endl;
                md.feat.amplitude_ant_time[IMMESPA][2]/=5;
			}

        }

		if (md.index==MUTANT_MESPAOVER){                     //calculate oscillation features for mespa mutant, including posterior amplitude, anterior amplitude and syncrony score for different species
			if (index==3){
				int time_660 = anterior_time(sd, (660)/sd.step_size);         //one hour after induction, 10 snapshot in 30 minutes
                md.feat.amplitude_ant_time[IMMESPB][0]=amp_10(sd, cl, CMMESPB, time_660);   //take only one time step since concentration may recover
                //cout<<"Mespa mmespb amp: "<<amp_10(sd, cl, CMMESPB, time_660)<<endl;
			}
		}
		
		if (md.index==MUTANT_MESPBOVER){                    //calculate oscillation features for mespb mutant, including posterior amplitude, anterior amplitude and syncrony score for different species
			if (index==2 ){
                int time_660 = anterior_time(sd, (660)/sd.step_size);         //one hour after induction, 10 snapshot in 30 minutes
                md.feat.amplitude_ant_time[IMMESPB][0]=amp_10(sd, cl, CMMESPB, time_660);   //take only one time step since concentration may recover
                //cout<<"Mespb mmespb amp: "<<amp_10(sd, cl, CMMESPB, time_660)<<endl;
			}
		}


		if (ip.ant_features) {
			int time_start = anterior_time(sd, sd.steps_til_growth + (sd.width_total - sd.width_initial - 1) * sd.steps_split);
			for (int col = start_col; col < end_col; col++) {
				if (ip.ant_features) {
					plot_ant_sync(sd, cl, time_start, &features_files[SYNC], col == start_col);
				}
				time_start += sd.steps_split;
			}
			features_files[SYNC].close();
		}
	}
    
	mfree(str_set_num);
    
}

void osc_features_post (sim_data& sd, input_params& ip, con_levels& cl, features& feat, features& wtfeat, mutant_data& md, char* filename_feats, int start, int end, int set_num) {   //151221:  we are only using the this for the peaktotrough condition in wildtype mutant. maybe you can delete some unnecessary lines
	/*
	 Calculates the oscillation features: period, amplitude, and peak to trough ratio for a set of concentration levels.
	 The values are calculated using the last peak and trough of the oscillations, since the amplitude of the first few oscillations can be slightly unstable.
	 For the wild type, the peak and trough at the middle of the graph are also calculated in order to ensure that the oscillations are sustained.
	*/
    if (md.index== MUTANT_WILDTYPE || md.index== MUTANT_DELTA){
        int strlen_set_num = INT_STRLEN(set_num); // How many bytes the ASCII representation of set_num takes
        //char* str_set_num = (char*)mallocate(sizeof(char) * (strlen_set_num + 1));
        char* str_set_num = (char*) mallocate(2);
        sprintf(str_set_num, "%d", set_num);

        int con[1] = {CMH1};
        int ind[1] = {IMH1};
        static const char* concs[2] = {"mh1", "deltac"};
        static const char* feat_names[NUM_FEATURES] = {"period", "amplitude", "sync"};
        ofstream features_files[NUM_FEATURES]; // Array that will hold the files in which to output the period and amplitude

        int num_genes = 1;
        for (int i = 0; i < num_genes; i++) {
            if (ip.post_features) {
                for (int j = 0; j < NUM_FEATURES; j++) {
                    char* filename = (char*)mallocate(sizeof(char) * strlen(filename_feats) + strlen("set_") + strlen_set_num + 1 + strlen(feat_names[j]) + 1 + strlen(concs[i]) + strlen("_post.feats") + 1);
                    sprintf(filename, "%sset_%s_%s_%s_post.feats", filename_feats, str_set_num, feat_names[j], concs[i]);
                    cout << "      ";
                    open_file(&(features_files[j]), filename, false);
                    mfree(filename);
                }
                
                features_files[PERIOD] << sd.height << "," << sd.width_initial << endl;
                features_files[AMPLITUDE] << sd.height << "," << sd.width_initial << endl;
            }
        
            int mr = con[i];
            int index = ind[i];

            double period_tot = 0;
            double peaktotrough_end = 0;
            double peaktotrough_mid = 0;
            for (int x = 0; x < sd.height; x++) {
                for (int y = 0; y < sd.width_current; y++) {
                    int cell = x * sd.width_total + y;
                    growin_array peaks(sd.steps_total / (20/sd.step_size)); 
                    growin_array troughs(sd.steps_total / (20/sd.step_size));
                    int num_peaks = 0;
                    int num_troughs = 0;
                    int peaks_period = 0;

                    double cell_period = 0;
                    bool calc_period = true;
                    concentration_level<double>::timespan conc = cl.cons[mr];
                    for (int j = start+101; j < end - 1; j++) {
                        if (abs(num_peaks - num_troughs) > 1) {
                            num_peaks = 0;
                            break;
                        }

                        //check if the current point is a peak
                        if (conc[j - 1][cell] < conc[j][cell] && conc[j][cell] > conc[j + 1][cell] ) {
                            peaks[num_peaks] = j;
                            num_peaks++;
                            if (calc_period) {
                                peaks_period++;
                            }
                            
                            // add the current period to the average calculation
                            if (num_peaks >= 2 && calc_period) {
                                double period = (peaks[num_peaks - 1] - peaks[num_peaks - 2]) * sd.step_size * sd.big_gran;
                                cell_period += period;
                                if (num_peaks >= 4) {
                                    features_files[PERIOD] << period << " ";
                                }
                            }
                        }
                        
                        //check if the current point is a trough
                        if (conc[j - 1][cell] > conc[j][cell] && conc[j][cell] < conc[j + 1][cell]) {
                            troughs[num_troughs] = j;
                            num_troughs++;
                        }
                    }
                    cell_period /= (peaks_period-1);

                    if (num_peaks >= 3) {
                        int peak_penult = peaks[num_peaks - 2];
                        int trough_ult = troughs[num_peaks - 2];
                        int peak_mid = peaks[num_peaks / 2];
                        int trough_mid = troughs[num_peaks / 2];

                        period_tot += cell_period;
                        peaktotrough_end += conc[trough_ult][cell] > 1 ? conc[peak_penult][cell] / conc[trough_ult][cell] : conc[peak_penult][cell];
                        peaktotrough_mid += conc[trough_mid][cell] > 1 ? conc[peak_mid][cell] / conc[trough_mid][cell] : conc[peak_mid][cell];

                    } else {
                        period_tot += (period_tot < INFINITY ? INFINITY : 0);
                        peaktotrough_end ++;
                        peaktotrough_mid ++;
                    }
                    features_files[PERIOD] << endl;
                    features_files[AMPLITUDE] << endl;
                }
            }
            features_files[PERIOD].close();
            features_files[AMPLITUDE].close();
            int cells = sd.height * sd.width_current;
            period_tot /= cells;
            peaktotrough_end /= cells;
            peaktotrough_mid /= cells;
            if (index == IMH1) {cout<< "period post: " << period_tot<<endl;}
            feat.period_post = period_tot;
            feat.peaktotrough_end = peaktotrough_end;
            feat.peaktotrough_mid = peaktotrough_mid;
        }

        mfree(str_set_num);

    }
}


/*
 * calculate the amplitude 
 * take 10 lowest and hightest amplitude and calculate average seperately 
 * return max- min for amplitude
 */

double amp_10 (sim_data& sd, con_levels& cl, int con, int time){
    double* tosort= new double[sd.cells_total];       // allocate an array
    
    for (int i=0 ; i< sd.cells_total; i++ ){           // put cons in one timestep into this array
        tosort[i]= cl.cons[con][time][i];
    }

    sort(tosort, tosort+sd.cells_total);            // sort the array
                                                    // NOTE: sorting is done differently with resgard to NaN in Linux and Mac
    int offset=0;
    if (con== 3 || con==2){
        offset=sd.cells_total/2;                    // if we are dealing with mespa and mespb, skip first half since they are 0
    }
    
    int point_cal= 10;                              // take 10 for usual simulations; 5 for 1D simulations; 3 for mespa/mespb in 1D simulations
    if (sd.cells_total <=100){
        if (con == CMMESPA || con == CMMESPB){point_cal = 3;}
        else {point_cal =5;}
    }
    double max=0;                                   //add min and max together
    double min=0;
    for (int i=0; i< point_cal; i++){
        min+= tosort[i+offset];
        max+= tosort[sd.cells_total-i-1];
    }
    min/=point_cal;                                 //divide to get average
    max/=point_cal;
    delete[] tosort;
    return max-min;
}

double avg_amp (sim_data& sd, con_levels& cl, int con, int time, int start , int end){              //151221: calculate the average concentration value, and use it as amplitude
	int pos_start = cl.active_start_record[time];
	int pos_cur = 0;
	double conslevel = 0;

	for (int x = 0; x < sd.height; x++) {
		for (int y=start; y <end; y++){
			if (pos_start - y < 0){
				pos_cur = pos_start -2*y + sd.width_total;
			} else {
				pos_cur = pos_start -2*y;
			}
			int cell = x * sd.width_total + y + pos_cur;
			conslevel += cl.cons[con][time][cell];
		}
	}
	return conslevel / (sd.height*(end-start));
}


/*
 * calculate complementary scores for mespa and mespb in wildtype
 */
double mesp_sync (sim_data& sd, con_levels& cl, int con, int time) {  //151221: calculate complementary score of mespa and mespb
    if (sd.height == 1) {
        return 1; // for 1d arrays there is no synchronization between rows
    }
    
    double her_row[sd.width_total];
    double mesp_row[sd.width_total];
    int pos_start = cl.active_start_record[time];
    int pos_first = 0;
    int pos_cur = 0;
    double pearson_sum = 0;
    
    //iterate through each row
    for (int r = 0; r < sd.height; r++){
        //find her concentration
        for (int y = 0; y < sd.width_total; y++) {
            if (con == 3 || con ==2){
                if (pos_start - y < 0){
                    pos_first = pos_start -2*y + sd.width_total;
                } else {
                    pos_first = pos_start -2*y;
                }
            }
            int cell= r*sd.width_total + y + pos_first;
            her_row[y] = cl.cons[1][time][cell];
        }
        //find mesp concentration
        for (int x = 0; x < sd.height; x++) {
            for (int y = 0; y < sd.width_total; y++) {
                if (con == 3 || con ==2){
                    if (pos_start - y < 0){
                        pos_cur = pos_start -2*y + sd.width_total;
                    } else {
                        pos_cur = pos_start -2*y;
                    }
                }
                int cell = x * sd.width_total + y + pos_cur;
                mesp_row[y] = cl.cons[con][time][cell];
            }
            pearson_sum += pearson_correlation(her_row, mesp_row, (int)(0.6*sd.width_total), sd.width_total );   //mespa and mespb only express in anterior
        }
    }
    return (pearson_sum / (sd.height ))/ (sd.height );
}

/*
 * calculating syncronization scores for anterior part
 */
double ant_sync (sim_data& sd, con_levels& cl, int con, int time) {  //151221: calculate syncronization score
	if (sd.height == 1) {
		return 1; // for 1d arrays there is no synchronization between rows 
	}

	double first_row[sd.width_total];
	double cur_row[sd.width_total];
	int pos_start = cl.active_start_record[time];
	int pos_first = 0;
	int pos_cur = 0;
	double pearson_sum = 0;
	for (int y = 0; y < sd.width_total; y++) {
		if (con == 3 || con ==2){
			if (pos_start - y < 0){
				pos_first = pos_start -2*y + sd.width_total;
			} else {
				pos_first = pos_start -2*y;
			}
		}
        int cell= y + pos_first;
        first_row[y] = cl.cons[con][time][cell];
		//first_row[y] = cl.cons[con][time][y + pos_first];
	}


	for (int x = 1; x < sd.height; x++) {
		for (int y = 0; y < sd.width_total; y++) {
			if (con == 3 || con ==2){
				if (pos_start - y < 0){
					pos_cur = pos_start -2*y + sd.width_total;
				} else {
					pos_cur = pos_start -2*y;
				}
			}
			int cell = x * sd.width_total + y + pos_cur;
			cur_row[y] = cl.cons[con][time][cell];
		}
		if (con == 3 || con ==2){
			
			pearson_sum += pearson_correlation(first_row, cur_row, (int)(0.6*sd.width_total), sd.width_total );   //mespa and mespb only express in anterior
		} else {
			
			pearson_sum += pearson_correlation(first_row, cur_row, 0, sd.width_total);
		}
	}

	return pearson_sum / (sd.height - 1); 
}

void plot_ant_sync (sim_data& sd, con_levels& cl, int time_start, ofstream* file_pointer, bool first_col) {
	int col = cl.active_start_record[time_start];
	
	double first_row[sd.width_total * sd.steps_split];
	double other_row[sd.width_total * sd.steps_split];
	memset(first_row, 0, sizeof(double) * sd.width_total * sd.steps_split);
	memset(other_row, 0, sizeof(double) * sd.width_total * sd.steps_split);
	
	first_row[0] = cl.cons[CMH1][time_start][col];	
	int time = time_start + 1;
	for (; cl.cons[BIRTH][time][col] == cl.cons[BIRTH][time - 1][col]; time++) {
		first_row[time - time_start] = cl.cons[CMH1][time][col];
	}
	int time_end = time;
	int interval = INTERVAL / sd.step_size;
	int num_points = (time_end - time_start - interval) / (interval / 2); 
	
	if (first_col) {
		*file_pointer << sd.height - 1 << "," << INTERVAL << "," << sd.steps_split * sd.small_gran << endl;
	}

	double sync_avg[num_points];
	memset(sync_avg, 0, sizeof(double) * num_points);
	for (int x = 1; x < sd.height; x++) {
		int cell = x * sd.width_total + col;
		for (int time = time_start + 1; cl.cons[BIRTH][time][cell] == cl.cons[BIRTH][time - 1][cell]; time++) {
			other_row[time - time_start] = cl.cons[CMH1][time][cell];
		}
		
		for (int time = time_start; time <= time_end - interval; time += interval / 2) {
			sync_avg[(time - time_start) / (interval / 2)] += pearson_correlation(first_row, other_row, time - time_start, time - time_start + interval);
		}
	}
	
	for (int i = 0; i < num_points; i++) {
		sync_avg[i] /= (sd.height - 1);
		*file_pointer << sync_avg[i] << ",";
	}
	*file_pointer << endl;
}

//only used in post to measure sync_score
//compare each cell in the posterior part to the cell in middle
double post_sync (sim_data& sd, con_levels& cl, int con, int start, int end) {
	double comp_cell[end - start + 1];
	double cur_cell[end - start + 1];
	
	int middle_cell = (sd.height / 2) * sd.width_total + (sd.width_current / 2);
	
	for (int j = start; j < end; j++) {
		comp_cell[j - start] = cl.cons[con][j][middle_cell];
	}
	
	double pearson_sum = 0;
	for (int x = 0; x < sd.height; x++) {
		for (int y = 0; y < sd.width_initial; y++) {
			int cell = x * sd.width_total + y;
			
			if (cell != middle_cell) {
				for (int j = start; j < end; j++) {
					cur_cell[j - start] = cl.cons[con][j][cell];
				}
				pearson_sum += pearson_correlation(comp_cell, cur_cell, 0, end - start);
			}
		}
	}
	return pearson_sum / ((sd.height * sd.width_initial) - 1);
}

double pearson_correlation (double* x, double* y, int start, int end) {
	double x_avg = 0;
	double y_avg = 0;	
	double sigma_x2 = 0;
	double sigma_y2 = 0;
	double sigma_xy = 0;
	
	for (int j = start; j < end; j++) {
		x_avg += x[j];
		y_avg += y[j];
	}
	x_avg /= (end - start);
	y_avg /= (end - start);
	
	for (int j = start; j < end; j++) {
		sigma_xy += (x[j] - x_avg) * (y[j] - y_avg);
		sigma_x2 += SQUARE(x[j] - x_avg);
		sigma_y2 += SQUARE(y[j] - y_avg);
	}
	sigma_x2 = sqrt(sigma_x2);
	sigma_y2 = sqrt(sigma_y2);
	
	if (sigma_x2 == 0 || sigma_y2 == 0) {
		return 1;
	} else {
		return sigma_xy / ((sigma_x2 * sigma_y2));
	}
}
