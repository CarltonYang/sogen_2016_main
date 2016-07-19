/*
 *  randnum.cpp
 *  
 *
 *  Created by Adam Lock on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;

int main(){
	
	ofstream ran;
	ran.open("randnums.params", fstream::out);
	/* Names of Random Numbers, and dXXX = delayXXX so dmh1 is really delaymh1 */
	
	/* Set the first line to Numbers to values that we know actually work */
	
	/* Create the variables for the multiple randoms */
    double RMSH1, RMSMESPA, RMSMESPB, RMSDELTA; //0-3
    double RMDH1, RMDMESPA, RMDMESPB, RMDDELTA; //4-7
	double RPSH1, RPSMESPA, RPSMESPB, RPSDELTA; //8-11
	double RPDH1, RPDMESPA, RPDMESPB, RPDDELTA; //12-15
    
	double RDAH1H1, RDAMESPAMESPA, RDAMESPAMESPB, RDAMESPBMESPB; //16-19
	double RDDIH1H1, RDDIMESPAMESPA, RDDIMESPAMESPB, RDDIMESPBMESPB; //20-23
	double RDDGH1H1, RDDGMESPAMESPA, RDDGMESPAMESPB, RDDGMESPBMESPB; //24-27
	double RDELAYMH1, RDELAYMMESPA, RDELAYMMESPB, RDELAYMDELTA; //28-31
	double RDELAYPH1, RDELAYPMESPA, RDELAYPMESPB, RDELAYPDELTA; //32-35
	double RCRITPH1H1, RCRITPDELTA, RCRITPMESPAMESPA, RCRITPMESPAMESPB, RCRITPMESPBMESPB; //36-40
	double NS1, NS2;
	
	/* Number of Random Numbers */
	int count = 200;
	
	srand((unsigned)time(0));
	double RMSH1a = 30, RMSH1b = 65;
	double RMSMESPAa = 30, RMSMESPAb = 65;
	double RMSMESPBa = 30, RMSMESPBb = 65;
	double RMSDELTAa = 30, RMSDELTAb = 65;
    
	double RMDH1a = 0.1, RMDH1b = 0.46;
	double RMDMESPAa = 0.1, RMDMESPAb = 0.46;
	double RMDMESPBa = 0.1, RMDMESPBb = 0.46;
	double RMDDELTAa = 0.1, RMDDELTAb = 0.46;
    
	double RPSH1a = 10, RPSH1b = 60;
	double RPSMESPAa = 10, RPSMESPAb = 60;
	double RPSMESPBa = 10, RPSMESPBb = 60;
	double RPSDELTAa = 20, RPSDELTAb = 60;
    
    double RPDH1a = 0.11, RPDH1b = 0.35;
    double RPDMESPAa = 0.11, RPDMESPAb = 0.35;
    double RPDMESPBa = 0.11, RPDMESPBb = 0.35;
    double RPDDELTAa = 0.07, RPDDELTAb = 0.35;
    
	double RDAH1H1a = 0.0003, RDAH1H1b = 0.03;
	double RDAMESPAMESPAa = 0.0003, RDAMESPAMESPAb = 0.03;
	double RDAMESPAMESPBa = 0.0003, RDAMESPAMESPBb = 0.03;
	double RDAMESPBMESPBa = 0.0003, RDAMESPBMESPBb = 0.03;
    
	double RDDIH1H1a = 0.0003, RDDIH1H1b = 0.03;
	double RDDIMESPAMESPAa = 0.0003, RDDIMESPAMESPAb = 0.03;
	double RDDIMESPAMESPBa = 0.0003, RDDIMESPAMESPBb = 0.03;
	double RDDIMESPBMESPBa = 0.0003, RDDIMESPBMESPBb = 0.03;
    
	double RDDGH1H1a = 0.11, RDDGH1H1b = 0.35;
	double RDDGMESPAMESPAa = 0.11, RDDGMESPAMESPAb = 0.35;
    double RDDGMESPAMESPBa = 0.11, RDDGMESPAMESPBb = 0.35;
    double RDDGMESPBMESPBa = 0.11, RDDGMESPBMESPBb = 0.35;
    
	double RDELAYMH1a = 8.0, RDELAYMH1b = 12.0;
	double RDELAYMMESPAa = 8.0, RDELAYMMESPAb = 12.0;
	double RDELAYMMESPBa = 8.0, RDELAYMMESPBb = 12.0;
	double RDELAYMDELTAa = 8.0, RDELAYMDELTAb = 12.0;
    
	double RDELAYPH1a = 0.4, RDELAYPH1b = 2.0;
	double RDELAYPMESPAa = 0.4, RDELAYPMESPAb = 2.0;
	double RDELAYPMESPBa = 0.4, RDELAYPMESPBb = 2.0;
	double RDELAYPDELTAa = 10.0, RDELAYPDELTAb = 25;
    
	double RCRITPH1H1a = 150.0, RCRITPH1H1b = 900.0;
	double RCRITPDELTAa = 200.0, RCRITPDELTAb = 800.0;
	double RCRITPMESPAMESPAa = 800.0, RCRITPMESPAMESPAb = 2400.0;
	double RCRITPMESPAMESPBa = 800.0, RCRITPMESPAMESPBb = 2400.0;
	double RCRITPMESPBMESPBa = 800.0, RCRITPMESPBMESPBb = 2400.0;

    double NS1a = 0.0, NS1b = 1.0;
	double NS2a = 0.0, NS2b = 1.0;
	
	while(count!=0){
		// Random doubles for every object
        
        RMSH1= double(((RMSH1b - RMSH1a)*rand()/(RAND_MAX + 1.0))+RMSH1a);
        RMSMESPA= double(((RMSMESPAb - RMSMESPAa)*rand()/(RAND_MAX + 1.0))+RMSMESPAa);
        RMSMESPB= double(((RMSMESPBb - RMSMESPBa)*rand()/(RAND_MAX + 1.0))+RMSMESPBa);
        RMSDELTA= double(((RMSDELTAb - RMSDELTAa)*rand()/(RAND_MAX + 1.0))+RMSDELTAa); //0-3
        
        
        RMDH1= double(((RMDH1b - RMDH1a)*rand()/(RAND_MAX + 1.0))+RMDH1a);
        RMDMESPA= double(((RMDMESPAb - RMDMESPAa)*rand()/(RAND_MAX + 1.0))+RMDMESPAa);
        RMDMESPB= double(((RMDMESPBb - RMDMESPBa)*rand()/(RAND_MAX + 1.0))+RMDMESPBa);
        RMDDELTA= double(((RMDDELTAb - RMDDELTAa)*rand()/(RAND_MAX + 1.0))+RMDDELTAa); //4-7
        
        RPSH1= double(((RPSH1b - RPSH1a)*rand()/(RAND_MAX + 1.0))+RPSH1a);
        RPSMESPA= double(((RPSMESPAb - RPSMESPAa)*rand()/(RAND_MAX + 1.0))+RPSMESPAa);
        RPSMESPB= double(((RPSMESPBb - RPSMESPBa)*rand()/(RAND_MAX + 1.0))+RPSMESPBa);
        RPSDELTA= double(((RPSDELTAb - RPSDELTAa)*rand()/(RAND_MAX + 1.0))+RPSDELTAa); //8-11
        
        RPDH1= double(((RPDH1b - RPDH1a)*rand()/(RAND_MAX + 1.0))+RPDH1a);
        RPDMESPA= double(((RPDMESPAb - RPDMESPAa)*rand()/(RAND_MAX + 1.0))+RPDMESPAa);
        RPDMESPB= double(((RPDMESPBb - RPDMESPBa)*rand()/(RAND_MAX + 1.0))+RPDMESPBa);
        RPDDELTA= double(((RPDDELTAb - RPDDELTAa)*rand()/(RAND_MAX + 1.0))+RPDDELTAa); //12-15
        
        RDAH1H1= double(((RDAH1H1b - RDAH1H1a)*rand()/(RAND_MAX + 1.0))+RDAH1H1a);
        RDAMESPAMESPA= double(((RDAMESPAMESPAb - RDAMESPAMESPAa)*rand()/(RAND_MAX + 1.0))+RDAMESPAMESPAa);
        RDAMESPAMESPB= double(((RDAMESPAMESPBb - RDAMESPAMESPBa)*rand()/(RAND_MAX + 1.0))+RDAMESPAMESPBa);
        RDAMESPBMESPB= double(((RDAMESPBMESPBb - RDAMESPBMESPBa)*rand()/(RAND_MAX + 1.0))+RDAMESPBMESPBa); // 16-19
        
        RDDIH1H1= double(((RDDIH1H1b - RDDIH1H1a)*rand()/(RAND_MAX + 1.0))+RDDIH1H1a);
        RDDIMESPAMESPA= double(((RDDIMESPAMESPAb - RDDIMESPAMESPAa)*rand()/(RAND_MAX + 1.0))+RDDIMESPAMESPAa);
        RDDIMESPAMESPB= double(((RDDIMESPAMESPBb - RDDIMESPAMESPBa)*rand()/(RAND_MAX + 1.0))+RDDIMESPAMESPBa);
        RDDIMESPBMESPB= double(((RDDIMESPBMESPBb - RDDIMESPBMESPBa)*rand()/(RAND_MAX + 1.0))+RDDIMESPBMESPBa); //20-23
        
        RDDGH1H1= double(((RDDGH1H1b - RDDGH1H1a)*rand()/(RAND_MAX + 1.0))+RDDGH1H1a);
        RDDGMESPAMESPA= double(((RDDGMESPAMESPAb - RDDGMESPAMESPAa)*rand()/(RAND_MAX + 1.0))+RDDGMESPAMESPAa);
        RDDGMESPAMESPB= double(((RDDGMESPAMESPBb - RDDGMESPAMESPBa)*rand()/(RAND_MAX + 1.0))+RDDGMESPAMESPBa);
        RDDGMESPBMESPB= double(((RDDGMESPBMESPBb - RDDGMESPBMESPBa)*rand()/(RAND_MAX + 1.0))+RDDGMESPBMESPBa); //24-27
        
        RDELAYMH1= double(((RDELAYMH1b - RDELAYMH1a)*rand()/(RAND_MAX + 1.0))+RDELAYMH1a);
        RDELAYMMESPA= double(((RDELAYMMESPAb - RDELAYMMESPAa)*rand()/(RAND_MAX + 1.0))+RDELAYMMESPAa);
        RDELAYMMESPB= double(((RDELAYMMESPBb - RDELAYMMESPBa)*rand()/(RAND_MAX + 1.0))+RDELAYMMESPBa);
        RDELAYMDELTA= double(((RDELAYMDELTAb - RDELAYMDELTAa)*rand()/(RAND_MAX + 1.0))+RDELAYMDELTAa); //28-31
        
        RDELAYPH1= double(((RDELAYPH1b - RDELAYPH1a)*rand()/(RAND_MAX + 1.0))+RDELAYPH1a);
        RDELAYPMESPA= double(((RDELAYPMESPAb - RDELAYPMESPAa)*rand()/(RAND_MAX + 1.0))+RDELAYPMESPAa);
        RDELAYPMESPB= double(((RDELAYPMESPBb - RDELAYPMESPBa)*rand()/(RAND_MAX + 1.0))+RDELAYPMESPBa);
        RDELAYPDELTA= double(((RDELAYPDELTAb - RDELAYPDELTAa)*rand()/(RAND_MAX + 1.0))+RDELAYPDELTAa); //32-35
        
        RCRITPH1H1= double(((RCRITPH1H1b - RCRITPH1H1a)*rand()/(RAND_MAX + 1.0))+RCRITPH1H1a);
        RCRITPDELTA= double(((RCRITPDELTAb - RCRITPDELTAa)*rand()/(RAND_MAX + 1.0))+RCRITPDELTAa);
        RCRITPMESPAMESPA= double(((RCRITPMESPAMESPAb - RCRITPMESPAMESPAa)*rand()/(RAND_MAX + 1.0))+RCRITPMESPAMESPAa);
        RCRITPMESPAMESPB= double(((RCRITPMESPAMESPBb - RCRITPMESPAMESPBa)*rand()/(RAND_MAX + 1.0))+RCRITPMESPAMESPBa);
        RCRITPMESPBMESPB= double(((RCRITPMESPBMESPBb - RCRITPMESPBMESPBa)*rand()/(RAND_MAX + 1.0))+RCRITPMESPBMESPBa); //36-40
        
        NS1= double(((NS1b - NS1a)*rand()/(RAND_MAX + 1.0))+NS1a);
        NS2= double(((NS2b - NS2a)*rand()/(RAND_MAX + 1.0))+NS2a);
        
		
		/* Place them all in the file separated by commas */
		ran<<RMSH1<<","<<RMSMESPA<<","<<RMSMESPB<<","<<RMSDELTA<<","<<RMDH1<<","<<RMDMESPA<<","<<RMDMESPB<<","<<RMDDELTA<<","<<RPSH1<<","<<RPSMESPA<<","<<RPSMESPB<<","<<RPSDELTA<<",";
		ran<<RPDH1<<","<<RPDMESPA<<","<<RPDMESPB<<","<<RPDDELTA<<","<<RDAH1H1<<","<<RDAMESPAMESPA<<","<<RDAMESPAMESPB<<","<<RDAMESPBMESPB<<","<<RDDGH1H1<<","<<RDDIMESPAMESPA<<","<<RDDIMESPAMESPB<<","<<RDDIMESPBMESPB<<",";
		ran<<RDDGH1H1<<","<<RDDGMESPAMESPA<<","<<RDDGMESPAMESPB<<","<<RDDGMESPBMESPB<<","<<RDELAYMH1<<","<<RDELAYMMESPA<<","<<RDELAYMMESPB<<","<<RDELAYMDELTA<<","<<RDELAYPH1<<","<<RDELAYPMESPA<<","<<RDELAYPMESPB<<","<<RDELAYPDELTA<<",";
		ran<<RCRITPH1H1<<","<<RCRITPDELTA<<","<<RCRITPMESPAMESPA<<","<<RCRITPMESPAMESPB<<","<<RCRITPMESPBMESPB<<","<<NS1<<","<<NS2<<endl;
		
		count--;
		
		if (count % 200 == 0)
			cout << count << " numbers left" << endl;
	}
	ran.close(); /* Close the file when done */
}
