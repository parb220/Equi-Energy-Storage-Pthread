/* Test: 
 * CEquiEnergy
 * CModel
 * CMixtureModel
 * CTransitionModel
 * CSimpleGaussianModel
 * CTransitionModel_SimpleGaussian
 * CUniformModel
 * CBoundedModel
 * CEES_Head
 * CEES_Node
 * CStorageHead
 * CPutGetBin
 * CSampleIDWeight
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include "equi_energy_setup_constant.h"
#include "CMixtureModel.h"
#include "CSimpleGaussianModel.h"
#include "CEquiEnergy.h"
#include "CUniformModel.h"
#include "CBoundedModel.h"
#include "CTransitionModel.h"
#include "CTransitionModel_SimpleGaussian.h"
#include "CEES_Pthread.h"
#include "CSampleIDWeight.h"
#include "CStorageHead.h"
#include "CStorageHeadPthread.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 
void *initialize_simulate(void*);
void *simulation(void*); 
void *tuning_simulation(void *); 
bool TuneEnergyLevels_UpdateStorage(CEES_Pthread *, double, double); 

void usage(int arc, char **argv)
{
        cout << "usage: " << argv[0] << " ";
        cout << "-d <dimension> \n";
        cout << "-f <files of the target model> \n";
        cout << "-p <probability of equi-energy jump> \n";
        cout << "-h <energy bound of the highest energy level>\n";
        cout << "-l <simulation length>\n";
	cout << "-c <C factor to determine temperature bounds according to energy bounds>\n"; 
	cout << "? this message\n";
}

int main(int argc, char ** argv)
{
	string filename_base = "../equi_energy_generic/gaussian_mixture_model.";
        int data_dimension = DATA_DIMENSION;
        double pee = PEE;
        double h_k_1 = HK_1;
        int simulation_length =SIMULATION_LENGTH;
	double c_factor = C; 
	double mh_target_acc = MH_TARGET_ACC; 

	int opt;
        while ( (opt = getopt(argc, argv, "d:f:p:h:l:c:?")) != -1)
        {
                switch (opt)
                {
                        case 'f':
                                filename_base = string(optarg); break;
                        case 'd':
                                data_dimension = atoi(optarg); break;
                        case 'p':
                                pee = atof(optarg); break;
                        case 'h':
                                h_k_1 = atof(optarg); break;
                        case 'l':
                                simulation_length = atoi(optarg); break;
			case 'c': 
				c_factor = atof(optarg); break; 
			case '?':
			{
				usage(argc, argv); 
				exit(-1); 
			}
                        default:
                        {
                                usage(argc, argv);
                                exit(-1);
                        }
                }
        }


	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, filename_base))
	{
		cout << "Error in configuring gaussian mixture model.\n"; 
		exit (-1);
	}

	/* Initializing CEES_Pthread */ 
	CEES_Pthread::SetEnergyLevelNumber(NUMBER_ENERGY_LEVEL); 	// Number of energy levels; 
	CEES_Pthread::SetEquiEnergyJumpProb(pee);			// Probability for equal energy jump
	CEES_Pthread::SetDataDimension(data_dimension); 		// Data dimension for simulation
	CEES_Pthread::ultimate_target = &target;	
	
	CEES_Pthread::SetEnergyLevels_GeometricProgression(H0, h_k_1);
	CEES_Pthread::InitializeMinMaxEnergy(ENERGY_TRACKING_NUMBER);	// For tuning energy levels based on newly identified min_energy 
	CEES_Pthread::SetTemperatures_EnergyLevels(T0, c_factor, true); 
	CEES_Pthread::SetTargetAcceptanceRate(mh_target_acc); 
	CEES_Pthread::SetPthreadParameters(NUMBER_ENERGY_LEVEL);	// Pthread condition and mutex
	
	if (MH_BLOCK)							// MH in blocks
		CEES_Pthread::SetBlockSize(NULL, CEES_Pthread::GetDataDimension()); 
	else 
		CEES_Pthread::SetBlockSize(NULL); 
	/*
 	Initialize the random_number_generator which will be used for all distributions to draw samples.
 	*/
	const gsl_rng_type *T; 
	gsl_rng *r;
	gsl_rng_env_setup(); 
	T = gsl_rng_default; 
	r = gsl_rng_alloc(T); 
	gsl_rng_set(r, (unsigned)time(NULL)); 	
	
	/* Initialize Storage  */
	int run_id = time(NULL); // by default, use current time as run_id; 
	int get_marker = 10000;	// Keep get_marker samples in memory for draw 
	int put_marker = 10000;	// Dump every put_marker samples
	int number_bins = NUMBER_ENERGY_LEVEL * NUMBER_ENERGY_LEVEL;   
	CSampleIDWeight::SetDataDimension(data_dimension);	// Data dimension for storage (put and get) 
	CStorageHeadPthread storage(run_id, get_marker, put_marker, number_bins, string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/")); 
	storage.makedir();
	
	CEES_Pthread::storage = &storage; 
	
	/*  Generate K CEES_Pthread objects */
	CEES_Pthread *simulator = new CEES_Pthread[CEES_Pthread::GetEnergyLevelNumber()]; 
	for (int i=0; i<CEES_Pthread::GetEnergyLevelNumber(); i++)
	{
		simulator[i].SetID_LocalTarget(i);
		simulator[i].r = r; 	// random number generator  
		if (i < CEES_Pthread::GetEnergyLevelNumber() -1)
			simulator[i].SetHigherNodePointer(simulator+i+1);
		else 
			simulator[i].SetHigherNodePointer(NULL);
	}

	/* Pthread */
	pthread_t *thread = new pthread_t[CEES_Pthread::GetEnergyLevelNumber()];

	/* Initializing - burn-in - MH-stepsize-Regression -- simulate*/
	cout << "Initialize, burn in, tune/estimate MH stepsize and simulate for " << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << " steps.\n"; 
	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
	{
		simulator[i].burnInL = BURN_IN_PERIOD; 
		simulator[i].mMH = MULTIPLE_TRY_MH; 
		//simulator[i].MHMaxTime = MH_STEPSIZE_TUNING_MAX_TIME; 
		//simulator[i].MHInitialL = MH_TRACKING_LENGTH; 
		simulator[i].MHInitialL = 20; 
		simulator[i].MHMaxTime = 10; 
		simulator[i].simulationL = ENERGY_LEVEL_TRACKING_WINDOW_LENGTH;
		simulator[i].depositFreq = DEPOSIT_FREQUENCY; 
		pthread_create(&(thread[i]), NULL, initialize_simulate, (void*)(simulator+i));
	}
	
	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
		pthread_join(thread[i], NULL);

	int nEnergyLevelTuning = 0;
	while (nEnergyLevelTuning < ENERGY_LEVEL_TUNING_MAX_TIME)
        {       
		cout << "Energy level tuning: " << nEnergyLevelTuning << " for " << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << " steps.\n"; 
		TuneEnergyLevels_UpdateStorage(simulator, c_factor, mh_target_acc);
		for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
		{
			simulator[i].simulationL = ENERGY_LEVEL_TRACKING_WINDOW_LENGTH; 
                       	pthread_create(&(thread[i]), NULL, tuning_simulation, (void*)(simulator+i));       
               	}
               	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)			
			pthread_join(thread[i], NULL);
		nEnergyLevelTuning ++;
	}
	// run through simulation
	cout << "Simulation for " << simulation_length << " steps.\n"; 
	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
        {
		simulator[i].simulationL = simulation_length ;
               	pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
        }
        for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
               	pthread_join(thread[i], NULL);
	
	storage.finalize(); 		// save to hard-disk of those unsaved data

	/* finalize: write CStorageHead and CEES_Pthread information into a file.
 	*/
	string file = storage.GetSummaryFileName(); 
	ofstream oFile; 
	oFile.open(file.c_str());
	if (!oFile)
	{
		cout << "Error in writing the summary file.\n"; 
		exit(-1); 
	}
	oFile << "Burn In:\t" << BURN_IN_PERIOD << endl;
        oFile << "Tune Energy Level Window Length:\t" << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << endl;
        oFile << "Tune Energy Level Number:\t" << ENERGY_LEVEL_TUNING_MAX_TIME << endl;
        oFile << "Deposit Frequency:\t" << DEPOSIT_FREQUENCY << endl;
        oFile << "MH Target Probability:\t" << mh_target_acc << endl;
        oFile << "MH Initial Window Length:\t" << MH_TRACKING_LENGTH << endl;
        oFile << "MH Window Number:\t" << MH_STEPSIZE_TUNING_MAX_TIME << endl;
        summary(oFile, storage);
        summary(oFile, simulator);
        for (int i=0; i<CEES_Node::GetEnergyLevelNumber(); i++)
                summary(oFile, simulator[i], i);

	oFile.close(); 
	
	/* Release dynamically allocated space */
	delete [] thread; 
	delete [] simulator; 
	/* Release random number generator */
	gsl_rng_free(r); 
}

bool Configure_GaussianMixtureModel_File(CMixtureModel &mixture_model, const string filename_base)
{
	/*weight */
	string filename = filename_base + "weight";
	int nComponent, nDim; 			// Number of components, dimension of variables
	bool equalComponent, equalDim;		// Whether parameters for different components are thes same; and whether parameters for different dimension are the same; 
	ifstream inputFile; 
	inputFile.open(filename.data());
	if (!inputFile)
		return false;
	inputFile >> nComponent; 
	double *weight = new double[nComponent]; 
	inputFile >> equalComponent; 
	if (!equalComponent)
	{
		for (int i=0; i<nComponent; i++)
			inputFile >> weight[i]; 
	} 
	else 
	{
		inputFile >> weight[0]; 
		for (int i=1; i<nComponent; i++)
			weight[i] = weight[0];
	}
	mixture_model.SetModelNumber(nComponent); 
	mixture_model.SetWeightParameter(weight, nComponent); 
	delete [] weight; 
	inputFile.close();

	/*sigma */
	filename = filename_base + "sigma"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** sigma = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		sigma[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> sigma[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> sigma[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				sigma[i][j] = sigma[i][0]; 
		} 
	}
	else
	{
		inputFile >> sigma[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][0];
		}
	}

	inputFile.close(); 

	/*mean */
	filename = filename_base + "mean"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** mean = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		mean[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> mean[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> mean[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				mean[i][j] = mean[i][0]; 
		} 
	}
	else
	{
		inputFile >> mean[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][0];
		}
	}

	inputFile.close(); 

	mixture_model.SetDataDimension(nDim); 
	for (int i=0; i<nComponent; i++)
	{
		mixture_model.Initialize(i, new CSimpleGaussianModel(nDim, mean[i], sigma[i])); 
		delete []mean[i]; 
		delete []sigma[i]; 
	}
	mixture_model.CalculateSetParameterNumber();
	return true;
}
