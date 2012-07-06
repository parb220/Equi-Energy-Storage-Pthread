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
#include "constant.h"
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
void *initialize(void*);
void *simulation(void*); 
void TuneEnergyLevels_UpdateStorage(CEES_Pthread *); 
int main()
{
	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
	string filename_base = "../equi_energy_generic/gaussian_mixture_model."; 
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, filename_base))
	{
		cout << "Error in configuring gaussian mixture model.\n"; 
		exit (-1);
	}

	/* Initializing CEES_Pthread */ 
	CEES_Pthread::SetEnergyLevelNumber(NUMBER_ENERGY_LEVEL); 		// Number of energy levels; 
	CEES_Pthread::SetEquiEnergyJumpProb(PEE);				// Probability for equal energy jump
	CEES_Pthread::SetPeriodBuildInitialRing(BUILD_INITIAL_ENERGY_SET_PERIOD);	// Period to build initial energy ring
	CEES_Pthread::SetDataDimension(DATA_DIMENSION); 		// Data dimension for simulation
	CEES_Pthread::SetDepositFreq(DEPOSIT_FREQUENCY); 		// Frequency of deposit
	CEES_Pthread::ultimate_target = &target;	
	CEES_Pthread::SetEnergyLevels_GeometricProgression(H0, HK_1);
	CEES_Pthread::InitializeMinEnergy(); 				// For tuning energy levels based on newly identified min_energy 
	CEES_Pthread::SetTemperatures_EnergyLevels(T0, TK_1, C);
	CEES_Pthread::SetPthreadParameters(NUMBER_ENERGY_LEVEL);	// Pthread condition and mutex
 
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
	CSampleIDWeight::SetDataDimension(DATA_DIMENSION);	// Data dimension for storage (put and get) 
	CStorageHeadPthread storage(run_id, get_marker, put_marker, number_bins, string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/")); 
	storage.makedir();
	
	CEES_Pthread::storage = &storage; 
	
	/*  Generate K CEES_Pthread objects */
	CEES_Pthread *simulator = new CEES_Pthread[CEES_Pthread::GetEnergyLevelNumber()]; 
	double *sigma = new double[CEES_Pthread::GetDataDimension()];
	for (int i=0; i<CEES_Pthread::GetEnergyLevelNumber(); i++)
	{
		simulator[i].SetID_LocalTarget(i);
		simulator[i].r = r; 	// random number generator  
		if (i < CEES_Pthread::GetEnergyLevelNumber() -1)
		{
			simulator[i].SetBurnInPeriod(0);
			simulator[i].SetHigherNodePointer(simulator+i+1);
		}
		else 
		{
			simulator[i].SetHigherNodePointer(NULL);
			simulator[i].SetBurnInPeriod(BURN_IN_PERIOD);  
		}
		/* Transition_SimpleGaussian with sigma=0.25sqrt(T) used for each energy level as the proposal model */	
                for (int j=0; j<CEES_Pthread::GetDataDimension(); j++)
                        sigma[j] = 0.25 * sqrt(simulator[i].GetTemperature());
                simulator[i].SetProposal(new CTransitionModel_SimpleGaussian(CEES_Pthread::GetDataDimension(), sigma));
	}
	delete [] sigma; 

	/* Pthread */
	pthread_t *thread = new pthread_t[CEES_Pthread::GetEnergyLevelNumber()];

	/* Initializing */
	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)	
		pthread_create(&(thread[i]), NULL, initialize, (void*)(simulator+i));

	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
		pthread_join(thread[i], NULL);

	/* simulation */
	if (!IF_ENERGY_LEVEL_TUNING)
	{
		for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
                {
			simulator[i].simulationL = SIMULATION_LENGTH; 
			pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
		}
		for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
			pthread_join(thread[i], NULL);
	}
	else 
	{
		int nEnergyLevelTuning = 0; 
		while (nEnergyLevelTuning < ENERGY_LEVEL_TUNING_MAX_TIME)
		{
			for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
                	{
				simulator[i].simulationL = ENERGY_LEVEL_TUNING_FREQUENCY;
				pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
			}
                	for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
                        	pthread_join(thread[i], NULL);
			if (CEES_Pthread::IfTuneEnergyLevel())
				TuneEnergyLevels_UpdateStorage(simulator); 
			nEnergyLevelTuning ++;
		}
		if (SIMULATION_LENGTH - ENERGY_LEVEL_TUNING_MAX_TIME*ENERGY_LEVEL_TUNING_FREQUENCY > 0)
		{
			for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
                        {
                                simulator[i].simulationL = SIMULATION_LENGTH - ENERGY_LEVEL_TUNING_MAX_TIME*ENERGY_LEVEL_TUNING_FREQUENCY;
                                pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
                        }
                        for (int i=CEES_Pthread::GetEnergyLevelNumber()-1; i>=0; i--)
                                pthread_join(thread[i], NULL);
		}
	}

	
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
	summary(oFile, storage); 
	summary(oFile, simulator[CEES_Pthread::GetEnergyLevelNumber()-1]);   
	oFile.close(); 
	
	/* Release dynamically allocated space */
	delete [] thread; 
	delete [] simulator; 
	/* Release random number generator */
	gsl_rng_free(r); 
}

void *initialize(void *node_void)
{
	CEES_Pthread *simulator = (CEES_Pthread *)node_void; 
	int id = simulator->GetID(); 

	/* Wait till the next-level's initial ring is built up */
	if ( id < CEES_Pthread::GetEnergyLevelNumber()-1)
	{
		CEES_Pthread::mutex_lock(id);
		CEES_Pthread::condition_wait(id); 
		CEES_Pthread::mutex_unlock(id);
		simulator->Initialize(); 
	}
	else
	{
		/* Uniform distribution in [0, 1]^2 used for initialization. */
		double *lB = new double [CEES_Pthread::GetDataDimension()]; 
		double *uB = new double [CEES_Pthread::GetDataDimension()]; 
		for (int i=0; i<CEES_Pthread::GetDataDimension(); i++)
		{
			lB[i] = 0.0; 
			uB[i] = 1.0; 
		}
		CModel *initial_model = new CUniformModel(CEES_Pthread::GetDataDimension(), lB, uB); 
		simulator->Initialize(initial_model); 
		delete initial_model;
		delete [] lB; 
		delete [] uB; 
	}

	bool not_check_yet = true;  
	while (!simulator->EnergyRingBuildDone())
		simulator->draw();
	if (id >0 && not_check_yet)
       	{
		not_check_yet = false;
               	/* Signal the previous level to start */
                CEES_Pthread::mutex_lock(id-1);
                CEES_Pthread::condition_signal(id-1);
                CEES_Pthread::mutex_unlock(id-1);
	}
}

void *simulation(void *node_void)
{
	CEES_Pthread *simulator = (CEES_Pthread *)node_void; 
	for (int n=0; n<simulator->simulationL; n++)
	{
		if ( (IF_MH_TRACKING && n%MH_TRACKING_FREQUENCY) == 0)
			simulator->MH_Tracking_Start(MH_TRACKING_LENGTH, 0.22, 0.32); 
		simulator->draw(); 
	}
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
