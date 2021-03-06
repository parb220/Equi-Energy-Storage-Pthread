#include <pthread.h>
#include "equi_energy_setup_constant.h"
#include "CEES_Pthread.h"
#include "CTransitionModel_SimpleGaussian.h"
#include "CUniformModel.h"

using namespace std;

void initialize_simulate(void *node_void)
{
        CEES_Pthread *simulator = (CEES_Pthread *)node_void;
        int id = simulator->GetID();
	cout << id << " ... Initializing" << endl; 
        /* Wait till the next-level's initial ring is built up */
        if ( id < CEES_Pthread::GetEnergyLevelNumber()-1)
        {
               	CEES_Pthread::mutex_lock(id);
		while (!CEES_Pthread::flag_status(id))
                	CEES_Pthread::condition_wait(id);
                CEES_Pthread::mutex_unlock(id);
		CEES_Pthread::flag_turn(id, false);
                if (!simulator->Initialize() )	// Initialize using a sample drawn from a higher energy level
		{
			CSampleIDWeight mode; 
			simulator->ultimate_target->GetMode(mode, 0); 
                        simulator->Initialize(mode); // Initialize using a sample directly drawn from target
		}
        }
        else
	{
		CSampleIDWeight mode;
                simulator->ultimate_target->GetMode(mode, 0);
                simulator->Initialize(mode); // Initialize using a sample directly drawn from target
	}

	// Set up proposal model
	int dim_cum_sum =0;  
	double *sigma;
	for (int iBlock =0; iBlock < CEES_Pthread::GetNumberBlocks(); iBlock++)
	{
		sigma = new double[CEES_Pthread::GetBlockSize(iBlock)]; 
		for (int j=0; j<CEES_Pthread::GetBlockSize(iBlock); j++)
		{
			if (id < CEES_Pthread::GetEnergyLevelNumber()-1)
				sigma[j] = simulator->GetNextLevel()->GetProposal(iBlock)->get_step_size(j); 
			else
				sigma[j] = simulator->MHProposalScale[dim_cum_sum+j]; 
		}
		dim_cum_sum += CEES_Pthread::GetBlockSize(iBlock); 
		simulator->SetProposal(new CTransitionModel_SimpleGaussian(CEES_Pthread::GetBlockSize(iBlock), sigma), iBlock); 
		delete [] sigma; 
	}
	
	cout << id << " ... Burn In" << endl; 
	simulator->BurnIn(); 
	cout << id << " ... Tune/Estimate MH Proposal StepSize" << endl; 
	simulator->MH_StepSize_Tune(); 		// Dan's adaptive strategy
	cout << id << " ... Simulate " << endl; 
	simulator->Simulate(); 
                
	//Signal the previous level to start  
	if (id > 0)
	{
        	CEES_Pthread::mutex_lock(id-1);
		CEES_Pthread::flag_turn(id-1, true);
        	CEES_Pthread::condition_signal(id-1);
       		CEES_Pthread::mutex_unlock(id-1);
	}
}

void tuning_simulation(void *node_void)
{
	CEES_Pthread *simulator = (CEES_Pthread *)node_void; 
	simulator->MH_StepSize_Tune(); 		// Dan's adaptive strategy
	simulator->Simulate(); 
}

void simulation(void *node_void)
{
        CEES_Pthread *simulator = (CEES_Pthread *)node_void;
	simulator->Simulate(); 
}
 
