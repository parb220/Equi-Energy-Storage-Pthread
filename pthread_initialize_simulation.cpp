#include <pthread.h>
#include "equi_energy_setup_constant.h"
#include "CEES_Pthread.h"
#include "CUniformModel.h"

using namespace std;

void *initialize_simulate(void *node_void)
{
        CEES_Pthread *simulator = (CEES_Pthread *)node_void;
        int id = simulator->GetID();
        /* Uniform distribution in [0, 1]^d used for initialization. */
        double *lB = new double [CEES_Pthread::GetDataDimension()];
        double *uB = new double [CEES_Pthread::GetDataDimension()];
        for (int i=0; i<CEES_Pthread::GetDataDimension(); i++)
        {
                lB[i] = 0.0;
                uB[i] = 1.0;
        }
        CModel *initial_model = new CUniformModel(CEES_Pthread::GetDataDimension(), lB, uB);

        /* Wait till the next-level's initial ring is built up */
        if ( id < CEES_Pthread::GetEnergyLevelNumber()-1)
        {
                CEES_Pthread::mutex_lock(id);
                CEES_Pthread::condition_wait(id);
                CEES_Pthread::mutex_unlock(id);
                if (!simulator->Initialize() )
                        simulator->Initialize(initial_model);
        }
        else
                simulator->Initialize(initial_model);

        delete initial_model;
        delete [] lB;
        delete [] uB;

        /* If burning-in and initial ring built-up is finished, then signal to wake up relevant threads. In this period, MH tuning is allowed. */
        bool not_check_yet = true;
        int n=0;
        while (!simulator->EnergyRingBuildDone())
        {
                if ( (IF_MH_TRACKING && n%MH_TRACKING_FREQUENCY) == 0)
                        simulator->MH_Tracking_Start(MH_TRACKING_LENGTH, 0.1, 0.6);

                simulator->draw(MULTIPLE_TRY_MH);
                n++;
        }
        if (id >0 && not_check_yet)
        {
                not_check_yet = false;
                /* Signal the previous level to start */
                CEES_Pthread::mutex_lock(id-1);
                CEES_Pthread::condition_signal(id-1);
                CEES_Pthread::mutex_unlock(id-1);

		cout << "ring " << id << ": initial ring built done.\n"; 
        }

	for (int n=0; n<simulator->simulationL; n++)
        {
                if ( (IF_MH_TRACKING && n%MH_TRACKING_FREQUENCY) == 0)
                        simulator->MH_Tracking_Start(MH_TRACKING_LENGTH, 0.1, 0.6);
                simulator->draw(MULTIPLE_TRY_MH);
        }
}

void *simulation(void *node_void)
{
        CEES_Pthread *simulator = (CEES_Pthread *)node_void;
        for (int n=0; n<simulator->simulationL; n++)
        {
                if ( (IF_MH_TRACKING && n%MH_TRACKING_FREQUENCY) == 0)
                        simulator->MH_Tracking_Start(MH_TRACKING_LENGTH, 0.1, 0.6);
                simulator->draw(MULTIPLE_TRY_MH);
        }
}
 
