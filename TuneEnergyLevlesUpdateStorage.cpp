#include <pthread.h>
#include "equi_energy_setup_constant.h"
#include "CEES_Pthread.h"
#include "CStorageHeadPthread.h"
#include "CBoundedModel.h"

using namespace std;
void *adjust(void *node_void)
{
	CEES_Pthread *simulator = (CEES_Pthread *) node_void; 
	simulator->AdjustLocalTarget(); 
	simulator->AssignSamplesGeneratedSoFar();
}
void TuneEnergyLevels_UpdateStorage(CEES_Pthread *simulator)
{
	CEES_Pthread::if_tune_energy_level = false; 
	double new_H0 = CEES_Pthread::min_energy; 

	// Re-determine and adjust energy level and temperature levels
	if (!CEES_Pthread::SetEnergyLevels_GeometricProgression(new_H0, HK_1))
        {
                cout << "Error in setting energy levels." << endl;
                exit(-1);
        }

	if (!CEES_Pthread::SetTemperatures_EnergyLevels(T0, TK_1, C) )
        {
                cout << "Error in setting temperature levels." << endl;
                exit(-1);
        }
	
	// Re-adjust local target distribution and process samples that have been generated; 
	pthread_t *thread = new pthread_t[CEES_Pthread::K]; 
	simulator->storage->CreateTemporaryBin(); 
	for (int i=0; i<CEES_Pthread::K; i++)
		pthread_create(&(thread[i]), NULL, adjust, (void *)(simulator+i)); 
	
	for (int i=0; i<CEES_Pthread::K; i++)
		pthread_join(thread[i], NULL); 
	
	simulator->storage->ClearTemporaryBin(); 
}
