#include <pthread.h>
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
bool TuneEnergyLevels_UpdateStorage(CEES_Pthread *simulator, double c_factor, double mh_target_acc)
{
	/*double new_H0_average = 0; 
	for (int i=0; i<(int)(CEES_Node::min_energy.size()); i++)
                new_H0_average += CEES_Node::min_energy[i];
        new_H0_average = new_H0_average/(int)(CEES_Node::min_energy.size());
        double new_HK_1_average = 0;
        for (int i=0; i<(int)(CEES_Node::max_energy.size()); i++)
                new_HK_1_average += CEES_Node::max_energy[i];
        new_HK_1_average=new_HK_1_average/(int)(CEES_Node::max_energy.size());*/

	double new_H0 = CEES_Node::min_energy[0] < CEES_Node::H[0] ? CEES_Node::min_energy[0] : CEES_Node::H[0]; 
	// double new_HK_1 = CEES_Node::max_energy[0] < 1.0e3 ?CEES_Node::max_energy[0]: 1.0e3;
	double new_HK_1 = CEES_Node::max_energy[0] < CEES_Node::H[CEES_Node::K-1] ? CEES_Node::max_energy[0]: CEES_Node::H[CEES_Node::K-1];

	if (new_H0 < CEES_Node::H[0] || new_HK_1 > CEES_Node::H[CEES_Node::K-1])
	{
		CEES_Pthread::SetEnergyLevels_GeometricProgression(new_H0, new_HK_1); 
		CEES_Pthread::SetTemperatures_EnergyLevels(CEES_Node::T[0], c_factor, true); 
		CEES_Pthread::SetTargetAcceptanceRate(mh_target_acc); 
	/*double new_TK_1 = new_HK_1_average * 100;
	if (new_TK_1 < CEES_Node::T[CEES_Node::K-1] *100)
        	CEES_Node::SetTemperatures_EnergyLevels(T0, new_TK_1);*/
	
	// Re-adjust local target distribution and process samples that have been generated; 
		pthread_t *thread = new pthread_t[CEES_Pthread::K]; 
		simulator->storage->CreateTemporaryBin(); 
		for (int i=0; i<CEES_Pthread::K; i++)
			pthread_create(&(thread[i]), NULL, adjust, (void *)(simulator+i)); 
	
		for (int i=0; i<CEES_Pthread::K; i++)
			pthread_join(thread[i], NULL); 
	
		simulator->storage->ClearTemporaryBin(); 
		delete [] thread; 
		return true; 
	}
	else
		return false; 
}
