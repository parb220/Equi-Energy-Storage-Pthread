#include "CEES_Pthread.h"

vector <pthread_cond_t> CEES_Pthread::condition; 
vector <pthread_mutex_t> CEES_Pthread::mutex; 
CStorageHeadPthread* CEES_Pthread::storage;

CEES_Pthread::CEES_Pthread(int _id) : CEES_Node(_id)
{
} 

CEES_Pthread::CEES_Pthread(int _id, CTransitionModel *transition, CEES_Node *next) : 
CEES_Node(_id, transition, next)
{
}

CEES_Pthread::~CEES_Pthread()
{
}

void CEES_Pthread::SetPthreadParameters(int n)
{
	condition.resize(n); 
	mutex.resize(n); 
	for (int i=0; i<n; i++)
	{
		pthread_cond_init(&(condition[i]), NULL); 
		pthread_mutex_init(&(mutex[i]), NULL); 
	}
}

int CEES_Pthread::mutex_lock(int _id)
{
	return pthread_mutex_lock(&(mutex[_id])); 
}

int CEES_Pthread::mutex_unlock(int _id)
{	
	return pthread_mutex_unlock(&(mutex[_id])); 
}

int CEES_Pthread::condition_wait(int _id)
{
	return pthread_cond_wait(&(condition[_id]), &(mutex[_id])); 
}

int CEES_Pthread::condition_signal(int _id)
{
	return pthread_cond_signal(&(condition[_id])); 
}

void CEES_Pthread::Initialize(CModel *initial)
{
	CEES_Node::Initialize(initial, r); 
}

bool CEES_Pthread::Initialize()
{
	int bin_id = next_level->BinID(id); 
	int id; 
	double weight; 
	if (storage->DrawSample(bin_id, x_new, CEES_Node::GetDataDimension(), id, weight, r)) 
	{
		CEES_Node::Initialize(x_new, CEES_Node::GetDataDimension()); 
		return true; 
	}
	else 
		return false;
}

void CEES_Pthread::draw(int mMH)
{
	CEES_Node::draw(r, (CStorageHead &)(*storage), mMH); 
}

void CEES_Pthread::SetHigherNodePointer(const CEES_Pthread *next)
{
	CEES_Node::SetHigherNodePointer((CEES_Node *)next);
}

void CEES_Pthread::UpdateMinEnergy(double _new_energy)
{
	pthread_mutex_t *local_mutex = new pthread_mutex_t; 
	pthread_mutex_init(local_mutex, NULL); 
	pthread_mutex_lock(local_mutex); 
	if (_new_energy < min_energy)
        {
                min_energy = _new_energy;
                if (min_energy < H[0])
                        if_tune_energy_level = true;
        }
	pthread_mutex_unlock(local_mutex); 
	pthread_mutex_destroy(local_mutex); 
}

void CEES_Pthread::AssignSamplesGeneratedSoFar()
{
	CEES_Node::AssignSamplesGeneratedSoFar((CStorageHead &)(*storage)); 
}
