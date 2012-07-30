#include "CEES_Pthread.h"

vector <pthread_cond_t> CEES_Pthread::condition; 
vector <pthread_mutex_t> CEES_Pthread::mutex; 
pthread_mutex_t CEES_Pthread::global_mutex; 
vector < bool> CEES_Pthread::flag; 
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
	flag.resize(n); 
	for (int i=0; i<n; i++)
	{
		pthread_cond_init(&(condition[i]), NULL); 
		pthread_mutex_init(&(mutex[i]), NULL); 
		flag[i] = false; 
	}
	pthread_mutex_init(&global_mutex, NULL); 
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
	int sample_id; 
	double sample_weight; 
	int level_id=0; 
	int bin_id = next_level->BinID(level_id); 
	while(level_id < K && ! storage->DrawSample(bin_id, x_new, GetDataDimension(), sample_id, sample_weight, r)) 
	{
		level_id ++; 
		bin_id = next_level->BinID(level_id); 
	}
	if (level_id < K)
	{
		CEES_Node::Initialize(x_new, GetDataDimension()); 
		return true; 
	}
	else 
		return false;
}

void CEES_Pthread::draw(int mMH)
{
	CEES_Node::draw(r, (CStorageHead &)(*storage), mMH); 
}

void CEES_Pthread::draw_block()
{
	CEES_Node::draw_block(r, (CStorageHead &)(*storage)); 
}

void CEES_Pthread::SetHigherNodePointer(const CEES_Pthread *next)
{
	CEES_Node::SetHigherNodePointer((CEES_Node *)next);
}

void CEES_Pthread::UpdateMinMaxEnergy(double _new_energy)
{
	pthread_mutex_lock(&global_mutex); 
	CEES_Node::UpdateMinMaxEnergy(_new_energy); 
	pthread_mutex_unlock(&global_mutex); 
}

void CEES_Pthread::AssignSamplesGeneratedSoFar()
{
	CEES_Node::AssignSamplesGeneratedSoFar((CStorageHead &)(*storage)); 
}
