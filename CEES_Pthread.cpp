#include "CEES_Pthread.h"

vector <pthread_cond_t> CEES_Pthread::condition; 
vector <pthread_mutex_t> CEES_Pthread::mutex; 
// pthread_mutex_t CEES_Pthread::global_mutex; 
vector < bool> CEES_Pthread::flag; 
CStorageHeadPthread* CEES_Pthread::storage;

CEES_Pthread::CEES_Pthread() : CEES_Node()
{
} 

CEES_Pthread::CEES_Pthread(int _id, CTransitionModel *transition, CEES_Node *next) :  CEES_Node(_id, transition, next)
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
	// pthread_mutex_init(&global_mutex, NULL); 
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
	CEES_Node::Initialize(r, initial); 
}

bool CEES_Pthread::Initialize()
{
	return CEES_Node::Initialize((CStorageHead &)(*storage), r); 
}

void CEES_Pthread::draw()
{
	CEES_Node::draw(r, (CStorageHead &)(*storage), mMH); 
}

void CEES_Pthread::BurnIn()
{
	CEES_Node::BurnIn(r, (CStorageHead &)(*storage), burnInL, mMH); 
}

void CEES_Pthread::MH_StepSize_Tune()
{
	CEES_Node::MH_StepSize_Tune(MHInitialL, MHMaxTime, r, mMH);
}

void CEES_Pthread::Simulate()
{
	CEES_Node::Simulate(r, (CStorageHead &)(*storage), simulationL, depositFreq, mMH); 
}

void CEES_Pthread::SetHigherNodePointer(const CEES_Pthread *next)
{
	CEES_Node::SetHigherNodePointer((CEES_Node *)next);
}

void CEES_Pthread::UpdateMinMaxEnergy(double _new_energy)
{
	// pthread_mutex_lock(&global_mutex); 
	CEES_Node::UpdateMinMaxEnergy(_new_energy); 
	// pthread_mutex_unlock(&global_mutex); 
}

