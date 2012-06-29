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

void CEES_Pthread::draw()
{
	CEES_Node::draw(r, (CStorageHead &)(*storage)); 
}

void CEES_Pthread::SetHigherNodePointer(const CEES_Pthread *next)
{
	CEES_Node::SetHigherNodePointer((CEES_Node *)next);
}
