#include "CStorageHeadPthread.h"

CStorageHeadPthread::CStorageHeadPthread(int _run_id, int _get_marker, int _put_marker, int _number_bins, string file) : CStorageHead(_run_id, _get_marker, _put_marker, _number_bins, file)
{
	mutex.resize(_number_bins); 
	for (int i=0; i<(int)(mutex.size()); i++)
		pthread_mutex_init(&(mutex[i]), NULL); 
}

CStorageHeadPthread::~CStorageHeadPthread() 
{
}

int CStorageHeadPthread::DepositSample(int _bin_index, const CSampleIDWeight &sample)
{
	pthread_mutex_lock(&(mutex[_bin_index]));
	int return_value = CStorageHead::DepositSample(_bin_index, sample); 
	pthread_mutex_unlock(&(mutex[_bin_index])); 
	return return_value;
}

int CStorageHeadPthread::DepositSample(int _bin_index, const double *_x, int _dim, int _id, double _weight)
{
	pthread_mutex_lock(&(mutex[_bin_index])); 
	int return_value = CStorageHead::DepositSample(_bin_index, _x, _dim, _id, _weight); 
	pthread_mutex_unlock(&(mutex[_bin_index])); 
	return return_value; 
}

CSampleIDWeight CStorageHeadPthread::DrawSample(int _bin_index, const gsl_rng *r)
{
	pthread_mutex_lock(&(mutex[_bin_index])); 
	CSampleIDWeight sample=CStorageHead::DrawSample(_bin_index, r); 
	pthread_mutex_unlock(&(mutex[_bin_index])); 
	return sample;
}

void CStorageHeadPthread::DrawSample(int _bin_index, double *_x, int _dim, int & _id, double & _weight, const gsl_rng *r)
{
	pthread_mutex_lock(&(mutex[_bin_index])); 
	CStorageHead::DrawSample(_bin_index, _x, _dim, _id, _weight, r); 
	pthread_mutex_unlock(&(mutex[_bin_index])); 
}

