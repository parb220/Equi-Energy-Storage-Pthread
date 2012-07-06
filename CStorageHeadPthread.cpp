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

vector <CSampleIDWeight> CStorageHeadPthread::RetrieveSamplesSequentially(bool if_clear_old_bin, int bin_id)
{	
	pthread_mutex_lock(&(mutex[bin_id])); 
	vector <CSampleIDWeight> samples = CStorageHead::RetrieveSamplesSequentially(if_clear_old_bin, bin_id);
	pthread_mutex_unlock(&(mutex[bin_id])); 
	return samples; 	
}

void CStorageHeadPthread::CreateTemporaryBin()
{
	CStorageHead::CreateTemporaryBin(); 
	if ((int)(mutex.size()) < 2*number_bins)
	{
		int old_size = (int)(mutex.size()); 
		mutex.resize(2*number_bins);
        	for (int i=old_size; i<(int)(mutex.size()); i++)
                	pthread_mutex_init(&(mutex[i]), NULL);
	}
}

int CStorageHeadPthread::DepositSample(bool if_new_bin, int bin_id, const double *x, int dX, int _id, double _weight)
{

	int return_value; 
	if (if_new_bin)
	{
		pthread_mutex_lock(&(mutex[bin_id+number_bins])); 
		return_value = bin[bin_id+number_bins].DepositSample(x, dX, _id, _weight); 
		pthread_mutex_unlock(&(mutex[bin_id+number_bins])); 
	}
	else 
	{
		pthread_mutex_lock(&(mutex[bin_id])); 
		return_value = bin[bin_id].DepositSample(x, dX, _id, _weight); 
		pthread_mutex_unlock(&(mutex[bin_id]));  
	}
	return return_value; 
}

int CStorageHeadPthread::DepositSample(bool if_new_bin, int bin_id, const CSampleIDWeight & sample) 
{
	int return_value;
        if (if_new_bin)
        {
                pthread_mutex_lock(&(mutex[bin_id+number_bins]));
                return_value = bin[bin_id+number_bins].DepositSample(sample);      
                pthread_mutex_unlock(&(mutex[bin_id+number_bins]));
        }
        else
        {
                pthread_mutex_lock(&(mutex[bin_id]));
                return_value = bin[bin_id].DepositSample(sample);
                pthread_mutex_unlock(&(mutex[bin_id]));
        }
        return return_value;
}

void CStorageHeadPthread::Consolidate(int bin_id)
{
	int old_bin_id = bin_id;
        int new_bin_id = bin_id + number_bins;
        if ((int)(bin.size()) > new_bin_id && bin[new_bin_id].GetNumberSamplesGeneratedByFar())
        {
		pthread_mutex_lock(&(mutex[old_bin_id])); 
		pthread_mutex_lock(&(mutex[new_bin_id])); 
                bin[new_bin_id].ChangeFileName(old_bin_id);
                bin[old_bin_id] = bin[new_bin_id];
                bin[old_bin_id].SetBinID(old_bin_id);
		pthread_mutex_unlock(&(mutex[new_bin_id])); 
		pthread_mutex_unlock(&(mutex[old_bin_id])); 
        }
        else
	{
		pthread_mutex_lock(&(mutex[old_bin_id])); 
                bin[old_bin_id].ClearDepositDrawHistory();
		pthread_mutex_unlock(&(mutex[old_bin_id])); 
	}
}


void CStorageHeadPthread::ClearTemporaryBin()
{

	CStorageHead::ClearTemporaryBin();
	for (int i=number_bins; i<(int)(mutex.size()); i++)
		pthread_mutex_destroy(&(mutex[i])); 
	mutex.erase(mutex.begin()+number_bins, mutex.end()); 	
}

