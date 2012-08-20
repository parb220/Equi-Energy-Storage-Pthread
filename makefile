SOURCES = CEES_Pthread.cpp CStorageHeadPthread.cpp test_gaussian_mixture_pthread.cpp TuneEnergyLevlesUpdateStorage.cpp pthread_initialize_simulation.cpp
OBJS = CEES_Pthread.o CStorageHeadPthread.o test_gaussian_mixture_pthread.o TuneEnergyLevlesUpdateStorage.o pthread_initialize_simulation.o
EXECUTABLE = test_gaussian_mixture_pthread

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall 
LIBS := $(LIBS) -lstdc++ -lpthread
LIBS_DIR := $(LIBS_DIR) -L/usr/lib64

EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw
INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
LIBS := $(LIBS) -lgsl -lgslcblas -lm
DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o $(DISTR_MODEL_DIR)/AddScaledLogs.o

SINGLE_CORE_VERSION_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_storage
SINGLE_CORE_VERSION_OBJS = $(SINGLE_CORE_VERSION_DIR)/CEES_Node.o $(SINGLE_CORE_VERSION_DIR)/CPutGetBin.o $(SINGLE_CORE_VERSION_DIR)/CSampleIDWeight.o $(SINGLE_CORE_VERSION_DIR)/CStorageHead.o $(SINGLE_CORE_VERSION_DIR)/MHAdaptive.o $(SINGLE_CORE_VERSION_DIR)/CParameterPackage.o

$(EXECUTABLE) : $(OBJS) $(DISTR_MODEL_OBJS) $(SINGLE_CORE_VERSION_OBJS) 
	$(CPP) $(CPPFLAGS) $(LIBS_DIR) $(LIBS) $(OBJS) $(DISTR_MODEL_OBJS) $(SINGLE_CORE_VERSION_OBJS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@

clean: 
	rm -f *.o $(EXECUTABLE) 
