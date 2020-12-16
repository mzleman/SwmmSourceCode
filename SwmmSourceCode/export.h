#pragma once
#define SWMM_EXPORT  __declspec(dllexport) __stdcall


#ifdef __cplusplus
extern "C" {
#endif 

	int  SWMM_EXPORT   swmm_open(char* f1, char* f2, char* f3);
	int  SWMM_EXPORT   swmm_start(int saveFlag);

	//int  SWMM_EXPORT   swmm_step(double* elapsedTime);                  
	int  SWMM_EXPORT   swmm_step(double TimeStep);                                                                                                                  //MZ-202007
	int  SWMM_EXPORT   swmm_end(void);
	int  SWMM_EXPORT   swmm_close(void);
	//int  SWMM_EXPORT   swmm_getVersion(void);
	//int  SWMM_EXPORT   swmm_getError(char* errMsg, int msgLen);
	//int  SWMM_EXPORT   swmm_getWarnings(void);

	int  SWMM_EXPORT   swmm_getCouplePoints(double blx, double bly, int row, int col, double delta, int* indexs, int* rows, int* cols, double* cqAs, double* weirBs, double* pondedA);  //MZ-202007
	int  SWMM_EXPORT   swmm_getCouplePointsN(double blx, double bly, int row, int col, double delta);                                                              //MZ-202007
	int  SWMM_EXPORT   swmm_getSWMMSimTime(void);//use after swmm_start                                                                                            //MZ-202007
	double SWMM_EXPORT swmm_getOverflow(int);                                                                                                                      //MZ-202007
	int  SWMM_EXPORT   swmm_setLatFlow(int, double);                                                                                                               //MZ-202007
	int  SWMM_EXPORT   swmm_setOption_allowPonding(int flag);                                                                                                      //MZ-202007
	double SWMM_EXPORT swmm_getNodeHead(int index, int* flag);																									   //MZ-202007
	double SWMM_EXPORT routing_getRoutingStep(int routingModel, double fixedStep);
	void SWMM_EXPORT swmm_nodeFlood(char*);

#ifdef __cplusplus 
}   // matches the linkage specification from above */ 
#endif