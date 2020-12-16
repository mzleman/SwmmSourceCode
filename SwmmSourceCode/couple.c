//-----------------
// 添加一些自定义函数并暴露
//
//-----------------
#include "headers.h"
#include "swmm5.h"
#include "export.h"


int swmm_getSWMMSimTime() 
{
	return (int)(TotalDuration/1000.0);
}

int  swmm_getCouplePoints(double blx, double bly, int row, int col, double delta, int* indexs, int* rows, int* cols, double* cqAs,double* weirBs,double* pondedA) 
{
	int n = -1;
	double left, right, top, bottom;
	int xIndex, yIndex;
	left = blx;
	bottom = bly;
	right = blx + col * delta;
	top = bly + row * delta;
	//calculate couple points count
	for (int i = 0; i < Nobjects[NODE]; i++) {
		if (Node[i].type == JUNCTION && Node[i].X != MISSING && Node[i].Y != MISSING) {
			if (Node[i].X > left && Node[i].X<right && Node[i].Y>bottom && Node[i].Y < top) {
				++n;
				indexs[n] = i;
				rows[n] = (top - Node[i].Y) / delta;
				cols[n] = (Node[i].X - left) / delta;
				cqAs[n] = Node[i].cqA*UCF(LENGTH)*UCF(LENGTH);				//后期提供在inp读取该参数时记得转化成ft2
				weirBs[n] = Node[i].weirB*UCF(LENGTH);                      //后期提供在inp读取该参数时记得转化成ft
				pondedA[n] = Node[i].pondedArea*UCF(LENGTH)*UCF(LENGTH);
			}
		}
	}
	return n+1;
}

int swmm_getCouplePointsN(double blx, double bly, int row, int col, double delta) 
{
	int n = 0;
	double left, right, top, bottom;
	int xIndex, yIndex;
	left = blx;
	bottom = bly;
	right = blx + col * delta;
	top = bly + row * delta;
	//calculate couple points count
	for (int i = 0; i < Nobjects[NODE]; i++) {
		if (Node[i].type == JUNCTION && Node[i].X != MISSING && Node[i].Y != MISSING) {
			if (Node[i].X > left && Node[i].X<right && Node[i].Y>bottom && Node[i].Y < top) {
				++n;
			}
		}
	}
	return n;
}

double swmm_getOverflow(int index) 
{
	if (index<0 || index>Nobjects[NODE]) {
		return -1.0;
	}
	return Node[index].overflow*UCF(FLOW);//由cfs换算成cms
}

int swmm_setLatFlow(int index, double q) 
{
	if (index<0 || index>Nobjects[NODE]) {
		return -1;
	}
	Node[index].oldLatFlow = Node[index].newLatFlow;
	Node[index].newLatFlow = q/UCF(FLOW);//常用的是CMS单位制，而swmm计算使用cfs单位制
	return 0;
}

int swmm_setOption_allowPonding(int flag) {
	AllowPonding = flag;
	return 0;
}

double swmm_getNodeHead(int index, int* overflow_flag)
{
	double head;
	*overflow_flag = 0;

	if (index<0 || index>Nobjects[NODE]) 
	{
		return -1.0;
	}

	head = Node[index].newDepth;

	if (head > Node[index].fullDepth) 
		*overflow_flag = 1;

	return (head+ Node[index].invertElev)*UCF(LENGTH);//返回位置水头
}

