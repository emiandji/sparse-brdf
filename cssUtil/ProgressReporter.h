
#pragma once
#include <omp.h>
#include <string.h>
#include <iostream>


class CProgressReporter
{

	size_t m_iNumTasks;
	size_t m_iCounter;
	int m_iPerc;
	int m_iNodeID;

public:
	CProgressReporter(size_t iNumTasks, int iNodeID);
	~CProgressReporter();

	void Update();
	void Done();
};

