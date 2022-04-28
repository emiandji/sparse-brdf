


#include "ProgressReporter.h"



CProgressReporter::CProgressReporter(size_t iNumTasks, int iNodeID)
{
	m_iNumTasks = iNumTasks;
	m_iCounter = 0;
	m_iPerc = 0;
	m_iNodeID = iNodeID;
	if(m_iNodeID == 0)
		std::cout << '\r' << m_iPerc << '/' << "100";
}


CProgressReporter::~CProgressReporter()
{

}

void CProgressReporter::Update()
{
#pragma omp atomic
		m_iCounter++;

	int iPerc = (m_iCounter * 100) / m_iNumTasks;
	if (iPerc > m_iPerc)
	{
#pragma omp critical
		{
			m_iPerc = iPerc;
			if (m_iNodeID == 0)
				std::cout << '\r' << m_iPerc << '/' << "100";
		}
	}
}

void CProgressReporter::Done()
{
	if (m_iNodeID == 0)
		std::cout << std::endl;
}
