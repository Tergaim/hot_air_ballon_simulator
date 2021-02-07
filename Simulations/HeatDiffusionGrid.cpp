#include "HeatDiffusionGrid.h"
#include <iostream>

////////////////////////////////////////////////////////
#pragma region GRID

HeatDiffusionGrid::Grid::Grid() {
}

HeatDiffusionGrid::Grid::Grid(unsigned int a_iWidth, unsigned int a_iHeight, unsigned int a_iDepth, float a_fDefaultTemperature)
	:
	m_iWidth{ a_iWidth },
	m_iHeight{ a_iHeight },
	m_iDepth{ a_iDepth },
	DATA_CAPACITY{ m_iWidth * m_iHeight * m_iDepth },
	GRID_MAX_INDEX{ DATA_CAPACITY - 1 },
	m_pArrData{ new float[DATA_CAPACITY] {} }
{
	int l_iXLimiter = m_iWidth - 1;
	int l_iYLimiter = m_iHeight - 1;

	for (int l_iZIndex = 0; l_iZIndex < m_iDepth; l_iZIndex++)
	{
		for (int l_iYIndex = 1; l_iYIndex < l_iYLimiter; l_iYIndex++)
		{
			for (int l_iXIndex = 1; l_iXIndex < l_iXLimiter; l_iXIndex++)
			{
				setVal(l_iXIndex, l_iYIndex, l_iZIndex, a_fDefaultTemperature);
			}
		}
	}
}

HeatDiffusionGrid::Grid::~Grid()
{
	delete[] m_pArrData;
	m_pArrData = nullptr;
}

float HeatDiffusionGrid::Grid::getVal(unsigned int a_iX, unsigned int a_iY, unsigned int a_iZ)
{
	if ((a_iX >= m_iWidth) || (a_iY >= m_iHeight) || (a_iZ >= m_iDepth))
	{
		std::cout << "ERROR:: Dimensions out of bounds \n";
	}

	unsigned int l_iIndex = (m_iWidth * m_iHeight * a_iZ) + (a_iY * m_iWidth) + a_iX;
	if (l_iIndex < 0 || l_iIndex > GRID_MAX_INDEX)
	{
		std::cout << "ERROR:: Grid index out of bounds\n";
		return 0.0f;
	}
	return m_pArrData[l_iIndex];
}

void HeatDiffusionGrid::Grid::setVal(unsigned int a_iX, unsigned int a_iY, unsigned int a_iZ, float a_fVal)
{
	if ((a_iX >= m_iWidth) || (a_iY >= m_iHeight) || (a_iZ >= m_iDepth))
	{
		std::cout << "ERROR:: Dimensions out of bounds \n";
	}

	unsigned int l_iIndex = (m_iWidth * m_iHeight * a_iZ) + (a_iY * m_iWidth) + a_iX;
	if (l_iIndex < 0 || l_iIndex > GRID_MAX_INDEX)
	{
		std::cout << "ERROR:: Grid index out of bounds\n";
		return;
	}
	m_pArrData[l_iIndex] = a_fVal;
}

#pragma endregion GRID
////////////////////////////////////////////////////////





HeatDiffusionGrid::HeatDiffusionGrid(unsigned int a_iGridX, unsigned int a_iGridY, unsigned int a_iGridZ, float a_fDefaultTemperature, float a_fMaxTemperature)
	: m_iGridX{ a_iGridX }, m_iGridY{ a_iGridY }, m_iGridZ{ a_iGridZ }, m_fDefaultTemperature{ a_fDefaultTemperature }, MAX_TEMPERATURE{ a_fMaxTemperature }
{
	Reset();
}

//Destructor
HeatDiffusionGrid::~HeatDiffusionGrid()
{
	if (m_pGrid1 != nullptr)
	{
		delete m_pGrid1;
		m_pGrid1 = nullptr;
	}

	if (m_pGrid2 != nullptr)
	{
		delete m_pGrid2;
		m_pGrid2 = nullptr;
	}
}

void HeatDiffusionGrid::Draw(DrawingUtilitiesClass* DUC, Vec3 a_v3Position, float a_fRadius)
{
	Vec3 l_v3SphereScale = { m_fSphereRadius ,m_fSphereRadius ,m_fSphereRadius };

	float l_fDeltaX = (float)m_fCubeDimension / (float)m_iGridX;
	float l_fDeltaY = (float)m_fCubeDimension / (float)m_iGridY;
	float l_fDeltaZ = (float)m_fCubeDimension / (float)m_iGridZ;

	int l_iGrixYLimiter = m_iGridY - 1;
	int l_iGrixXLimiter = m_iGridX - 1;

	for (int l_iIndexZ = 0; l_iIndexZ < m_iGridZ; l_iIndexZ++)
	{
		float l_fXPoZ = (m_fCubeDimension * -0.5f) + l_fDeltaZ * l_iIndexZ;
		for (int l_iIndexY = 1; l_iIndexY < l_iGrixYLimiter; l_iIndexY++)
		{
			float l_fYPos = (m_fCubeDimension * -1.0f) + l_fDeltaY * l_iIndexY*2;
			for (int l_iIndexX = 1; l_iIndexX < l_iGrixXLimiter; l_iIndexX++)
			{
				Vec3 p = a_fRadius*Vec3{ (m_fCubeDimension * -1.0f) + l_fDeltaX * l_iIndexX * 2, l_fYPos, l_fXPoZ };
				if (p.x*p.x + p.y*p.y + p.z*p.z < a_fRadius && p.y + a_v3Position.y > -0.95) {
					float l_fNewValue = m_pNewGrid->getVal(l_iIndexX, l_iIndexY, l_iIndexZ);
					float l_fNormalizedTemperature = l_fNewValue / MAX_TEMPERATURE;


					DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, Vec3(1.0f, l_fNormalizedTemperature, l_fNormalizedTemperature));
					DUC->drawSphere(p + a_v3Position, l_v3SphereScale);
				}
			}
		}
	}
}

float HeatDiffusionGrid::getTemperature()
{
	float l_fReturnAvg = 0;
	int l_iGridYLimiter = m_iGridY - 1;
	int l_iGridXLimiter = m_iGridX - 1;
	int l_iCellCount = 0;

	for (int l_iIndexZ = 0; l_iIndexZ < m_iGridZ; l_iIndexZ++)
		for (int l_iIndexY = 1; l_iIndexY < l_iGridYLimiter; l_iIndexY++)
			for (int l_iIndexX = 1; l_iIndexX < l_iGridXLimiter; l_iIndexX++, l_iCellCount++)
			{
				l_fReturnAvg += m_pNewGrid->getVal(l_iIndexX, l_iIndexY, l_iIndexZ);
			}

	return  l_fReturnAvg / l_iCellCount;
}

void HeatDiffusionGrid::increaseTemperature(const float a_fTemperatureAdd)
{
	if (a_fTemperatureAdd > 0.0f)
	{
		int l_iXLimiter = m_iGridX - 1;
		for (int l_iXIndex = 1; l_iXIndex < l_iXLimiter; l_iXIndex++)
		{
			float l_fTempValue = m_pNewGrid->getVal(l_iXIndex, 1, 0) + a_fTemperatureAdd;
			m_pNewGrid->setVal(l_iXIndex, 1, 0, std::fmin(l_fTempValue + a_fTemperatureAdd, MAX_TEMPERATURE));
		}
	}
}

void HeatDiffusionGrid::Reset()
{

	if (m_pGrid1 != nullptr)
	{
		delete m_pGrid1;
		m_pGrid1 = nullptr;
	}

	if (m_pGrid2 != nullptr)
	{
		delete m_pGrid2;
		m_pGrid2 = nullptr;
	}

	//Grid Setup
	m_pGrid1 = new Grid(m_iGridX, m_iGridY, m_iGridZ, m_fDefaultTemperature);
	m_pGrid2 = new Grid(m_iGridX, m_iGridY, m_iGridZ, m_fDefaultTemperature);

	m_pNewGrid = m_pGrid1; //set with values on start or default on start
	m_pOldGrid = m_pGrid2; // Is set to 0 on start

}


void HeatDiffusionGrid::simulateTimestep(const float& a_fTimeStep, float a_fExternalTemprature)
{
	// Use the old to store the prev frames new grid values
	//Use the new grid to store the calculated new values from the old grid
	Grid* l_pTemp = m_pOldGrid;
	m_pOldGrid = m_pNewGrid;
	m_pNewGrid = l_pTemp;


	int l_iXLimiter = m_iGridX - 1;
	int l_iYLimiter = m_iGridY - 1;

	//Set the boundary cells with the external temperature
	/*for (int l_iZIndex = 0; l_iZIndex < m_iGridZ; l_iZIndex++)
		for (int l_iYIndex = 0; l_iYIndex <= l_iYLimiter; l_iYIndex++)
			for (int l_iXIndex = 0; l_iXIndex <= l_iXLimiter; l_iXIndex++)
				if (l_iYIndex == 0 || l_iYIndex == l_iYLimiter || l_iXIndex == 0 || l_iXIndex == l_iXLimiter)
					m_pOldGrid->setVal(l_iXIndex, l_iYIndex, l_iZIndex, a_fExternalTemprature);*/

	//std::cout << "Before\n";
	//print(m_pOldGrid); 
	diffuseTemperatureImplicit(a_fTimeStep);

	
	//std::cout << "After\n";
	//print(m_pNewGrid);
	//std::cout << "\n\n\n";
}


void HeatDiffusionGrid::setupB(std::vector<Real>& b) {//add your own parameters
	// to be implemented
	//set vector B[sizeX*sizeY]
	int l_iGridSize = m_iGridX * m_iGridY * m_iGridZ;
	for (int l_iIndexZ = 0; l_iIndexZ < m_iGridZ; l_iIndexZ++)
	{
		for (int l_iIndexY = 0; l_iIndexY < m_iGridY; l_iIndexY++)
		{
			int l_iCurrentRowIndex = (l_iIndexZ * m_iGridY * m_iGridX) + (l_iIndexY * m_iGridX);
			for (int l_iIndexX = 0; l_iIndexX < m_iGridX; l_iIndexX++)
			{
				b.at(l_iCurrentRowIndex + l_iIndexX) = m_pOldGrid->getVal(l_iIndexX, l_iIndexY, l_iIndexZ);
			}
		}
	}
}

void HeatDiffusionGrid::fillT(std::vector<Real>& x) {//add your own parameters
	// to be implemented
	//fill T with solved vector x
	//make sure that the temperature in boundary cells stays zero

	int l_iXLimiter = m_iGridX - 1;
	int l_iYLimiter = m_iGridY - 1;
	for (int l_iIndexZ = 0; l_iIndexZ < m_iGridZ; l_iIndexZ++)
	{
		for (int l_iIndexY = 1; l_iIndexY < l_iYLimiter; l_iIndexY++)
		{
			for (int l_iIndexX = 1; l_iIndexX < l_iXLimiter; l_iIndexX++)
			{
				float l_fTemperature = x.at((l_iIndexZ * m_iGridX * m_iGridY) + (l_iIndexY * m_iGridX) + l_iIndexX);
				m_pNewGrid->setVal(l_iIndexX, l_iIndexY, l_iIndexZ, l_fTemperature);
			}
		}
	}
}

void HeatDiffusionGrid::setupA(SparseMatrix<Real>& A, const float& a_fTimeStep)
{	//add your own parameters
	// to be implemented
	//setup Matrix A[sizeX*sizeY*sizeZ, sizeX*sizeY*sizeZ]
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);
	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0

	float l_fCoefficient = m_fDiffusionAlpa * a_fTimeStep;
	int l_iGridXLimiter = m_iGridX - 1;
	int l_iGridYLimiter = m_iGridY - 1;

	int l_iGridXY = m_iGridX * m_iGridY;
	int l_iGridZLimiterIndex = (m_iGridZ - 1) * l_iGridXY;

	float l_fDeltaSqX = (float)m_fCubeDimension / (float)m_iGridX;
	l_fDeltaSqX *= l_fDeltaSqX;
	float l_fCoefficientX = l_fCoefficient / l_fDeltaSqX;

	float l_fDeltaSqY = (float)m_fCubeDimension / (float)m_iGridY;
	l_fDeltaSqY *= l_fDeltaSqY;
	float l_fCoefficientY = l_fCoefficient / l_fDeltaSqY;

	float l_fDeltaSqZ = (float)m_fCubeDimension / (float)m_iGridZ;
	l_fDeltaSqZ *= l_fDeltaSqZ;
	float l_fCoefficientZ = l_fCoefficient / l_fDeltaSqZ;

	for (int l_iIndexZ = 0; l_iIndexZ < m_iGridZ; l_iIndexZ++)
	{
		for (int l_iIndexY = 0; l_iIndexY < m_iGridY; l_iIndexY++)
		{
			for (int l_iIndexX = 0; l_iIndexX < m_iGridX; l_iIndexX++)
			{
				//Check if border index
				if ((l_iIndexX == 0 ||
					l_iIndexX == l_iGridXLimiter ||
					l_iIndexY == 0 ||
					l_iIndexY == l_iGridYLimiter))
				{
					// If border and diagonal
					if (l_iIndexY == l_iIndexX)
					{
						A.set_element(l_iIndexX, l_iIndexY, 1); // set diagonal of boundary values to 1.0
					}
					// If border and not diagonal
					else
					{
						A.set_element(l_iIndexX, l_iIndexY, 0);
					}
				}
				else
				{
					int l_iCurrentIndex = (l_iIndexZ * l_iGridXY) + (l_iIndexY * m_iGridX) + l_iIndexX;

					A.set_element(l_iCurrentIndex, l_iCurrentIndex - 1, -l_fCoefficientX); //x - 1 affects X
					A.set_element(l_iCurrentIndex, l_iCurrentIndex + 1, -l_fCoefficientX); //x + 1 affects X
					A.set_element(l_iCurrentIndex, l_iCurrentIndex - m_iGridX, -l_fCoefficientY); // y - 1 affects Y
					A.set_element(l_iCurrentIndex, l_iCurrentIndex + m_iGridX, -l_fCoefficientY); // y + 1 affects Y

					if (l_iCurrentIndex > l_iGridXY && l_iCurrentIndex < l_iGridZLimiterIndex)
					{
						A.set_element(l_iCurrentIndex, l_iCurrentIndex - l_iGridXY, -l_fCoefficientZ); // z - 1 affects Z
						A.set_element(l_iCurrentIndex, l_iCurrentIndex + l_iGridXY, -l_fCoefficientZ); // z + 1 affects Z
					}
					A.set_element(l_iCurrentIndex, l_iCurrentIndex, 1.0f + 2.0f * (l_fCoefficientX + l_fCoefficientY + l_fCoefficientZ)); // x, y, z affects all dimensions
				}
			}
		}
	}
}


void HeatDiffusionGrid::diffuseTemperatureImplicit(const float& a_fTimeStep) {//add your own parameters
	// solve A T = b
	// to be implemented
	int l_iGridElementCount = m_iGridX * m_iGridY * m_iGridZ;
	//const int N = l_iGridElementCount;//N = sizeX*sizeY*sizeZ
	SparseMatrix<Real> A(l_iGridElementCount);
	std::vector<Real> b(l_iGridElementCount);

	setupA(A, a_fTimeStep);
	setupB(b);

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(l_iGridElementCount);
	for (int j = 0; j < l_iGridElementCount; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);
	// x contains the new temperature values
	fillT(x);//copy x to T
}

void HeatDiffusionGrid::print(Grid* a_pGrid) {
	for (int l_iIndexY = 0; l_iIndexY < m_iGridY; l_iIndexY++) {
		for (int l_iIndexX = 0; l_iIndexX < m_iGridX; l_iIndexX++)
			std::cout << a_pGrid->getVal(l_iIndexX, l_iIndexY, 0) << '\t';
		std::cout << std::endl;
	}
}
