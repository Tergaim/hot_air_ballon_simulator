#pragma once
#include "Simulator.h"
#include "vectorbase.h"
#include "pcgsolver.h"

class HeatDiffusionGrid
{
private:
	//impement your own grid class for saving grid data
	class Grid {
	public:
		// Construtors
		Grid();
		Grid(unsigned int a_iWidth, unsigned int a_iHeight, unsigned int a_iDepth, float a_fDefaultTemperature);
		~Grid();

		float getVal(unsigned int a_iX, unsigned int a_iY, unsigned int a_iZ);
		void setVal(unsigned int a_iX, unsigned int a_iY, unsigned int a_iZ, float a_fVal);

		inline unsigned int getWidth() { return m_iWidth; }
		inline unsigned int getHeight() { return m_iHeight; }
		inline unsigned int getDepth() { return m_iDepth; }

	private:
		// Attributes
		unsigned int m_iWidth = 0;
		unsigned int m_iHeight = 0;
		unsigned int m_iDepth = 0;
		const unsigned int DATA_CAPACITY = 0;
		const unsigned int GRID_MAX_INDEX = 0;
		float* m_pArrData = nullptr;
	};


	/// Member variables
	unsigned int m_iGridX = 16;
	unsigned int m_iGridY = 16;
	unsigned int m_iGridZ = 1;
	float m_fSphereRadius = 0.04f;
	float m_fCubeDimension = 1.0f;
	Grid* m_pGrid1 = nullptr; //save results of every time step
	Grid* m_pGrid2 = nullptr; //save results of every time step
	Grid* m_pOldGrid = nullptr;
	Grid* m_pNewGrid = nullptr;
	float m_fDiffusionAlpa = 0.9f;
	float m_fDefaultTemperature = 0.0f;

	const float MAX_TEMPERATURE = 100.0f;

	void setupB(std::vector<Real>& b);
	void fillT(std::vector<Real>& x);
	void setupA(SparseMatrix<Real>& A, const float& a_fTimeStep);
	void diffuseTemperatureImplicit(const float& a_fTimeStep);

public:
	HeatDiffusionGrid(unsigned int a_iGridX, unsigned int a_iGridY, unsigned int a_iGridZ, float a_fDefaultTemperature, float a_fMaxTemperature);
	~HeatDiffusionGrid();

	//Draws the grid at given position
	void Draw(DrawingUtilitiesClass* DUC, Vec3 a_v3Position);

	//Returns the temperature of the upper row of the diffusion grid
	float getTemperature();

	//Increases the temperature of the lower row of the grid by a_fTemperatureAdd 
	void increaseTemperature(const float a_fTemperatureAdd);

	//Recreates the grid with the set size and sets all the temperature values to 0
	void Reset();

	//Calculates the new temperature of each cell in the grid
	void simulateTimestep(const float& a_fTimeStep, float a_fExternalTemprature);

	void print(Grid* a_pGrid);
};