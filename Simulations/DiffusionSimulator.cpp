#include "DiffusionSimulator.h"
#include "pcgsolver.h"
using namespace std;

Grid::Grid(int length, int height, float diff_coef) {
	m = length;
	n = height;
	alpha = diff_coef;
	temperature = new float[m*n];
	for (int i = 0; i < m*n; i++)
		temperature[i] = 0;
}

Grid::Grid() {
	m = 10;
	n = 10;
	alpha = 0.1;
	temperature = new float[m*n];
	for (int i = 0; i < m*n; i++)
		temperature[i] = 0;
}


DiffusionSimulator::DiffusionSimulator()
{
	m_iTestCase = 0;
	m_vfMovableObjectPos = Vec3();
	m_vfMovableObjectFinalPos = Vec3();
	m_vfRotate = Vec3();
	// to be implemented

	T = Grid(10, 10, 0.1);
}

const char * DiffusionSimulator::getTestCasesStr(){
	return "Explicit_solver, Implicit_solver";
}

void DiffusionSimulator::reset(){
		m_mouse.x = m_mouse.y = 0;
		m_trackmouse.x = m_trackmouse.y = 0;
		m_oldtrackmouse.x = m_oldtrackmouse.y = 0;

}

void DiffusionSimulator::initUI(DrawingUtilitiesClass * DUC)
{
	this->DUC = DUC;
	// to be implemented
}

void DiffusionSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	m_vfMovableObjectPos = Vec3(0, 0, 0);
	m_vfRotate = Vec3(0, 0, 0);
	//
	//to be implemented
	//
	switch (m_iTestCase)
	{
	case 0:
		cout << "Explicit solver!\n";
		break;
	case 1:
		cout << "Implicit solver!\n";
		break;
	default:
		cout << "Empty Test!\n";
		break;
	}
}

void DiffusionSimulator::diffuseTemperatureExplicit(float dt) {
	int m = T.m, n = T.n;

	float **tempGrid;
	tempGrid = new float*[m];
	for (int k = 0; k < m; k++)
		tempGrid[k] = new float[n];
	float update_coef = T.alpha * dt; // here dx = dy = 1

	for (int i = 1; i < m - 1; i++)
		for (int j = 1; j < n - 1; j++)
			tempGrid[i][j] = update_coef * (tempGrid[i + 1][j] + tempGrid[i - 1][j] + tempGrid[i][j + 1] + tempGrid[i][j - 1] - 4 * tempGrid[i][j]);
	for (int i = 1; i < m - 1; i++)
		for (int j = 1; j < n - 1; j++)
			T.temperature[n*i + j] = tempGrid[i][j];
	std::cout << T.alpha;
}

void setupB(std::vector<Real>& b, const Grid &T) {
	for (int i = 0; i < T.m*T.n; i++) {
		b.at(i) = T.temperature[i];
	}
}

void fillT(std::vector<Real>& x, Grid &T) {
	for (int i = 0; i < T.m; i++) 
		for(int j = 0; j < T.n; j++){
			T.temperature[T.n*i+j] = x.at(T.n*i+j);
			if (i == 0 || i == T.m - 1 || j == 0 || j == T.n - 1)
				T.temperature[T.n*i+j] = 0;
		}
}

void setupA(SparseMatrix<Real>& A, int m, int n, float update_coef) {//add your own parameters
	// to be implemented
	//setup Matrix A[sizeX*sizeY*sizeZ, sizeX*sizeY*sizeZ]
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);
	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0

	for (int k = 0; k < m*n; k++) {
			A.set_element(k, k, 1); // set diagonal
			// get coordinates for T(k)
			int i = k / n;
			int j = k - n * i;
			if (i == 0 || i == m - 1 || j == 0 || j == n - 1)
				continue;
			A.set_element(k, n * (i+1) + j, update_coef); // T[i+1][j]
			A.set_element(k, n * (i-1) + j, update_coef); // T[i-1][j]
			A.set_element(k, n * i + j+1, update_coef); // T[i][j+1]
			A.set_element(k, n * i + j-1, update_coef); // T[i][j-1]
			A.set_element(k, k, 1+4*update_coef); // T[i][j]
	}
}


void DiffusionSimulator::diffuseTemperatureImplicit(float dt) {//add your own parameters
	// solve A T = b
	// to be implemented
	const int N = T.m*T.n;
	SparseMatrix<Real> *A = new SparseMatrix<Real> (N);
	std::vector<Real> *b = new std::vector<Real>(N);

	setupA(*A, T.m, T.n, T.alpha*dt); // again, here dx = dy = 1
	setupB(*b, T);

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(N);
	for (int j = 0; j < N; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(*A, *b, x, ret_pcg_residual, ret_pcg_iterations, 0);
	// x contains the new temperature values
	fillT(x, T);//copy x to T
}



void DiffusionSimulator::simulateTimestep(float timeStep)
{
	// to be implemented
	// update current setup for each frame
	switch (m_iTestCase)
	{
	case 0:
		diffuseTemperatureExplicit(timeStep);
		break;
	case 1:
		diffuseTemperatureImplicit(timeStep);
		break;
	}
}

void DiffusionSimulator::drawObjects()
{
	// to be implemented
	//visualization
}


void DiffusionSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	drawObjects();
}

void DiffusionSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void DiffusionSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}
