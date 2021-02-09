#include "BalloonFleetSimulator.h"


BalloonFleetSimulator::BalloonFleetSimulator()
{
	std::cout << "Create simulator" << std::endl;
	n_balloons = 1;
	srand(time(NULL));
}

// UI Functions
const char * BalloonFleetSimulator::getTestCasesStr() { return ""; }

void BalloonFleetSimulator::initUI(DrawingUtilitiesClass * DUC) {
	this->DUC = DUC;

	TwAddVarRW(DUC->g_pTweakBar, "N_balloons", TW_TYPE_INT32, &n_balloons, "min=1 max=1000");
}

void BalloonFleetSimulator::reset() {
	std::cout << "Start reset" << std::endl;
	fleet_size = 0;
	int n = getNumberOfBalloons();
	bool game = false;
	if (n_balloons == 1)
		 game = true;
	while (n < n_balloons) {
		std::cout << "Add balloon n'" << n << std::endl;
		n++;
		addBalloon(game);
		//std::cout << "Done" << std::endl;
	}
}

void BalloonFleetSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext) {
	int n = getNumberOfBalloons();
	while (n < n_balloons) {
		n++;
		addBalloon(false);
	}
	while (n > n_balloons) {
		n--;
		fleet_size -= 1;
	}
	for (int i = 0; i < getNumberOfBalloons(); i++){
		//std::cout << "Drawing balloon " << i << " at position " << balloons[i]->getCenter() << std::endl;
		balloons[i]->draw(DUC);
	}
}

void BalloonFleetSimulator::notifyCaseChanged(int testCase) {
	reset();
	std::cout << "Welcome to test case !" << std::endl;
}

void BalloonFleetSimulator::externalForcesCalculations(float timeElapsed) {
	for (int i = 0; i < getNumberOfBalloons(); i++)
		balloons[i]->externalForcesCalculations(timeElapsed);
}

void BalloonFleetSimulator::simulateTimestep(float timeStep) {
	for (int i = 0; i < getNumberOfBalloons(); i++)
		balloons[i]->simulateTimestep(timeStep);
}


void BalloonFleetSimulator::onClick(int x, int y) {}
void BalloonFleetSimulator::onMouse(int x, int y) {}

// Specific Functions
void BalloonFleetSimulator::addBalloon(bool game) { 
	int n = getNumberOfBalloons();
	if (game) {
		balloons[n] = new BalloonSimulator();
	}

	else {
		float x_temp = rand() % (20) -10;
		float x = x_temp;// 10.0 - 5.0;
		std::cout << x << " " << x_temp << std::endl;
		float y = 6;// (rand() % (80)) / 80.0 + 1;
		float radius = 1; // (rand() % (50)) / 100.0 + 0.5;

		balloons[n] = new BalloonSimulator(Vec3(x, y, 0), radius);
		std::cout << "created balloon at position " << balloons[n]->getCenter() << std::endl;
	}
	fleet_size += 1;
}


int BalloonFleetSimulator::getNumberOfBalloons() { return fleet_size; }
Vec3 BalloonFleetSimulator::getPositionOfBalloon(int index) { return balloons[index]->getCenter(); }


void BalloonFleetSimulator::onMouseBtnDown()
{
	for (int i = 0; i < getNumberOfBalloons(); i++)
		balloons[i]->onMouseBtnDown();
}