#include "BalloonSimulator.h"


BalloonSimulator::BalloonSimulator() 
	:m_HeatDiffusionGrid{ 6,6,1 }
{
	m_fMass = 10;
	m_fStiffness = 40;
	m_fDamping = 0.1;
	res_envelope = 10;
	create_envelope(Vec3(0,0,0), 1.0);
}

void BalloonSimulator::create_envelope(Vec3 center, float radius) {
	std::cout << "Create balloon envelope with " << res_envelope << " points." << std::endl;
	double angle_increment = 2 * 3.14159 / res_envelope;
	double angle = angle_increment;
	float initial_length = sqrt(2-2*cos(angle_increment));
	int start = addMassPoint(Vec3(1, 0, 0), Vec3(1, 0, 0));
	for (int i = 1; i < res_envelope; i++) {
		addMassPoint(Vec3(cos(angle), sin(angle), 0), Vec3(1, 0, 0));
		addSpring(start + i - 1, start + i, initial_length);
		angle += angle_increment;
	}
	addSpring(start, start + res_envelope-1, initial_length);
}

// UI Functions
const char * BalloonSimulator::getTestCasesStr() { return "One-step,Simple Spring,Complex scene"; }

void BalloonSimulator::initUI(DrawingUtilitiesClass * DUC) {
	this->DUC = DUC;

	TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &gravity, "min=0 max=10 step=0.1");
	TwAddVarRW(DUC->g_pTweakBar, "Wind Force", TW_TYPE_FLOAT, &wind, "min=0 step=0.1");
	TwAddVarRW(DUC->g_pTweakBar, "Spring stiffness", TW_TYPE_FLOAT, &m_fStiffness, "min=10 step=10");
	TwAddVarRW(DUC->g_pTweakBar, "Damping", TW_TYPE_FLOAT, &m_fDamping, "min=0 max=1 step=0.1");
}

void BalloonSimulator::reset() {
	envelope_points.clear();
	envelope_springs.clear();
	m_HeatDiffusionGrid.Reset();

	create_envelope(Vec3(0,0,0), 1);
}

void BalloonSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext) {
	for (int i = 0; i < getNumberOfSprings(); i++) {
		DUC->beginLine();
		DUC->drawLine(envelope_points[envelope_springs[i].p1].position, Vec3(1, 1, 1), envelope_points[envelope_springs[i].p2].position, Vec3(1, 1, 1));
		DUC->endLine();
	}

	DUC->setUpLighting(Vec3(0.5, 0.5, 0.5), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
	for (int i = 0; i < getNumberOfPoints(); i++) 
		DUC->drawSphere(envelope_points[i].position, 0.01);


	m_HeatDiffusionGrid.Draw(DUC, {0.0f, 0.0f, 0.0f});
}

void BalloonSimulator::notifyCaseChanged(int testCase) {
	reset();
	std::cout << "Welcome to test case !" << std::endl;
}

void BalloonSimulator::externalForcesCalculations(float timeElapsed) {
	for (int i = 0; i < getNumberOfPoints(); i++) {
		// gravity
		envelope_points[i].Velocity[1] += -timeElapsed * gravity;
		// wind
		envelope_points[i].Velocity[0] += timeElapsed * wind / m_fMass;
		//Damping
		envelope_points[i].Velocity *= (1 - m_fDamping);
	}
}

void BalloonSimulator::simulateTimestep(float timeStep) {
	eulerImplicitIntegrator(timeStep);
	collisionCheck();
	m_HeatDiffusionGrid.simulateTimestep(timeStep);
}

void BalloonSimulator::collisionCheck() {
	for (int i = 0; i < res_envelope; i++)
		if (envelope_points[i].position[1] < -0.95) {
			envelope_points[i].position[1] = -0.95;
			envelope_points[i].Velocity[1] = 0;
		}
}

void BalloonSimulator::onClick(int x, int y) 
{
	m_HeatDiffusionGrid.increaseTemperature(1000.0f);
}


void BalloonSimulator::onMouse(int x, int y) {}

// Specific Functions
void BalloonSimulator::setMass(float mass) { m_fMass = mass; }
void BalloonSimulator::setStiffness(float stiffness) { m_fStiffness = stiffness; }
void BalloonSimulator::setDampingFactor(float damping) { m_fDamping = damping; }

int BalloonSimulator::addMassPoint(Vec3 position, Vec3 Velocity) {
	int i = getNumberOfPoints();
	MassPoint p = MassPoint(position, Velocity);
	envelope_points.push_back(p);
	return i;
}

void BalloonSimulator::addSpring(int masspoint1, int masspoint2, float initialLength) {
	envelope_springs.push_back(Spring(masspoint1, masspoint2, initialLength));
}

int BalloonSimulator::getNumberOfPoints() { return envelope_points.size(); }
int BalloonSimulator::getNumberOfSprings() { return envelope_springs.size(); }
Vec3 BalloonSimulator::getPositionOfMassPoint(int index) { return envelope_points[index].position; }
Vec3 BalloonSimulator::getVelocityOfMassPoint(int index) { return envelope_points[index].Velocity; }
void BalloonSimulator::applyExternalForce(Vec3 force) {}

void BalloonSimulator::eulerImplicitIntegrator(float h) { // Implements semi-implicit Euler method

	// 1. Update Velocity
	for (int i = 0; i < getNumberOfSprings(); i++) {
		MassPoint p1 = envelope_points[envelope_springs[i].p1];
		MassPoint p2 = envelope_points[envelope_springs[i].p2];
		float il = envelope_springs[i].initial_length;
		float l = sqrt(p1.position.squaredDistanceTo(p2.position));
		float force = m_fStiffness * (l - il) / m_fMass;
		p1.Velocity = p1.Velocity - h * force*(p1.position - p2.position) / l;
		envelope_points[envelope_springs[i].p1] = p1;
		
		p2.Velocity = p2.Velocity - h * force*(p2.position - p1.position) / l;
		envelope_points[envelope_springs[i].p2] = p2;
	}
	// 2. Update Position (after velocity since semi-implicit)
	for (int i = 0; i < getNumberOfPoints(); i++) {
		envelope_points[i].position = envelope_points[i].position + h * envelope_points[i].Velocity;
	}
}