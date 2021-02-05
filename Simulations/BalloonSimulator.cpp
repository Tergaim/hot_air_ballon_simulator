#include "BalloonSimulator.h"


BalloonSimulator::BalloonSimulator() 
	:m_HeatDiffusionGrid{ 6,6,1,T0, 500.0f }
{
	m_fMass = 10;
	m_fStiffness = 100;
	m_fDamping = 0.1;
	res_envelope = 10;
	create_envelope(Vec3(0,0,0), 1.0);

	srand(time(NULL));

	generate_pickup();
}

void BalloonSimulator::create_envelope(Vec3 center, float radius) {
	// create net
	std::cout << "Create net with " << res_net << " * " << res_net << " points." << std::endl;
	float baseLength = 1.0 / res_net;
	float mass = m_fMass / res_net;
	start_net = getNumberOfPoints();
	for (int i = 0; i < res_net; i++)
		for (int j = 0; j < res_net; j++)
			addMassPoint(Vec3(0, -baseLength * i, baseLength * j - 0.5), Vec3(0, 0, 0), mass);

	for (int i = 0; i < res_net; i++)
		for (int j = 0; j < res_net; j++) {
			// Structural
			if (j < res_net - 1) addSpring(start_net + res_net * i + j, start_net + res_net * i + j + 1, baseLength, m_fStiffness);
			if (i < res_net - 1) addSpring(start_net + res_net * i + j, start_net + res_net * (i + 1) + j, baseLength, m_fStiffness);
			// Flexion
			if (j < res_net - 2) addSpring(start_net + res_net * i + j, start_net + res_net * i + j + 2, 2 * baseLength, m_fStiffness);
			if (i < res_net - 2) addSpring(start_net + res_net * i + j, start_net + res_net * (i + 2) + j, 2 * baseLength, m_fStiffness);
			// Shear
			if (i < res_net - 1 && j < res_net - 1) addSpring(start_net + res_net * i + j, start_net + res_net * (i + 1) + j + 1, 1.414*baseLength, m_fStiffness);
			if (i < res_net - 1 && j > 0) addSpring(start_net + res_net * i + j, start_net + res_net * (i + 1) + j - 1, 1.414*baseLength, m_fStiffness);
		}

	std::cout << "Create balloon envelope with " << res_envelope << " points." << std::endl;
	double angle_increment = 2 * 3.14159 / res_envelope;
	double angle = 3.14159/2 + angle_increment;
	float initial_length = sqrt(2-2*cos(angle_increment));
	mass = m_fMass / res_envelope;
	start_envelope = addMassPoint(Vec3(0, 2, 0) + center, Vec3(0, 0, 0), mass);
	addMassPoint(Vec3(cos(angle), 1+sin(angle), 0) + center, Vec3(0, 0, 0), mass);
	spring_top_a = addSpring(start_envelope, start_envelope + 1, initial_length, m_fStiffness);
	angle += angle_increment;

	for (int i = 2; i < res_envelope; i++) {
		addMassPoint(Vec3(cos(angle), 1+sin(angle), 0) + center, Vec3(1, 0, 0), mass);
		addSpring(start_envelope + i - 1, start_envelope + i, initial_length, m_fStiffness);
		angle += angle_increment;
	}
	spring_top_b = addSpring(start_envelope, start_envelope + res_envelope-1, initial_length, m_fStiffness);
}

void BalloonSimulator::generate_pickup() {
	int y = (rand() % (2 * m_fVerticalBoundary + 1)) - m_fVerticalBoundary;

	m_pickupPosition = Vec3(m_fHorizontalBoundary, y, 0);
}

// UI Functions
const char * BalloonSimulator::getTestCasesStr() { return ""; }

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

	m_fMass = 100;
	m_fStiffness = 40;
	m_fDamping = 0.1;
	broken = false;
	create_envelope(Vec3(0,0,0), 1);
}

void BalloonSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext) {
	for (int i = 0; i < getNumberOfSprings(); i++) {
		if (broken && (i == spring_top_a || i == spring_top_b))
			continue;
		DUC->beginLine();
		DUC->drawLine(envelope_points[envelope_springs[i].p1].position, Vec3(1, 1, 1), envelope_points[envelope_springs[i].p2].position, Vec3(1, 1, 1));
		DUC->endLine();
	}

	DUC->setUpLighting(Vec3(0.5, 0.5, 0.5), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
	for (int i = 0; i < getNumberOfPoints(); i++) 
		if(!broken || i != start_envelope )
			DUC->drawSphere(envelope_points[i].position, 0.01);


	m_HeatDiffusionGrid.Draw(DUC, {0.0f, 0.0f, 0.0f});

	// Draw Pickup
	DUC->drawSphere(m_pickupPosition, 0.1);
}

void BalloonSimulator::notifyCaseChanged(int testCase) {
	reset();
	collisionCheck();
	std::cout << "Welcome to test case !" << std::endl;
}

void BalloonSimulator::externalForcesCalculations(float timeElapsed) {
	float mass = cargo / (res_envelope+res_net);
	for (int i = 0; i < getNumberOfPoints(); i++) {
		// gravity
		if(i<start_envelope || i >= start_envelope+res_envelope)
			envelope_points[i].Velocity[1] += -timeElapsed * gravity * (envelope_points[i].m + mass);
		// wind
		envelope_points[i].Velocity[0] -= timeElapsed * wind / (envelope_points[i].m + mass);
		// Damping
		envelope_points[i].Velocity *= (1 - m_fDamping);
	}
}

void BalloonSimulator::simulateTimestep(float timeStep) {
	externalForcesCalculations(timeStep);
	eulerImplicitIntegrator(timeStep);
	collisionCheck();
	
	if (m_pickupPosition.x < -m_fHorizontalBoundary) {
		generate_pickup();
	}

	m_pickupPosition -= Vec3(0.1, 0, 0);
}

void BalloonSimulator::collisionCheck() {
	for (int i = 0; i < getNumberOfPoints(); i++)
		if (envelope_points[i].position.y < -0.95) {
			envelope_points[i].position.y = -0.95;
			envelope_points[i].Velocity.y = 0;
		}
	for (int i = start_envelope; i < start_envelope + res_envelope; i++) {
		if (envelope_points[i].position.y > 10) {
			broken = true;
		}
	}
}

void BalloonSimulator::onClick(int x, int y) 
{
	
}


void BalloonSimulator::onMouse(int x, int y) {}

// Specific Functions
void BalloonSimulator::setMass(float mass) { m_fMass = mass; }
void BalloonSimulator::setStiffness(float stiffness) { m_fStiffness = stiffness; }
void BalloonSimulator::setDampingFactor(float damping) { m_fDamping = damping; }

int BalloonSimulator::addMassPoint(Vec3 position, Vec3 Velocity, float Mass) {
	int i = getNumberOfPoints();
	MassPoint p = MassPoint(position, Velocity, Mass);
	envelope_points.push_back(p);
	return i;
}

int BalloonSimulator::addSpring(int masspoint1, int masspoint2, float initialLength, float stiffness) {
	int i = getNumberOfSprings();
	envelope_springs.push_back(Spring(masspoint1, masspoint2, initialLength, stiffness));
	return i;
}

int BalloonSimulator::getNumberOfPoints() { return envelope_points.size(); }
int BalloonSimulator::getNumberOfSprings() { return envelope_springs.size(); }
Vec3 BalloonSimulator::getPositionOfMassPoint(int index) { return envelope_points[index].position; }
Vec3 BalloonSimulator::getVelocityOfMassPoint(int index) { return envelope_points[index].Velocity; }
void BalloonSimulator::applyExternalForce(Vec3 force) {}

void BalloonSimulator::eulerImplicitIntegrator(float h) { // Implements semi-implicit Euler method

	// 1. Update Velocity

	// compute center envelope and mean speed
	Vec3 center = Vec3(0, 0, 0);
	Vec3 velocity = Vec3(0, 0, 0);
	for (int i = start_envelope; i < res_envelope + start_envelope; i++) {
		center += envelope_points[i].position;
		velocity += envelope_points[i].Velocity;
	}
	center /= res_envelope;
	velocity /= res_envelope;


	// Springs
	for (int i = 0; i < getNumberOfSprings(); i++) {
		if (broken && (i == spring_top_a || i == spring_top_b))
			continue;
		MassPoint p1 = envelope_points[envelope_springs[i].p1];
		MassPoint p2 = envelope_points[envelope_springs[i].p2];
		float il = envelope_springs[i].initial_length;
		float l = sqrt(p1.position.squaredDistanceTo(p2.position));
		float force = m_fStiffness * (l - il) / p1.m;
		p1.Velocity = p1.Velocity - h * force*(p1.position - p2.position) / l;
		envelope_points[envelope_springs[i].p1] = p1;
		
		p2.Velocity = p2.Velocity - h * force*(p2.position - p1.position) / l;
		envelope_points[envelope_springs[i].p2] = p2;
	}

	if (!broken) {
		// pression
		Vec3 radius = envelope_points[start_envelope].position - center;
		float radius_length = sqrt(radius.x*radius.x + radius.y*radius.y);
		float V = 2 * 3.14159 * radius_length * radius_length; // 2D so volume = area * 1 m3
		float l_fAltitude = center.y * 100; // In this game world, altitude is (center.y*100) m
		float T_ext = T0 - a * center.y; //temperature
		m_HeatDiffusionGrid.simulateTimestep(h, T_ext);
		float T_int = m_HeatDiffusionGrid.getTemperature();
		float exponent = M * gravity / (8.314 *a);
		float P_ext = P0 * pow(1 - a * center.y / T0, exponent); // pression from https://en.wikipedia.org/wiki/Barometric_formula
		float P_int = (2 * 3.14159 / V) * (r0*8.314*T_int / M); // Ideal gas law: PV=nRT and n = r0*V0/M

		float angle_increment = 2 * 3.14159 / res_envelope;
		for (int i = start_envelope; i < res_envelope + start_envelope; i++) {
			//get vector from center
			Vec3 dir = envelope_points[i].position - center;
			float l = sqrtf(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
			dir /= l;
			float surface = l * angle_increment; // element of surface: area = Arc length * 1 m2
			//envelope_points[i].Velocity += 1e-4 * (h * (P_int - P_ext) / surface) * dir; // 1e-4: scaling factor, TBD
		}

		// Archimedes
		float r = P_ext * M / (8.314 * T_ext); // density in kg.m-3, Ideal gas law
		for (int i = start_envelope; i < res_envelope + start_envelope; i++) {
			//envelope_points[i].Velocity.y += h * (r * V - r0 * 2 * 3.14159) * gravity / envelope_points[i].m;
		}

		// Update air resistance
		float S = sqrt(4 * V / 3.14159); // we are in 2D
		float damping = 0.25 * r * S;
		setDampingFactor(damping/4); //scaling
	}


	// 2. Update Position (after velocity since semi-implicit)
	for (int i = 0; i < getNumberOfPoints(); i++) {
		envelope_points[i].position = envelope_points[i].position + h * envelope_points[i].Velocity;
		if (i > start_envelope && i < start_envelope + res_envelope) {
			envelope_points[i].position.x -= center.x; // Keep balloon at center of screen
			envelope_points[i].position.z = 0; // Keep the balloon 2D
		}
	}

	// Make net stick to balloon
	for (int i = start_net; i < start_net + res_net; i++) {
		int bottom = start_envelope + res_envelope / 2;
		envelope_points[i].position.x = envelope_points[bottom].position.x;
		envelope_points[i].position.y = envelope_points[bottom].position.y;
	}
}

void BalloonSimulator::onMouseBtnDown()
{
	m_HeatDiffusionGrid.increaseTemperature(100.0f);
}