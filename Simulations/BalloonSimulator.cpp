#include "BalloonSimulator.h"


BalloonSimulator::BalloonSimulator() 
	:m_HeatDiffusionGrid{ 11,11,1,0, 1000.0f }
{
	game = true;
	m_fMass = 100;
	m_fStiffness = 100;
	m_fDamping = 0.1;
	res_envelope = 100;
	res_net = 1;

	gravity = 9.81;

	original_center = Vec3(0, 6, 0);
	original_radius = 1.0;
	create_envelope(original_center, original_radius);

	srand(time(NULL));

	generate_pickup();
}

BalloonSimulator::BalloonSimulator(Vec3 position, float radius)
	:m_HeatDiffusionGrid{ 11,11,1,0, 1000.0f }
{
	m_fMass = 100;
	m_fStiffness = 100;
	m_fDamping = 0.1;
	res_envelope = 100;
	res_net = 1;

	gravity = 9.81;

	original_center = position;
	original_radius = radius;
	create_envelope(original_center, original_radius);
	std::cout << original_center << std::endl ;

	//srand(time(NULL));

	generate_pickup();
}

void BalloonSimulator::create_envelope(Vec3 a_center, float radius) {
	//std::cout << "Create balloon envelope with " << res_envelope << " points." << std::endl;
	center = a_center;
	double angle_increment = 2 * 3.14159 / res_envelope;
	double angle = 3.14159/2 + angle_increment;
	float initial_length = sqrt(2-2*cos(angle_increment));
	float mass = m_fMass / res_envelope;
	start_envelope = addMassPoint(radius * Vec3(0, 1, 0) + center, Vec3(0, 0, 0), mass);
	mid_envelope = start_envelope + (res_envelope / 2) - 1;
	addMassPoint(radius*Vec3(cos(angle), sin(angle), 0) + center, Vec3(0, 0, 0), mass);
	spring_top_a = addSpring(start_envelope, start_envelope + 1, initial_length, m_fStiffness);
	angle += angle_increment;

	for (int i = 2; i < res_envelope; i++) {
		addMassPoint(radius*Vec3(cos(angle), sin(angle), 0) + center, Vec3(1, 0, 0), mass);
		addSpring(start_envelope + i - 1, start_envelope + i, initial_length, m_fStiffness);
		angle += angle_increment;
	}
	spring_top_b = addSpring(start_envelope, start_envelope + res_envelope-1, initial_length, m_fStiffness);
}

void BalloonSimulator::generate_pickup() {
	int y = (rand() % (m_fTopBoundary - m_fBotBoundary + 1)) + m_fBotBoundary;

	m_pickupPosition = Vec3(m_fHorizontalBoundary, y, 0);
}

// UI Functions
const char * BalloonSimulator::getTestCasesStr() { return ""; }

void BalloonSimulator::initUI(DrawingUtilitiesClass * DUC) {
	if (game) {
		this->DUC = DUC;
		TwAddVarRW(DUC->g_pTweakBar, "Gravity", TW_TYPE_FLOAT, &gravity, "min=0 max=10 step=0.1");
		TwAddVarRW(DUC->g_pTweakBar, "Damping", TW_TYPE_FLOAT, &m_fDamping, "min=0 max=1 step=0.1");
		TwAddVarRW(DUC->g_pTweakBar, "Altitude", TW_TYPE_FLOAT, &l_fAltitude, "");
		TwAddVarRW(DUC->g_pTweakBar, "T_ext", TW_TYPE_FLOAT, &T_ext, "");
		TwAddVarRW(DUC->g_pTweakBar, "T_int", TW_TYPE_FLOAT, &T_int, "");
		TwAddVarRW(DUC->g_pTweakBar, "Score", TW_TYPE_FLOAT, &cargo, "");
	}

}

void BalloonSimulator::reset() {
	envelope_points.clear();
	envelope_springs.clear();
	m_HeatDiffusionGrid.Reset();

	m_fMass = 100;
	m_fStiffness = 40;
	m_fDamping = 0.1;
	broken = false;
	create_envelope(original_center, original_radius);
}

void BalloonSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext) {
	for (int i = 0; i < getNumberOfSprings(); i++) {
		if (broken && (i == spring_top_a || i == spring_top_b))
			continue;
		DUC->beginLine();
		//std::cout << "Began Line" << std::endl;
		DUC->drawLine(envelope_points[envelope_springs[i].p1].position, Vec3(1, 1, 1), envelope_points[envelope_springs[i].p2].position, Vec3(1, 1, 1));
		//std::cout << "Ended Line" << std::endl;
		DUC->endLine();
	}

	DUC->setUpLighting(Vec3(0.5, 0.5, 0.5), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
	//std::cout << "Set up lighting" << std::endl;
	for (int i = 0; i < getNumberOfPoints(); i++) 
		if(!broken || i != start_envelope )
			DUC->drawSphere(envelope_points[i].position, 0.01);

	if(!broken)
		m_HeatDiffusionGrid.Draw(DUC, center, m_fRadius);

	// Draw Pickup
	DUC->setUpLighting(Vec3(0, 0, 1), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
	DUC->drawSphere(m_pickupPosition, 0.1);
}

void BalloonSimulator::draw(DrawingUtilitiesClass* DUC) {
	for (int i = 0; i < getNumberOfSprings(); i++) {
		if (broken && (i == spring_top_a || i == spring_top_b))
			continue;
		DUC->beginLine();
		DUC->drawLine(envelope_points[envelope_springs[i].p1].position, Vec3(1, 1, 1), envelope_points[envelope_springs[i].p2].position, Vec3(1, 1, 1));
		DUC->endLine();
	}

	DUC->setUpLighting(Vec3(0.5, 0.5, 0.5), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
	for (int i = 0; i < getNumberOfPoints(); i++)
		if (!broken || i != start_envelope)
			DUC->drawSphere(envelope_points[i].position, 0.01);

	if (!broken)
		m_HeatDiffusionGrid.Draw(DUC, center, m_fRadius);

	// Draw Pickup
	if (game) {
		DUC->setUpLighting(Vec3(0, 0, 1), Vec3(0.5, 0.5, 0.5), 0.1, Vec3(0.5, 0.5, 0.5));
		DUC->drawSphere(m_pickupPosition, 0.1);
	}
}

void BalloonSimulator::notifyCaseChanged(int testCase) {
	reset();
	collisionCheck();
	std::cout << "Welcome to test case !" << std::endl;
}

void BalloonSimulator::externalForcesCalculations(float timeElapsed) {
	float mass = cargo / res_envelope;
	for (int i = 0; i < getNumberOfPoints(); i++) {
		// gravity
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

	m_pickupPosition -= Vec3(timeStep * 10, 0, 0);
}

void BalloonSimulator::collisionCheck() {
	if (envelope_points[start_envelope].position.y >= 11) {
		broken = true;
	}

	if (envelope_points[mid_envelope].position.y <= -0.8) {
		for (int i = 0; i < getNumberOfPoints(); i++)
			if (envelope_points[i].position.y < -0.95) {
				broken = true;
				envelope_points[i].position.y = -0.95;
				envelope_points[i].Velocity.y = 0;
			}
	}

	m_fRadius = (envelope_points[start_envelope].position.y - envelope_points[mid_envelope].position.y) / 2;

	if (center.squaredDistanceTo(m_pickupPosition) <= (m_fRadius*m_fRadius)) {
		cargo += 10;

		generate_pickup();
	}

	/*
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
	*/
}

void BalloonSimulator::onClick(int x, int y) {}
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
Vec3 BalloonSimulator::getCenter() { return center; }

void BalloonSimulator::applyExternalForce(Vec3 force) {}

void BalloonSimulator::eulerImplicitIntegrator(float h) { // Implements semi-implicit Euler method

	// 1. Update Velocity
	center = (envelope_points[start_envelope].position + envelope_points[mid_envelope].position) / 2;


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

	//pression
	float radius_length = envelope_points[start_envelope].position.y - center.y;
	float V = 2 * 3.14159 * radius_length * radius_length; // 2D so volume = area * 1 m3
	l_fAltitude = center.y * 1000; // In this game world, altitude is (center.y*100) m
	T_ext = T0 - a * l_fAltitude; // we use a linear model for temperature
	float exponent = M * gravity / (8.314 *a);
	float P_ext = P0 * pow(1 - a * l_fAltitude / T0, exponent); // pression from French version of https://en.wikipedia.org/wiki/Barometric_formula
	float P_int = P_ext;

	if (!broken) {
		m_HeatDiffusionGrid.simulateTimestep(h, T_ext);
		T_int = T_ext + m_HeatDiffusionGrid.getTemperature();
		P_int = (2 * 3.14159 / V) * (r0*8.314*T_int / M); // Ideal gas law: PV=nRT and n = r0*V0/M
	}

	float surface = radius_length * 2 * 3.14159; // element of surface: area = Arc length * 1 m2
	for (int i = start_envelope; i < res_envelope + start_envelope; i++) {
		//get vector from center
		Vec3 dir = envelope_points[i].position - center;
		dir /= radius_length;
		envelope_points[i].Velocity += (h * (P_int - P_ext) / (surface*res_envelope)) * dir; 
	}

	if (!broken) {
		// Archimedes
		float r = P_ext * M / (8.314 * T_ext); // density in kg.m-3, Ideal gas law
		float archimedes = (1000/res_envelope) * (r * V - r0 * 2 * 3.14159) * gravity; // 1000/res_envelope: scaling factor for gameplay
		//std::cout << "Archimedes : " << archimedes << std::endl;
		for (int i = start_envelope; i < res_envelope + start_envelope; i++) {
			envelope_points[i].Velocity.y += h * archimedes / envelope_points[i].m;
		}

		// Update air resistance
		float S = sqrt(4 * V / 3.14159); // we are in 2D
		float damping = 0.25 * r * S;
		setDampingFactor(damping/4); //scaling
	}


	// 2. Update Position (after velocity since semi-implicit)

	for (int i = 0; i < getNumberOfPoints(); i++) {
		envelope_points[i].position = envelope_points[i].position + 10 * h * envelope_points[i].Velocity;
		if (i > start_envelope && i < start_envelope + res_envelope) {
			envelope_points[i].position.x += original_center.x - center.x; // Keep balloon at center of screen
			envelope_points[i].position.z = 0; // Keep the balloon 2D
		}
	}
}

void BalloonSimulator::onMouseBtnDown()
{
	m_HeatDiffusionGrid.increaseTemperature(100.0f);
}