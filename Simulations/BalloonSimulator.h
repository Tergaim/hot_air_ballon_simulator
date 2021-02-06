#pragma once
#ifndef BALLOONSIMULATOR_h
#define BALLOONSIMULATOR_h
#include <cmath>
#include "Simulator.h"
#include "HeatDiffusionGrid.h"

struct MassPoint {
	Vec3 position;
	Vec3 Velocity;
	float m;

	MassPoint(Vec3 position, Vec3 Velocity, float Mass) : position(position), Velocity(Velocity), m(Mass) {}
};

struct Spring {
	int p1, p2;
	float initial_length;
	float stiffness;

	Spring(int p1, int p2, float initial_length, float stiffness) : p1(p1), p2(p2), initial_length(initial_length), stiffness(stiffness) {}
};

class BalloonSimulator :public Simulator {
public:
	// Construtors
	BalloonSimulator();

	// UI Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void collisionCheck();
	void onClick(int x, int y);
	void onMouse(int x, int y);

	// Specific Functions
	void setMass(float mass);
	void setStiffness(float stiffness);
	void setDampingFactor(float damping);
	int addMassPoint(Vec3 position, Vec3 Velocity, float Mass);
	int addSpring(int masspoint1, int masspoint2, float initialLength, float stiffness);
	int getNumberOfPoints();
	int getNumberOfSprings();
	Vec3 getPositionOfMassPoint(int index);
	Vec3 getVelocityOfMassPoint(int index);
	void applyExternalForce(Vec3 force);

	void onMouseBtnDown();

private:
	void eulerImplicitIntegrator(float h);
	void create_envelope(Vec3 center, float radius);

	void generate_pickup();

	// Data Attributes
	int start_net;
	int start_envelope;
	int spring_top_a;
	int spring_top_b;
	bool broken = false;

	float m_fMass;
	float m_fStiffness;
	float m_fDamping;
	float gravity = 0;
	float wind = 0;

	float a = 0.0065; // in K/m 
	float T0 = 298.15; // in K
	float P0 = 101325; // in Pa
	float M = 0.029; // in kg/mol
	float r0 = P0 * M / (8.314 * T0);

	int res_envelope = 10;
	int res_net = 10;
	vector<MassPoint> envelope_points;
	vector<Spring> envelope_springs;
	float cargo = 0;

	// Pickup Attributes
	const int m_fHorizontalBoundary = 10;
	const int m_fTopBoundary = 10;
	const int m_fBotBoundary = 0;

	Vec3 m_pickupPosition;

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;


	HeatDiffusionGrid m_HeatDiffusionGrid;
};
#endif