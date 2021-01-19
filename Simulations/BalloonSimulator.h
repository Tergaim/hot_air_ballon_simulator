#pragma once
#ifndef BALLOONSIMULATOR_h
#define BALLOONSIMULATOR_h
#include "Simulator.h"
#include <cmath>

struct MassPoint {
	Vec3 position;
	Vec3 Velocity;

	MassPoint(Vec3 position, Vec3 Velocity) : position(position), Velocity(Velocity) {}
};

struct Spring {
	int p1, p2;
	float initial_length;

	Spring(int p1, int p2, float initial_length) : p1(p1), p2(p2), initial_length(initial_length) {}
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
	int addMassPoint(Vec3 position, Vec3 Velocity);
	void addSpring(int masspoint1, int masspoint2, float initialLength);
	int getNumberOfPoints();
	int getNumberOfSprings();
	Vec3 getPositionOfMassPoint(int index);
	Vec3 getVelocityOfMassPoint(int index);
	void applyExternalForce(Vec3 force);

private:
	void eulerImplicitIntegrator(float h);
	void create_envelope(Vec3 center, float radius);

	// Data Attributes
	float m_fMass;
	float m_fStiffness;
	float m_fDamping;
	float gravity = 0;
	float wind = 0;

	int res_envelope = 3;
	vector<MassPoint> envelope_points;
	vector<Spring> envelope_springs;

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
};
#endif