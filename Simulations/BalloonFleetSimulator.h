#pragma once
#ifndef BALLOONSIMULATOR_h
#define BALLOONSIMULATOR_h
#include <cmath>
#include "Simulator.h"
#include "BalloonSimulator.h"



class BalloonFleetSimulator :public Simulator {
public:
	// Construtors
	BalloonFleetSimulator();

	// UI Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void onClick(int x, int y);
	void onMouse(int x, int y);

	void addBalloon(bool game);
	int getNumberOfBalloons();
	Vec3 getPositionOfBalloon(int index);

	void onMouseBtnDown();

private:
	// Data Attributes
	int n_balloons;
	int fleet_size;
	BalloonSimulator* balloons[1000];

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
};
#endif