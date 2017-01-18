//
//  main.cpp
//  788MsSprClthSimulation
//
//  Created by  on 11/10/30.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//
#include <glut.h>
#include <GL/GL.h>
#include<gl/GLU.h>

#include <math.h>
#include <vector>
#include <iostream>

int eximFlag = 2; //1:ex  2:im
float circlePosX = 2, circlePosY = -7, circlePosZ = 5;
//float circlePosX = 2-1, circlePosY = -7, circlePosZ = 100;
float circleRad = 1.0;
bool circleMove = false;
int cclStrMove = 4000; // when step does the circle start to move

float DEFAULT_DAMPING = -0.125f;

class Mat3
{
public:
	float f[3][3];

	Mat3(float e00, float e01, float e02, float e10, float e11, float e12, float e20, float e21, float e22)
	{
		f[0][0] = e00;
		f[0][1] = e01;
		f[0][2] = e02;
		f[1][0] = e10;
		f[1][1] = e11;
		f[1][2] = e12;
		f[2][0] = e20;
		f[2][1] = e21;
		f[2][2] = e22;
	}

	Mat3() {}
	/*
	void operator+= (const Mat3 &v)
	{
	return Mat3(f[0][0]+v.f[0][0],f[0][1]+v.f[0][1],f[0][2]+v.f[0][2],
	f[1][0]+v.f[1][0],f[1][1]+v.f[1][1],f[1][2]+v.f[1][2],
	f[2][0]+v.f[2][0],f[2][1]+v.f[2][1],f[2][2]+v.f[2][2]);
	}*/

	Mat3 operator/ (const float &a)
	{
		return Mat3(f[0][0] / a, f[0][1] / a, f[0][2] / a,
			f[1][0] / a, f[1][1] / a, f[1][2] / a,
			f[2][0] / a, f[2][1] / a, f[2][2] / a);
	}

	Mat3 operator- (const Mat3 &v)
	{
		return Mat3(f[0][0] - v.f[0][0], f[0][1] - v.f[0][1], f[0][2] - v.f[0][2],
			f[1][0] - v.f[1][0], f[1][1] - v.f[1][1], f[1][2] - v.f[1][2],
			f[2][0] - v.f[2][0], f[2][1] - v.f[2][1], f[2][2] - v.f[2][2]);
	}

	Mat3 operator+ (const Mat3 &v)
	{
		return Mat3(f[0][0] + v.f[0][0], f[0][1] + v.f[0][1], f[0][2] + v.f[0][2],
			f[1][0] + v.f[1][0], f[1][1] + v.f[1][1], f[1][2] + v.f[1][2],
			f[2][0] + v.f[2][0], f[2][1] + v.f[2][1], f[2][2] + v.f[2][2]);
	}

	Mat3 operator* (const float &a)
	{
		return Mat3(f[0][0] * a, f[0][1] * a, f[0][2] * a,
			f[1][0] * a, f[1][1] * a, f[1][2] * a,
			f[2][0] * a, f[2][1] * a, f[2][2] * a);
	}

	/*
	Vec3 operator-()
	{
	return Vec3(-f[0],-f[1],-f[2]);
	}*/
};


class Vec3
{
public:
	float f[3];

	Vec3(float x, float y, float z)
	{
		f[0] = x;
		f[1] = y;
		f[2] = z;
	}

	Vec3() {}

	float length()
	{
		return sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
	}

	Vec3 normalized()
	{
		float l = length();
		return Vec3(f[0] / l, f[1] / l, f[2] / l);
	}

	void operator+= (const Vec3 &v)
	{
		f[0] += v.f[0];
		f[1] += v.f[1];
		f[2] += v.f[2];
	}

	Vec3 operator/ (const float &a)
	{
		return Vec3(f[0] / a, f[1] / a, f[2] / a);
	}

	Vec3 operator- (const Vec3 &v)
	{
		return Vec3(f[0] - v.f[0], f[1] - v.f[1], f[2] - v.f[2]);
	}

	Vec3 operator+ (const Vec3 &v)
	{
		return Vec3(f[0] + v.f[0], f[1] + v.f[1], f[2] + v.f[2]);
	}

	Vec3 operator* (const float &a)
	{
		return Vec3(f[0] * a, f[1] * a, f[2] * a);
	}

	Vec3 operator-()
	{
		return Vec3(-f[0], -f[1], -f[2]);
	}

	Vec3 cross(const Vec3 &v)
	{
		return Vec3(f[1] * v.f[2] - f[2] * v.f[1], f[2] * v.f[0] - f[0] * v.f[2], f[0] * v.f[1] - f[1] * v.f[0]);
	}

	float dot(const Vec3 &v)
	{
		return f[0] * v.f[0] + f[1] * v.f[1] + f[2] * v.f[2];
	}

	Mat3 outerProduct(const Vec3 &v)
	{
		return Mat3(f[0] * v.f[0], f[0] * v.f[1], f[0] * v.f[2],
			f[1] * v.f[0], f[1] * v.f[1], f[1] * v.f[2],
			f[2] * v.f[0], f[2] * v.f[1], f[2] * v.f[2]);
	}

	Vec3 operator* (const Mat3 &a)//Vec3 = a*Vec3
	{
		return Vec3(a.f[0][0] * f[0] + a.f[0][1] * f[1] + a.f[0][2] * f[2],
			a.f[1][0] * f[0] + a.f[1][1] * f[1] + a.f[1][2] * f[2],
			a.f[2][0] * f[0] + a.f[2][1] * f[1] + a.f[2][2] * f[2]);
	}
};


class Particle
{
private:

public:
	bool bMovable;
	Vec3 v3Pos;                                       //粒子位置
	Vec3 v3Velocity;								  //粒子速率
	Vec3 v3PosTemp;									  //粒子临时位置
	Vec3 v3VelocityTemp;							  //粒子临时速率
	Vec3 v3Force;									  //粒子受力
	float fMass;
	Vec3 v3Gravity;
	Mat3 A;                                           //质量矩阵
	Vec3 b;                                           //b=h(f0+h ∂f/∂xv0)
	Vec3 P;                                            //M的对角阵
	Vec3 Pinv;										  //1/P
	Vec3 r;                                            //误差
	Vec3 d;
	Vec3 q;
	Vec3 Pr;
	Vec3 dV;
	Mat3 m3dFdX;									 //力对位置的偏导
	Mat3 m3dFdV;									//力对速度的偏导
	Vec3 v3ballF;
	Vec3 y;
	Vec3 z;

	Particle(Vec3 pos) : v3Pos(pos), bMovable(true){ v3Velocity = Vec3(0, 0, 0); v3Force = Vec3(0, 0, 0); fMass = 1.0; v3Gravity = Vec3(0, -0.0098, 0); y = Vec3(0, 0, 0); z = Vec3(0, 0, 0); }
	Particle() {}

	Vec3& getPos() { return v3Pos; }

	void setUnMovable()
	{
		bMovable = false;
	}

	void timeStep()
	{
		if (bMovable)
		{
			Vec3 temp = v3Pos;
			//Vec3 addElement = Vec3(0.000001,0000001,0.000001);
			//v3Pos = v3Pos + addElement;
		}
	}

};

class Spring
{
private:

public:
	Particle *pP1, *pP2;
	int iSprType;   //1:structure, 2:sheer, 3:bend
	float fKd, fKs;
	float fSprLength;
	float fInvLength;
	float fC;
	Vec3 v3dCdP;
	float fCDot;
	Vec3 v3DltPosition2;
	Mat3 m3dFdX;
	Mat3 m3dFdV;

	Spring(Particle* p1, Particle* p2, int sprType) : pP1(p1), pP2(p2), iSprType(sprType)
	{
		if (sprType == 1)
		{
			fKs = 2.0f;
			fKd = -0.5f;
		}
		else if (sprType == 2)
		{
			fKs = 2.0f;
			fKd = -0.5f;
		}
		else if (sprType == 3)
		{
			fKs = 1.4f;
			fKd = -0.8f;
		}
		Vec3 vec = pP1->getPos() - pP2->getPos();
		fSprLength = vec.length();
	}

	void getPtcOfSpring(Particle** p1, Particle** p2)
	{
		*p1 = pP1;
		*p2 = pP2;
	}
};

class Cloth
{
private:
	int iNmPtcWidth;
	int iNmPtcHeight;

	std::vector<Particle> particles;
	std::vector<Spring> springs;

	Particle* getParticle(int x, int y) { return &particles[y*iNmPtcWidth + x]; }
	void makeSpring(Particle *p1, Particle *p2, int sprType) { springs.push_back(Spring(p1, p2, sprType)); }
public:
	int testNum;
	Vec3 v3CclCenter;
	float fCclR;
	Vec3 v3OldCclCenter;

	Cloth(int width, int height, int nmPtcWidth, int nmPtcHeight) : iNmPtcWidth(nmPtcWidth), iNmPtcHeight(nmPtcHeight)
	{
		testNum = 0;
		v3CclCenter = Vec3(circlePosX, circlePosY, circlePosZ);
		fCclR = circleRad;

		particles.resize(iNmPtcWidth*iNmPtcHeight);

		for (int x = 0; x<iNmPtcWidth; x++)
		{
			for (int z = 0; z<iNmPtcHeight; z++)
			{
				Vec3 pos = Vec3(width * (x / (float)iNmPtcWidth),
					0,
					height * (z / (float)iNmPtcHeight));
				particles[z*iNmPtcWidth + x] = Particle(pos); // insert particle in column x at y'th row
				if ((x == 0 && z == 0) || (x == iNmPtcWidth - 1 && z == 0))
					particles[z*iNmPtcWidth + x].setUnMovable();
			}
		}

		// Connecting immediate neighbor particles with constraints (distance 1 and sqrt(2) in the grid)
		for (int x = 0; x<iNmPtcWidth; x++)
		{
			for (int y = 0; y<iNmPtcHeight; y++)
			{
				if (x<iNmPtcWidth - 1) makeSpring(getParticle(x, y), getParticle(x + 1, y), 1);
				if (y<iNmPtcHeight - 1) makeSpring(getParticle(x, y), getParticle(x, y + 1), 1);
				if (x<iNmPtcWidth - 1 && y<iNmPtcHeight - 1)
				{
					makeSpring(getParticle(x, y), getParticle(x + 1, y + 1), 2);
					makeSpring(getParticle(x + 1, y), getParticle(x, y + 1), 2);
				}
			}
		}


		// Connecting secondary neighbors with constraints (distance 2 and sqrt(4) in the grid)
		for (int x = 0; x<iNmPtcWidth; x++)
		{
			for (int y = 0; y<iNmPtcHeight; y++)
			{
				if (x<iNmPtcWidth - 2) makeSpring(getParticle(x, y), getParticle(x + 2, y), 3);
				if (y<iNmPtcHeight - 2) makeSpring(getParticle(x, y), getParticle(x, y + 2), 3);
				//if (x<num_particles_width-2 && y<num_particles_height-2) makeConstraint(getParticle(x,y),getParticle(x+2,y+2));
				//if (x<num_particles_width-2 && y<num_particles_height-2) makeConstraint(getParticle(x+2,y),getParticle(x,y+2));			
			}
		}


		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->dV.f[0] = 0;
			particle->dV.f[1] = 0;
			particle->dV.f[2] = 0;
		}

	}

	void implicitTimeStep()
	{
		computeImplicitForce();
		implicitIntegration(5.0 / 100.0);
		v3OldCclCenter = v3CclCenter;
		testNum++;
		if (circleMove && testNum > cclStrMove){
			v3CclCenter.f[1] += 0.001;
		}
	}

	void computeImplicitForce()
	{
		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->v3Force = particle->v3Gravity;
			particle->v3Force = particle->v3Force + particle->v3Velocity*DEFAULT_DAMPING;//f=-kv
			//particle->v3Force = particle->v3Force + particle->v3ballF;
		}

		std::vector<Spring>::iterator spring;
		/* for(spring = springs.begin(); spring != springs.end(); spring++ )
		{
		Particle *p1, *p2;
		(*spring).getPtcOfSpring(&p1, &p2);
		if( spring->iSprType == 1)
		{
		Vec3 p1p2 = p2->getPos() - p1->getPos();
		float nowL = p1p2.length();
		float oriL = spring->fSprLength;
		Vec3 ratio = p1p2*(1-oriL/nowL)*0.2;
		if( p1->bMovable )
		p1->v3Pos = p1->v3Pos + ratio;
		if( p2->bMovable )
		p2->v3Pos = p2->v3Pos - ratio;
		}

		}*/

		for (spring = springs.begin(); spring != springs.end(); spring++)
		{
			Particle *p1, *p2;
			(*spring).getPtcOfSpring(&p1, &p2);

			Vec3 dltPosition = p1->getPos() - p2->getPos();             //deltaX
			Vec3 dltVelocity = p1->v3Velocity - p2->v3Velocity;         //deltaV

			float nowLength = dltPosition.length();                     
			spring->fInvLength = 1.0f / nowLength;                      //弹力系数
			spring->fC = nowLength - spring->fSprLength;
			spring->v3dCdP = dltPosition / nowLength;
			spring->fCDot = p1->v3Velocity.dot(spring->v3dCdP * -1) + p2->v3Velocity.dot(spring->v3dCdP);
			spring->v3DltPosition2 = Vec3(dltPosition.f[0] * dltPosition.f[0], dltPosition.f[1] * dltPosition.f[1], dltPosition.f[2] * dltPosition.f[2]);

			float leftTerm = -spring->fKs * (nowLength - spring->fSprLength);
			float rightTerm = spring->fKd * (dltVelocity.dot(dltPosition) / nowLength);
			//float rightTerm = spring->fKd * (dltVelocity.dot(dltVelocity)/nowLength);
			Vec3 springForce = dltPosition.normalized()*(leftTerm + rightTerm);

			p1->v3Force = p1->v3Force + springForce;
			p2->v3Force = p2->v3Force - springForce;

			if (p1->bMovable == false)p1->v3Force = Vec3(0, 0, 0);
			if (p2->bMovable == false)p2->v3Force = Vec3(0, 0, 0);
		}
	}

	void implicitIntegration(float deltaTime)
	{
		computeForceDerivative();

		std::vector<Particle>::iterator particle;
		std::vector<Spring>::iterator spring;

		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			Mat3 M = Mat3(particle->fMass, 0.0, 0.0, 0.0, particle->fMass, 0.0, 0.0, 0.0, particle->fMass);
			particle->A = M - (particle->m3dFdV + particle->m3dFdX*deltaTime)*deltaTime;
			particle->b = (particle->v3Force + ((particle->v3Velocity)*(particle->m3dFdX*deltaTime)) + particle->y*particle->m3dFdX)*deltaTime;
			particle->b = particle->b + particle->z*M;
			particle->P = Vec3(particle->A.f[0][0], particle->A.f[1][1], particle->A.f[2][2]);
			particle->Pinv = Vec3(1.0 / particle->A.f[0][0], 1.0 / particle->A.f[1][1], 1.0 / particle->A.f[2][2]);
			particle->y = Vec3(0, 0, 0);
			particle->z = Vec3(0, 0, 0);
		}

		solveConjugateGradient();

		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->v3VelocityTemp = particle->v3Velocity + (particle->dV*deltaTime);
			particle->v3PosTemp = particle->v3Pos + particle->v3Velocity*deltaTime;
		}

		handleObjectCollision();

		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->v3Velocity = particle->v3VelocityTemp;
			particle->v3Pos = particle->v3PosTemp;
		}
	}

	void computeForceDerivative()
	{
		std::vector<Particle>::iterator particleClear;
		for (particleClear = particles.begin(); particleClear != particles.end(); particleClear++)
		{
			particleClear->m3dFdV = Mat3(0, 0, 0, 0, 0, 0, 0, 0, 0);
			particleClear->m3dFdX = Mat3(0, 0, 0, 0, 0, 0, 0, 0, 0);
		}

		std::vector<Spring>::iterator spring;
		for (spring = springs.begin(); spring != springs.end(); spring++)
		{
			float c = spring->fC;
			Mat3 d2CdP200 = Mat3(-c*spring->v3DltPosition2.f[0] + c, 0.0, 0.0,
				0.0, -c*spring->v3DltPosition2.f[1] + c, 0.0,
				0.0, 0.0, -c*spring->v3DltPosition2.f[2] + c);
			Mat3 d2CdP201 = Mat3(c*spring->v3DltPosition2.f[0] - c, 0.0, 0.0,
				0.0, c*spring->v3DltPosition2.f[1] - c, 0.0,
				0.0, 0.0, c*spring->v3DltPosition2.f[2] - c);
			Mat3 d2CdP211 = d2CdP200;

			Mat3 dp1 = spring->v3dCdP.outerProduct(spring->v3dCdP);
			Mat3 dp2 = spring->v3dCdP.outerProduct(spring->v3dCdP*-1);
			Mat3 dp3 = (spring->v3dCdP*-1).outerProduct(spring->v3dCdP*-1);

			/*
			spring->m3dFdV = spring->m3dFdV + dp1 * -spring->fKd;
			spring->m3dFdV = spring->m3dFdV + dp2 * -spring->fKd;
			spring->m3dFdV = spring->m3dFdV + dp3 * -spring->fKd;
			*/
			Particle *p1, *p2;
			(*spring).getPtcOfSpring(&p1, &p2);

			p1->m3dFdX = p1->m3dFdX + (dp1 + d2CdP200*c) * -spring->fKs - (d2CdP200*spring->fCDot)*spring->fKd;
			p1->m3dFdX = p1->m3dFdX + (dp2 + d2CdP201*c) * -spring->fKs - (d2CdP201*spring->fCDot)*spring->fKd;
			p1->m3dFdX = p1->m3dFdX + (dp3 + d2CdP211*c) * -spring->fKs - (d2CdP211*spring->fCDot)*spring->fKd;

			p1->m3dFdV = p1->m3dFdV + dp1 * -spring->fKd;
			p1->m3dFdV = p1->m3dFdV + dp2 * -spring->fKd;
			p1->m3dFdV = p1->m3dFdV + dp3 * -spring->fKd;

			p2->m3dFdX = p2->m3dFdX + (dp1 + d2CdP200*c) * -spring->fKs - (d2CdP200*spring->fCDot)*spring->fKd;
			p2->m3dFdX = p2->m3dFdX + (dp2 + d2CdP201*c) * -spring->fKs - (d2CdP201*spring->fCDot)*spring->fKd;
			p2->m3dFdX = p2->m3dFdX + (dp3 + d2CdP211*c) * -spring->fKs - (d2CdP211*spring->fCDot)*spring->fKd;

			p2->m3dFdV = p2->m3dFdV + dp1 * -spring->fKd;
			p2->m3dFdV = p2->m3dFdV + dp2 * -spring->fKd;
			p2->m3dFdV = p2->m3dFdV + dp3 * -spring->fKd;
		}
	}

	void solveConjugateGradient()
	{
		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->r = particle->b - particle->dV*particle->A;
			particle->d.f[0] = particle->Pinv.f[0] * particle->r.f[0];
			particle->d.f[1] = particle->Pinv.f[1] * particle->r.f[1];
			particle->d.f[2] = particle->Pinv.f[2] * particle->r.f[2];
			particle->Pr.f[0] = particle->P.f[0] * particle->r.f[0];
			particle->Pr.f[1] = particle->P.f[1] * particle->r.f[1];
			particle->Pr.f[2] = particle->P.f[2] * particle->r.f[2];
		}
		float alpha_new = 0;
		float alpha = 0;
		float beta = 0;
		float delta_old = 0;
		float delta_new = 0;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			delta_new += particle->r.f[0] * particle->Pr.f[0];
			delta_new += particle->r.f[1] * particle->Pr.f[1];
			delta_new += particle->r.f[2] * particle->Pr.f[2];
		}
		float delta0 = delta_new;
		float i = 0;
		const float EPS = 0.001f;
		const float EPS2 = EPS*EPS;
		const int i_max = 10;
		while (i<i_max && delta_new> EPS2*delta0) {
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				particle->q = particle->d*particle->A;
			}
			float dotdq = 0;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				dotdq += particle->d.f[0] * particle->q.f[0];
				dotdq += particle->d.f[1] * particle->q.f[1];
				dotdq += particle->d.f[2] * particle->q.f[2];
			}
			alpha = delta_new / dotdq;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				particle->dV = particle->dV + particle->d*alpha;
				particle->r = particle->r - particle->q*alpha;
			}
			delta_old = delta_new;
			delta_new = 0;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				delta_new += particle->r.f[0] * particle->r.f[0];
				delta_new += particle->r.f[1] * particle->r.f[1];
				delta_new += particle->r.f[2] * particle->r.f[2];
			}
			beta = delta_new / delta_old;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				particle->d = particle->r + particle->d*beta;
			}
			i++;
		}
	}

	void handleObjectCollision()
	{
		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->v3ballF = Vec3(0, 0, 0);
			if ((particle->v3PosTemp - v3CclCenter).length() < fCclR)
			{
				Vec3 p;
				if ((particle->v3Pos - v3CclCenter).length() >= (particle->v3PosTemp - v3CclCenter).length()){
					p = calParticleBallIntersection(particle->v3Pos, particle->v3PosTemp, v3CclCenter, fCclR);
					particle->v3ballF = Vec3(0, 0, 0);
				}
				else{
					//p = calParticleBallIntersection( particle->v3PosTemp, particle->v3Pos, v3CclCenter, fCclR );
					//p = particle->v3PosTemp + (v3CclCenter-v3OldCclCenter);
					p = particle->v3PosTemp + (v3CclCenter - v3OldCclCenter);
					//    particle->v3ballF = (particle->v3Pos-v3OldCclCenter).normalized()*(v3CclCenter-v3OldCclCenter).dot(particle->v3Pos-v3OldCclCenter);
				}
				Vec3 n = (p - v3CclCenter).normalized();

				float nRatio = particle->v3VelocityTemp.dot(n) - (particle->v3PosTemp - p).dot(n);
				nRatio = nRatio<0.0000001 ? 0 : nRatio;
				Vec3 newV = n*nRatio;

				Vec3 tmp = particle->v3VelocityTemp - n*particle->v3VelocityTemp.dot(n);
				if (tmp.f[0] < 0.00001 || tmp.f[1] < 0.00001 || tmp.f[2] < 0.00001){
					//particle->v3VelocityTemp = newV;
					particle->z = newV - particle->v3VelocityTemp;
				}
				else{
					Vec3 vtn = (particle->v3VelocityTemp - n*particle->v3VelocityTemp.dot(n)).normalized();
					float a = particle->v3VelocityTemp.dot(vtn);
					float b = (particle->v3PosTemp - p).dot(vtn);
					a = a<0.0000001 ? 0 : a;
					b = b<0.0000001 ? 0 : b;
					float tRatio = (a - b);
					tRatio = tRatio<0.0000001 ? 0 : tRatio;
					//particle->v3VelocityTemp = newV + vtn*tRatio;
					particle->z = (newV + vtn*tRatio) - particle->v3VelocityTemp;
				}
				particle->y = p - particle->v3PosTemp;
			}

		}
	}

	Vec3 calParticleBallIntersection(Vec3 p1, Vec3 p2, Vec3 c, float r)
	{
		Vec3 p;
		if ((p1 - c).length()<r && (p2 - c).length()<r){
			return (p1 - c).normalized()*(r + 0) + c;
		}
		while (1)
		{
			p = Vec3((p1.f[0] + p2.f[0]) / 2, (p1.f[1] + p2.f[1]) / 2, (p1.f[2] + p2.f[2]) / 2);
			//printf("%f\n", (p-c).length() - r );
			if (fabs((p - c).length() - r) <0.000001)
			{
				//printf("gg\n");
				break;
			}
			else{
				if ((p - c).length() < r) p2 = p;
				else p1 = p;
			}
		}
		return p;
	}

	void explicitTimeStemp()
	{
		computeExplicitForce();
		explicitIntegration(5 / 100.0);
		v3OldCclCenter = v3CclCenter;
		// testNum++;
		if (testNum > 4000){
			v3CclCenter.f[1] += 0.0025;
			//v3CclCenter.f[2]+=0.001;
		}
	}

	void computeExplicitForce()
	{
		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			particle->v3Force = particle->v3Gravity / particle->fMass;
			particle->v3Force = particle->v3Force + particle->v3Velocity*DEFAULT_DAMPING;
		}

		std::vector<Spring>::iterator spring;
		/*for(spring = springs.begin(); spring != springs.end(); spring++ )
		{
		Particle *p1, *p2;
		(*spring).getPtcOfSpring(&p1, &p2);
		if( spring->iSprType == 1)
		{
		Vec3 p1p2 = p2->getPos() - p1->getPos();
		float nowL = p1p2.length();
		float oriL = spring->fSprLength;
		Vec3 ratio = p1p2*(1-oriL/nowL)*0.1;
		if( p1->bMovable )
		p1->v3Pos = p1->v3Pos + ratio;
		if( p2->bMovable )
		p2->v3Pos = p2->v3Pos - ratio;
		}
		}*/
		for (spring = springs.begin(); spring != springs.end(); spring++)
		{
			Particle *p1, *p2;
			(*spring).getPtcOfSpring(&p1, &p2);

			Vec3 dltPosition = p1->getPos() - p2->getPos();
			Vec3 dltVelocity = p1->v3Velocity - p2->v3Velocity;

			float nowLength = dltPosition.length();
			float theH = (nowLength - spring->fSprLength) * spring->fKs;
			float theD = (dltVelocity.dot(dltPosition) * spring->fKd) / nowLength;
			Vec3 sprForce = dltPosition.normalized()*(theD - theH);
			p1->v3Force = p1->v3Force + sprForce;
			p2->v3Force = p2->v3Force - sprForce;
		}
	}

	void explicitIntegration(float deltaTime)
	{
		std::vector<Particle>::iterator particle;
		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			if (particle->bMovable == false)continue;
			float dltTmMass = deltaTime / particle->fMass;

			Vec3 oldV = particle->v3Velocity;
			particle->v3VelocityTemp = particle->v3Velocity + (particle->v3Force*dltTmMass);
			particle->v3PosTemp = particle->v3Pos + oldV*deltaTime;
		}

		handleObjectCollision();

		for (particle = particles.begin(); particle != particles.end(); particle++)
		{
			if (particle->bMovable == false)continue;
			particle->v3Velocity = particle->v3Velocity + (particle->v3Force*(deltaTime / particle->fMass));
			particle->v3Pos = particle->v3PosTemp;
		}
	}

	void drawCloth()
	{
		std::vector<Particle>::iterator particle;

		glColor3f(1.0, 0.0, 0.0);
		glPointSize(3.0);
		glBegin(GL_POINTS);

		for (int x = 0; x<iNmPtcWidth; x++)
		{
			for (int y = 0; y<iNmPtcHeight; y++)
			{
				Particle *p = getParticle(x, y);
				glVertex3fv((GLfloat*)&(p->getPos()));
			}
		}
		glEnd();

		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);
		std::vector<Spring>::iterator spring;
		for (spring = springs.begin(); spring != springs.end(); spring++)
		{
			Particle *p1, *p2;
			(*spring).getPtcOfSpring(&p1, &p2);
			glVertex3fv((GLfloat*)&(p1->getPos()));
			glVertex3fv((GLfloat*)&(p2->getPos()));
		}
		glEnd();



		glBegin(GL_TRIANGLES);
		for (int x = 0; x<iNmPtcWidth - 1; x++)
		{
			for (int y = 0; y<iNmPtcHeight - 1; y++)
			{
				Vec3 color(0, 0, 0);
				color = Vec3(0.6f, 0.2f, 0.2f);


				drawTriangle(getParticle(x + 1, y), getParticle(x, y), getParticle(x, y + 1), color);
				color = Vec3(1.0f, 1.0f, 1.0f);
				drawTriangle(getParticle(x + 1, y + 1), getParticle(x + 1, y), getParticle(x, y + 1), color);



			}
		}
		glEnd();
	}
	void drawTriangle(Particle *p1, Particle *p2, Particle *p3, const Vec3 color)
	{
		glColor3fv((GLfloat*)&color);
		Vec3 n1 = (p1->getPos() - p2->getPos()).cross((p3->getPos() - p2->getPos()));
		glNormal3fv((GLfloat *)&(n1));
		glVertex3fv((GLfloat *)&(p1->getPos()));

		glNormal3fv((GLfloat *)&(n1));
		glVertex3fv((GLfloat *)&(p2->getPos()));

		glNormal3fv((GLfloat *)&(n1));
		glVertex3fv((GLfloat *)&(p3->getPos()));
	}

};

Cloth myCloth(14,10,55,45);
//Cloth myCloth(7, 10, 30, 40);



void init()//(GLvoid)
{
	glClearColor(0.2f, 0.2f, 0.4f, 0.5f);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_COLOR_MATERIAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat lightPos[4] = { -1.0, 1.0, 0.5, 0.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, (GLfloat *)&lightPos);

	glEnable(GL_LIGHT1);

	GLfloat lightAmbient1[4] = { 0.0, 0.0, 0.0, 0.0 };
	GLfloat lightPos1[4] = { 1.0, 0.0, -0.2, 0.0 };
	GLfloat lightDiffuse1[4] = { 0.5, 0.5, 0.3, 0.0 };

	glLightfv(GL_LIGHT1, GL_POSITION, (GLfloat *)&lightPos1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, (GLfloat *)&lightAmbient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, (GLfloat *)&lightDiffuse1);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

}

//OpenGL display is refereed from cg.alexandra.dk/2009/06/02/mosegaards-cloth-simulation-coding-tutorial/
void display(void)
{

	if (eximFlag == 1)
		myCloth.explicitTimeStemp();
	else
		myCloth.implicitTimeStep();


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glDisable(GL_LIGHTING); 	glBegin(GL_POLYGON);
	//glColor3f(0.8f,0.8f,1.0f);
	glColor3f(0.0f, 0.0f, 0.0f);
	glVertex3f(-200.0f, -100.0f, -100.0f);
	glVertex3f(200.0f, -100.0f, -100.0f);
	//glColor3f(0.4f,0.4f,0.8f);	
	glColor3f(0.0f, 0.0f, 0.0f);
	glVertex3f(200.0f, 100.0f, -100.0f);
	glVertex3f(-200.0f, 100.0f, -100.0f);
	glEnd();
	glEnable(GL_LIGHTING);


	glTranslatef(-6.5, 6, -9.0f);
	glRotatef(35, 0, 1, 0);
	//glRotatef(45, 1, 0, 0);
	myCloth.drawCloth();

	glPushMatrix();
	glTranslatef(myCloth.v3CclCenter.f[0], myCloth.v3CclCenter.f[1], myCloth.v3CclCenter.f[2]);
	glColor3f(1.0f, 1.0f, 0.0f);
	glutSolidSphere(myCloth.fCclR - 0.15, 50, 50);
	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (h == 0)
		gluPerspective(80, (float)w, 1.0, 5000.0);
	else
		gluPerspective(80, (float)w / (float)h, 1.0, 5000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 27:
		exit(0);
		break;
	default:
		break;
	}
}


int main(int argc, char * argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(1280, 720);

	glutCreateWindow("Mass Spring Cloth Simulation - 788");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);

	glutMainLoop();
}
