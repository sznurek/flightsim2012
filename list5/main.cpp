#include <cstdio>
#include "../lib/vector.h"

using namespace std;
using namespace linalg;

const double GRAV_CONST = 9.81;

struct state {
	vector x;
	vector v;

	state(const vector& x, const vector& v) : x(x), v(v) {}
	state operator*(double h) const {
		return state(this->x * h, this->v * h);
	}

	state operator+(const state& s) const {
		return state(this->x + s.x, this->v + s.v);
	}
};

state F(const state& s) {
	return state(s.v, -0.05 * s.v * s.v.length() + vector(0, 0, -GRAV_CONST));
}

state euler_step(const state& s, double h) {
	state s_acc = F(s);

	return s_acc * h + s;
}

state heun_step(const state& s, double h) {
	state s_acc = F(s);
	state s_euler = F(euler_step(s, h));

	return s + (s_acc + s_euler) * (0.5 * h);	
}

state euler(const state& start, double h) {
	double t = 0;
	
	state s = start;
	do {
		s = euler_step(s, h);
		t += h;
	} while(s.x.z > 0);

	return s;
}

state heun(const state& start, double h) {
	double t = 0;
	
	state s = start;
	do {
		s = heun_step(s, h);
		t += h;
	} while(s.x.z > 0);

	return s;
}

int main(int argc, char** argv) {
	state start(vector(0, 0, 0), vector(1000, 0, 1000));

	double h = 1;
	for(int i = 1; i < 7; i++) {
		h /= 10.0;

		state stop = euler(start, h);
		printf("Euler [%lf]: (%lf, %lf, %lf)\n", h, stop.x.x, stop.x.y, stop.x.z);
	}

	h = 1;
	for(int i = 1; i < 7; i++) {
		h /= 10.0;

		state stop = heun(start, h);
		printf("Heun [%lf]: (%lf, %lf, %lf)\n", h, stop.x.x, stop.x.y, stop.x.z);
	}

	return 0;
}

