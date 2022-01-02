#ifndef LEAF_HPP
#define LEAF_HPP

#include <cmath>
#include "gl.hpp"

class leaf {
public:
	instance* graphic = nullptr;
	vec3 velocity;
	vec3 limit;
	vec3 delta;
	double range = 0;
	double speed = 0;
	
	leaf() {}
	
	leaf(instance* g, double r, double s, const vec3& l = {}, const vec3& d = {}) : graphic(g), range(r), speed(s),
		limit(l), delta(d) {}
	
	void reset(const vec3& p) {
		graphic->position = p + vec3{random01(), random01(), random01()} * (range * 2) - vec3{1, 1, 1} * range;
		graphic->rotation = vec3{random01(), random01(), random01()} * M_PI * 2;
		graphic->scale = {.1, 1, .1};
	}
	
	void update(const vec3& p) {
		vec3 acceleration = vec3{random01(), random01(), random01()} * (speed * 2) + vec3{1, 1, 1} * -speed + delta;
		if (fabs(velocity.x + acceleration.x) < limit.x) velocity.x += acceleration.x;
		if (fabs(velocity.y + acceleration.y) < limit.y) velocity.y += acceleration.y;
		if (fabs(velocity.z + acceleration.z) < limit.z) velocity.z += acceleration.z;
		graphic->position += velocity;
		vec3 relative = graphic->position - p;
		if (fabs(relative.x) > range || fabs(relative.y) > range || fabs(relative.z > range)) reset(p);
	}
};

#endif
