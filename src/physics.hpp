#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <vector>
#include "gl.hpp"

template <class t>
struct listed {
	t object;
	bool valid;
};

template <class t>
class list {
public:
	std::vector<listed<t>> objects;
	
	bool* add(t o) {
		objects.push_back({o, true});
		return &objects.back().valid;
	}
	
	void remove() {
		for (auto o = objects.begin(); o != objects.end();) o->valid ? ++o : o = objects.erase(o);
	}
	
	size_t size() const {
		return objects.size();
	}
	
	const t& operator[](size_t k) const {
		return objects[k].object;
	}
};

class cuboid {
public:
	vec3 v1;
	vec3 v2;
	
	cuboid() : v1(vec3()), v2(vec3()) {}
	
	cuboid(const vec3& v1, const vec3& v2) : v1(v1), v2(v2) {}
	
	void reset(vec3 v, double w, double h, double d) {
		v1 = v;
		v2.x = v.x + w;
		v2.y = v.y + h;
		v2.z = v.z + d;
	}
	
	bool contain(const vec3& v) const {
		return v1.x < v.x && v.x < v2.x && v1.y < v.y && v.y < v2.y && v1.z < v.z && v.z < v2.z;
	}
	
	bool contain(double x, double y, double z) const {
		return v1.x < x && x < v2.x && v1.y < y && y < v2.y && v1.z < z && z < v2.z;
	}
	
	bool hittest(const cuboid& c) {
		return contain(c.v1.x, c.v1.y, c.v1.z) || contain(c.v1.x, c.v1.y, c.v2.z) ||
			contain(c.v1.x, c.v2.y, c.v1.z) || contain(c.v1.x, c.v2.y, c.v2.z) ||
			contain(c.v2.x, c.v1.y, c.v1.z) || contain(c.v2.x, c.v1.y, c.v2.z) ||
			contain(c.v2.x, c.v2.y, c.v1.z) || contain(c.v2.x, c.v2.y, c.v2.z) ||
			c.contain(v1.x, v1.y, v1.z) || c.contain(v1.x, v1.y, v2.z) ||
			c.contain(v1.x, v2.y, v1.z) || c.contain(v1.x, v2.y, v2.z) ||
			c.contain(v2.x, v1.y, v1.z) || c.contain(v2.x, v1.y, v2.z) ||
			c.contain(v2.x, v2.y, v1.z) || c.contain(v2.x, v2.y, v2.z);
	}
};

class solid {
public:
	cuboid box;
	vec3 position;
	double width = 0;
	double height = 0;
	double depth = 0;
	bool* valid = nullptr;
	static list<solid*> world;
	
	solid() {
		valid = world.add(this);
	}
	
	solid(const vec3& p, double w, double h, double d) : position(p), width(w), height(h), depth(d) {
		valid = world.add(this);
		reset();
	}
	
	~solid() {
		*valid = false;
	}
	
	void reset() {
		box.reset(position, width, height, depth);
	}
	
	void collide(solid* s, double x, double y, double z) {
		cuboid& collided = s->box;
		if (s != this && box.hittest(collided)) {
			if (x != 0) position.x = x > 0 ? collided.v1.x - width - buffer : collided.v2.x + buffer;
			if (y != 0) position.y = y > 0 ? collided.v1.y - height - buffer : collided.v2.y + buffer;
			if (z != 0) position.z = z > 0 ? collided.v1.z - depth - buffer : collided.v2.z + buffer;
			reset();
		}
	}
	
	void move(const vec3& d) {
		world.remove();
		position.x += d.x;
		reset();
		for (size_t i = world.size(); i --> 0;) collide(world[i], d.x, 0, 0);
		position.y += d.y;
		reset();
		for (size_t i = world.size(); i --> 0;) collide(world[i], 0, d.y, 0);
		position.z += d.z;
		reset();
		for (size_t i = world.size(); i --> 0;) collide(world[i], 0, 0, d.z);
	}
	
private:
	static double buffer;
};

inline double solid::buffer = 1e-4;
inline list<solid*> solid::world;

class fragility {
public:
	cuboid box;
	vec3 position;
	double width = 0;
	double height = 0;
	double depth = 0;
	
	fragility() {}
	
	fragility(const vec3& p, double w, double h, double d) : position(p), width(w), height(h), depth(d) {
		reset();
	}
	
	void reset() {
		box.reset(position, width, height, depth);
	}
	
	void update() {
		solid::world.remove();
		for (size_t i = solid::world.size(); i --> 0;) if (box.hittest(solid::world[i]->box)) hit(*solid::world[i]);
	}
	
	virtual void hit(const solid& s) {}
};

#endif
