#ifndef FLAT_HPP
#define FLAT_HPP

#include "gl.hpp"
#include "obj.hpp"

class flatshader : public shader {
public:
	bool vextexshader(const model& m, int f, mat<4, 1>* v) override {
		vec3 normal = vec3{v[0][0][0] - v[1][0][0], v[0][1][0] - v[1][1][0], v[0][2][0] - v[1][2][0]} ^
			vec3{v[0][0][0] - v[2][0][0], v[0][1][0] - v[2][1][0], v[0][2][0] - v[2][2][0]};
		if (normal.magnitude() < eps) return false;
		normal = normal.normalize();
		for (int i = 0; i < 3; ++i) v[i] = transform * v[i];
		const objmodel& obj = static_cast<const objmodel&>(m);
		const vec3& rgb = obj.diffusergb[obj.facediffuse[f / 3]];
		if (rgb.x == 0.3569 && normal.y == -1) return false;
		if (rgb.x == 0.0980 && normal.y == 1 && obj.diffusergb.size() == 1) return false;
		if (rgb.x == 0.0588 && normal.y == 1 && obj.diffusergb.size() == 1) return false;
		color = floor(rgb.z * 0xff) + floor(rgb.y * 0xff) * 0x100 + floor(rgb.x * 0xff) * 0x10000;
		return true;
	}
	
	void geometryshader(const model& m, int f, mat<4, 1>* p) override {}
	
	void fragmentshader(const model& m, int f, const vec3& b, vec3& n, uint32_t& c) override {
		c = color;
	}
	
private:
	double color = 0;
};

#endif
