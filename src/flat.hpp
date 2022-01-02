#ifndef FLAT_HPP
#define FLAT_HPP

#include <cmath>
#include "gl.hpp"
#include "obj.hpp"

class flatcamera : public camera {
public:
	int cursortype = 0;
	int dashtime = 10;
	int loadtime = 0;
	double distances[16];
	mat<4, 10> cursor[2];
	
	flatcamera() {}
	
	flatcamera(int w, int h, double fy, double zn, double zf) : camera(w, h, fy, zn, zf) {}
	
	void draw(const std::vector<instance*>& is, shader& s, std::vector<uint32_t>& c) {
		int total = 0;
		for (auto& i : is) total += i->data->facevertices.size();
		std::vector<vec3> world(total);
		deferredshading(is, s, world, c);
		for (int i = width * height; i --> 0;) {
			int x = i % width;
			int y = i / width;
			if (loadtime < 15 && zbuffer[i] > distances[loadtime])
				c[i] = zbuffer[i] < distances[loadtime + 1] ? 0xc9ffd5 : 0xa9d9cb;
			double distance = sqrt((x - width / 2) * (x - width / 2) + (y - height / 2) * (y - height / 2));
			double intensity = fmax(0, distance - height / 2) * (10 - dashtime) * .0004;
			c[i] = floor(fmax(c[i] % 0x100 - intensity * 0xff, 0)) +
				floor(fmax(c[i] / 0x100 % 0x100 - intensity * 0xff, 0)) * 0x100 +
				floor(fmax(c[i] / 0x10000 - intensity * 0xff, 0)) * 0x10000;
		}
		if (dashtime < 10) ++dashtime;
		if (loadtime < 15) ++loadtime;
		for (int x = 0; x < 10; ++x)
			for (int y = 0; y < 4; ++y)
				if (cursor[cursortype][y][x] == 1) c[x - 4 + width / 2 + (y - 4 + height / 2) * width] = 0x012840;
	}
};

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
