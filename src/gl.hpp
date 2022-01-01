#ifndef GL_HPP
#define GL_HPP

#include <cmath>
#include <string>
#include <vector>

const double eps = 1e-8;

inline double random01() {
	return rand() / (RAND_MAX + 1.);
}

class vec2 {
public:
	double x = 0;
	double y = 0;
	
	vec2() {}
	
	vec2(double x, double y) : x(x), y(y) {}
	
	vec2 operator-() const {
		return {-x, -y};
	}
	
	void operator+=(const vec2& v) {
		x += v.x;
		y += v.y;
	}
	
	void operator-=(const vec2& v) {
		x -= v.x;
		y -= v.y;
	}
	
	void operator*=(double v) {
		x *= v;
		y *= v;
	}
	
	void operator/=(double v) {
		x /= v;
		y /= v;
	}
	
	double magnitude() const {
		return sqrt(x * x + y * y);
	}
	
	double distance(const vec2& v) const {
		return sqrt((x - v.x) * (x - v.x) + (y - v.y) * (y - v.y));
	}
	
	vec2 normalize() const {
		double d = 1 / sqrt(x * x + y * y);
		return {x * d, y * d};
	}
	
	vec2 rotate(double a) const {
		return {x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)};
	}
	
	static vec2 randomunit() {
		double angle = random01() * M_PI * 2;
		return {cos(angle), sin(angle)};
	}
};

inline vec2 operator+(const vec2& v1, const vec2& v2) {
	return {v1.x + v2.x, v1.y + v2.y};
}

inline vec2 operator-(const vec2& v1, const vec2& v2) {
	return {v1.x - v2.x, v1.y - v2.y};
}

inline vec2 operator*(const vec2& v1, double v2) {
	return {v1.x * v2, v1.y * v2};
}

inline vec2 operator*(double v1, const vec2& v2) {
	return {v2.x * v1, v2.y * v1};
}

inline vec2 operator/(const vec2& v1, double v2) {
	return {v1.x / v2, v1.y / v2};
}

inline double operator*(const vec2& v1, const vec2& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

inline double operator^(const vec2& v1, const vec2& v2) {
	return v1.x * v2.y - v1.y * v2.x;
}

class vec3 {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	
	vec3() {}
	
	vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	
	vec3 operator-() const {
		return {-x, -y, -z};
	}
	
	void operator+=(const vec3& v) {
		x += v.x;
		y += v.y;
		z += v.z;
	}
	
	void operator-=(const vec3& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}
	
	void operator*=(double v) {
		x *= v;
		y *= v;
		z *= v;
	}
	
	void operator/=(double v) {
		x /= v;
		y /= v;
		z /= v;
	}
	
	double magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}
	
	double distance(const vec3& v) const {
		return sqrt((x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z));
	}
	
	vec3 normalize() const {
		double d = 1 / sqrt(x * x + y * y + z * z);
		return {x * d, y * d, z * d};
	}
	
	vec3 rotate(const vec3& v, double a) const {
		return {(cos(a) + (1 - cos(a)) * v.x * v.x) * x +
			((1 - cos(a)) * v.x * v.y - sin(a) * v.z) * y +
			((1 - cos(a)) * v.x * v.z + sin(a) * v.y) * z,
			((1 - cos(a)) * v.x * v.y + sin(a) * v.z) * x +
			(cos(a) + (1 - cos(a)) * v.y * v.y) * y +
			((1 - cos(a)) * v.y * v.z - sin(a) * v.x) * z,
			((1 - cos(a)) * v.x * v.z - sin(a) * v.y) * x +
			((1 - cos(a)) * v.y * v.z + sin(a) * v.x) * y +
			(cos(a) + (1 - cos(a)) * v.z * v.z) * z};
	}
	
	static vec3 randomunit() {
		double angle1 = random01() * M_PI * 2;
		double angle2 = random01() * M_PI * 2;
		return {cos(angle1) * cos(angle2), sin(angle2), sin(angle1) * cos(angle2)};
	}
};

inline vec3 operator+(const vec3& v1, const vec3& v2) {
	return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

inline vec3 operator-(const vec3& v1, const vec3& v2) {
	return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

inline vec3 operator*(const vec3& v1, double v2) {
	return {v1.x * v2, v1.y * v2, v1.z * v2};
}

inline vec3 operator*(double v1, const vec3& v2) {
	return {v2.x * v1, v2.y * v1, v2.z * v1};
}

inline vec3 operator/(const vec3& v1, double v2) {
	return {v1.x / v2, v1.y / v2, v1.z / v2};
}

inline double operator*(const vec3& v1, const vec3& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline vec3 operator^(const vec3& v1, const vec3& v2) {
	return {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
}

template <int r, int c>
class mat {
public:
	double m[r][c];
	
	double* operator[](size_t k) {
		return m[k];
	}
	
	void operator+=(const mat<r, c>& v) {
		for (int i = 0; i < r * c; ++i) m[0][i] += v.m[0][i];
	}

	void operator-=(const mat<r, c>& v) {
		for (int i = 0; i < r * c; ++i) m[0][i] -= v.m[0][i];
	}

	void operator*=(double v) {
		for (int i = 0; i < r * c; ++i) m[0][i] *= v;
	}

	void operator/=(double v) {
		for (int i = 0; i < r * c; ++i) m[0][i] /= v;
	}
	
	static mat<r, c> identity(int v) {
		mat<r, c> matrix;
		for (int i = 0; i < v; ++i) matrix[i][i] = 1;
		return matrix;
	}
	
	mat<c, r> transpose() const {
		mat<c, r> matrix;
		for (int i = 0; i < r; ++i) for (int j = 0; j < c; ++j) matrix[j][i] = m[i][j];
		return matrix;
	}
};

template <int r, int c>
mat<r, c> operator+(const mat<r, c>& v1, const mat<r, c>& v2) {
	mat<r, c> matrix;
	for (int i = 0; i < r * c; ++i) matrix[0][i] = v1.m[0][i] + v2.m[0][i];
	return matrix;
}

template <int r, int c>
mat<r, c> operator-(const mat<r, c>& v1, const mat<r, c>& v2) {
	mat<r, c> matrix;
	for (int i = 0; i < r * c; ++i) matrix[0][i] = v1.m[0][i] - v2.m[0][i];
	return matrix;
}

template <int r, int c>
mat<r, c> operator*(const mat<r, c>& v1, double v2) {
	mat<r, c> matrix;
	for (int i = 0; i < r * c; ++i) matrix[0][i] = v1.m[0][i] * v2;
	return matrix;
}

template <int r, int c>
mat<r, c> operator*(double v1, const mat<r, c>& v2) {
	mat<r, c> matrix;
	for (int i = 0; i < r * c; ++i) matrix[0][i] = v2.m[0][i] * v1;
	return matrix;
}

template <int l1, int l2, int l3>
mat<l1, l3> operator*(const mat<l1, l2>& v1, const mat<l2, l3>& v2) {
	mat<l1, l3> matrix;
	std::fill(matrix[0], matrix[0] + l1 * l3, 0);
	for (int i = 0; i < l1; ++i)
		for (int j = 0; j < l2; ++j)
			for (int k = 0; k < l3; ++k) matrix[i][k] += v1.m[i][j] * v2.m[j][k];
	return matrix;
}

template <int r, int c>
mat<r, c> operator/(const mat<r, c>& v1, double v2) {
	mat<r, c> matrix;
	for (int i = 0; i < r * c; ++i) matrix[0][i] = v1.m[0][i] / v2;
	return matrix;
}

class image {
public:
	int width = 0;
	int height = 0;
	int bytes = 0;
	std::vector<uint8_t> data;
};

class model {
public:
	std::vector<vec3> vertices;
	std::vector<vec2> uvs;
	std::vector<vec3> normals;
	std::vector<int> facevertices;
	std::vector<int> faceuvs;
	std::vector<int> facenormals;
};

class instance {
public:
	vec3 position;
	vec3 rotation;
	vec3 scale;
	model* data = nullptr;
	
	instance(model* d) : data(d), position(vec3()), rotation(vec3()), scale(vec3(1, 1, 1)) {}
	
	instance(model* d, const vec3& p, const vec3& r, const vec3& s) : data(d), position(p), rotation(r), scale(s) {}
};

class light {
public:
	vec3 color = {1, 1, 1};
	
	static bool intersect(const vec3& a, const vec3& b, const vec3& c, const vec3& o, const vec3& d, double l) {
		vec3 ab = b - a;
		vec3 ac = c - a;
		vec3 ao = o - a;
		vec3 p = d ^ ac;
		vec3 q = ao ^ ab;
		double inverse = 1 / (ab * p);
		double u = d * q * inverse;
		double v = ao * p * inverse;
		double t = ac * q * inverse;
		return t > eps && t < l && u > -eps && v > -eps && u + v < 1 + eps;
	}
	
	static bool intersect(const vec3& a, const vec3& b, const vec3& c, const vec3& o, const vec3& d, double l,
		vec3& bc) {
		vec3 ab = b - a;
		vec3 ac = c - a;
		vec3 ao = o - a;
		vec3 p = d ^ ac;
		vec3 q = ao ^ ab;
		double inverse = 1 / (ab * p);
		double u = d * q * inverse;
		double v = ao * p * inverse;
		double t = ac * q * inverse;
		if (t < eps || t > l || u < -eps || v < -eps || u + v > 1 + eps) return false;
		bc.x = 1 - u - v;
		bc.y = v;
		bc.z = u;
		return true;
	}
	
	virtual vec3 lighting(const vec3& p, const vec3& n) const {
		return {};
	}
	
	virtual vec3 cast(const vec3& p, const vec3& n, const std::vector<vec3>& w) const {
		return {};
	};
};

class pointlight : public light {
public:
	vec3 position;
	double intensity = 1;
	double decay = .1;
	
	pointlight(const vec3& p, double i = 1, double d = .1) : position(p), intensity(i), decay(d) {}
	
	vec3 lighting(const vec3& p, const vec3& n) const override {
		double distance = position.distance(p);
		vec3 direction = (position - p).normalize();
		return color * fmax(direction * n * intensity / (distance * distance * decay), 0);
	}
	
	vec3 cast(const vec3& p, const vec3& n, const std::vector<vec3>& w) const override {
		double distance = position.distance(p);
		vec3 direction = (position - p).normalize();
		for (int f = 0; f < w.size(); f += 3)
			if (intersect(w[f], w[f + 1], w[f + 2], p, direction, distance)) return {};
		return color * fmax(direction * n * intensity / (distance * distance * decay), 0);
	}
};

class directionallight : public light {
public:
	vec3 direction;
	double intensity = 1;
	
	directionallight(const vec3& d, double i = 1) : direction(-d), intensity(i) {}
	
	vec3 lighting(const vec3& p, const vec3& n) const override {
		return color * fmax(direction * n * intensity, 0);
	}
	
	vec3 cast(const vec3& p, const vec3& n, const std::vector<vec3>& w) const override {
		for (int f = 0; f < w.size(); f += 3)
			if (intersect(w[f], w[f + 1], w[f + 2], p, direction, limit)) return {};
		return color * fmax(direction * n * intensity, 0);
	}
	
private:
	static double limit;
};

inline double directionallight::limit = 1000;

class arealight : public light {
public:
	vec3 position;
	vec3 toright;
	vec3 tolower;
	vec3 normal;
	double intensity = 1;
	double area = 0;
	
	arealight(const vec3& p, const vec3& tr, const vec3& tl, double i = 1) : position(p), toright(tr), tolower(tl),
		normal((tl ^ tr).normalize()), area((tl ^ tr).magnitude()), intensity(i) {}
	
	static bool emit(const vec3& o, const vec3& d, const std::vector<vec3>& w, vec3& p, vec3& n) {
		double minimum = limit;
		vec3 barycenter;
		bool collided = false;
		for (int f = 0; f < w.size(); f += 3) {
			if (intersect(w[f], w[f + 1], w[f + 2], o, d, limit, barycenter)) {
				vec3 position = w[f] * barycenter.x + w[f + 1] * barycenter.y + w[f + 2] * barycenter.z;
				if (o.distance(position) < minimum) {
					minimum = o.distance(position);
					p = position;
					n = (w[f + 1] - w[f]) ^ (w[f + 2] - w[f]);
				}
				collided = true;
			}
		}
		if (collided) n = n.normalize();
		return collided;
	}
	
	vec3 cast(const vec3& p, const vec3& n, const std::vector<vec3>& w) const override {
		vec3 origin = position + toright * random01() + tolower * random01();
		double distance = origin.distance(p);
		vec3 direction = (origin - p).normalize();
		for (int f = 0; f < w.size(); f += 3)
			if (intersect(w[f], w[f + 1], w[f + 2], p, direction, distance)) return {};
		return color * fmax(intensity * (direction * n) * -(direction * normal) * area / (distance * distance), 0);
	}
	
private:
	static double limit;
};

inline double arealight::limit = 1000;

class shader {
public:
	std::vector<light*> lights;
	mat<4, 4> transform;
	mat<4, 4> translation;
	mat<4, 4> rotation;
	mat<4, 4> scaling;
	
	virtual bool vextexshader(const model& m, int f, mat<4, 1>* v) = 0;
	virtual void geometryshader(const model& m, int f, mat<4, 1>* p) = 0;
	virtual void fragmentshader(const model& m, int f, const vec3& b, vec3& n, uint32_t& c) = 0;
	
protected:
	static const uint8_t* mapping(const image& i, double u, double v) {
		return i.data.data() + (static_cast<int>(u * i.width) + static_cast<int>(v * i.height) * i.width) * i.bytes;
	}
};

struct ginfo {
	bool vaild = false;
	vec3 position;
	vec3 normal;
};

class camera {
public:
	int width = 0;
	int height = 0;
	double fovy = 0;
	double znear = 0;
	double zfar = 0;
	double brdf = .2;
	std::string type = "vertex";
	
	static void lighting(const vec3& r, uint32_t& c) {
		c = floor(c % 0x100 * fmin(r.x, 1)) +
			floor(c / 0x100 % 0x100 * fmin(r.y, 1)) * 0x100 +
			floor(c / 0x10000 * fmin(r.z, 1)) * 0x10000;
	}
	
	static void modeltransform(const instance& i, mat<4, 4>& t, mat<4, 4>& r, mat<4, 4>& s) {
		t = {
			1, 0, 0, i.position.x,
			0, 1, 0, i.position.y,
			0, 0, 1, i.position.z,
			0, 0, 0, 1,
		};
		mat<4, 4> rotationx = {
			1, 0, 0, 0,
			0, cos(i.rotation.x), -sin(i.rotation.x), 0,
			0, sin(i.rotation.x), cos(i.rotation.x), 0,
			0, 0, 0, 1,
		};
		mat<4, 4> rotationy = {
			cos(i.rotation.y), 0, -sin(i.rotation.y), 0,
			0, 1, 0, 0,
			sin(i.rotation.y), 0, cos(i.rotation.y), 0,
			0, 0, 0, 1,
		};
		mat<4, 4> rotationz = {
			cos(i.rotation.z), -sin(i.rotation.z), 0, 0,
			sin(i.rotation.z), cos(i.rotation.z), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1,
		};
		r = rotationx * rotationy * rotationz;
		s = {
			i.scale.x, 0, 0, 0,
			0, i.scale.y, 0, 0,
			0, 0, i.scale.z, 0,
			0, 0, 0, 1,
		};
	}
	
	camera() {}
	
	camera(int w, int h, double fy, double zn, double zf) : width(w), height(h), fovy(fy), znear(zn), zfar(zf) {
		zbuffer = std::vector<double>(width * height);
		gbuffer = std::vector<ginfo>(width * height);
		projection = {
			1 / (tan(fovy / 2) * width / height), 0, 0, 0,
			0, 1 / tan(fovy / 2), 0, 0,
			0, 0, (znear + zfar) / (znear - zfar), 2 * zfar * znear / (znear - zfar),
			0, 0, -1, 0,
		};
	}
	
	void lookat(const vec3& p, const vec3& d, const vec3& u) {
		position = p;
		direction = d;
		up = u;
		vec3 r = d ^ u;
		viewing = {
			r.x, r.y, r.z, -(p * r),
			u.x, u.y, u.z, -(p * u),
			-d.x, -d.y, -d.z, p * d,
			0, 0, 0, 1,
		};
	}
	
	void draw(const std::vector<instance*>& is, const std::vector<light*>& ls, shader& s, std::vector<uint32_t>& c) {
		s.lights = ls;
		int total = 0;
		for (auto& i : is) total += i->data->facevertices.size();
		std::vector<vec3> world(total);
		deferredshading(is, s, world, c);
		if (type == "vertex")
			vertexlighting(ls, c);
		else if (type == "classic")
			classicraytracing(ls, world, c);
		else if (type == "path")
			pathtracing(ls, world, c);
	}
	
	void deferredshading(const std::vector<instance*>& is, shader& s, std::vector<vec3>& w, std::vector<uint32_t>& c) {
		int worldvertex = 0;
		for (int i = width * height; i --> 0;) {
			zbuffer[i] = 1;
			gbuffer[i].vaild = false;
		}
		mat<4, 4> transform = projection * viewing;
		s.transform = transform;
		vec3 barycenters[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
		for (auto& i : is) {
			model& modeldata = *i->data;
			modeltransform(*i, s.translation, s.rotation, s.scaling);
			mat<4, 4> affine = s.translation * s.rotation * s.scaling;
			for (int f = 0; f < modeldata.facevertices.size(); f += 3) {
				mat<4, 1> vertexmatrices[3];
				for (int i = 0; i < 3; ++i) {
					vec3& vertex = modeldata.vertices[modeldata.facevertices[f + i]];
					vertexmatrices[i] = affine * mat<4, 1>{vertex.x, vertex.y, vertex.z, 1.};
					w[worldvertex++] = {vertexmatrices[i][0][0], vertexmatrices[i][1][0], vertexmatrices[i][2][0]};
				}
				// vertex shader
				if (!s.vextexshader(modeldata, f, vertexmatrices)) continue;
				// clipping
				vec3 fixed[4];
				mat<4, 1> primitive[4];
				int endpoints = 0;
				for (int l = 0; l < 3; ++l) {
					mat<4, 1>& point1 = vertexmatrices[l];
					mat<4, 1>& point2 = vertexmatrices[(l + 1) % 3];
					if (point1[3][0] < znear && point2[3][0] < znear) continue;
					const vec3& barycenter1 = barycenters[l];
					const vec3& barycenter2 = barycenters[(l + 1) % 3];
					if (point1[3][0] > znear && point2[3][0] > znear) {
						fixed[endpoints] = barycenter2;
						primitive[endpoints++] = point2;
						continue;
					}
					double weight1 = fabs(point1[3][0] - znear);
					double weight2 = fabs(point2[3][0] - znear);
					double inverse = 1 / (weight1 + weight2);
					if (point1[3][0] < znear && point2[3][0] > znear) {
						fixed[endpoints] = (barycenter1 * weight2 + barycenter2 * weight1) * inverse;
						primitive[endpoints++] = (point1 * weight2 + point2 * weight1) * inverse;
						fixed[endpoints] = barycenter2;
						primitive[endpoints++] = point2;
					} else {
						fixed[endpoints] = (barycenter1 * weight2 + barycenter2 * weight1) * inverse;
						primitive[endpoints++] = (point1 * weight2 + point2 * weight1) * inverse;
					}
				}
				// perspective devision
				vec3 vertices[4];
				for (int i = 0; i < endpoints; ++i)
					vertices[i] = vec3{primitive[i][0][0], primitive[i][1][0], primitive[i][2][0]} / primitive[i][3][0];
				// geomerty shader
				s.geometryshader(modeldata, f, primitive);
				// viewport transform
				for (int i = 0; i < endpoints; ++i) {
					vertices[i].x = vertices[i].x * width / 2 + width / 2;
					vertices[i].y = -vertices[i].y * height / 2 + height / 2;
				}
				// rasterization
				for (int p = 2; p < endpoints; ++p) {
					vec3& vertexa = vertices[0];
					vec3& vertexb = vertices[p - 1];
					vec3& vertexc = vertices[p];
					vec3& fixed1 = fixed[0];
					vec3& fixed2 = fixed[p - 1];
					vec3& fixed3 = fixed[p];
					vec3& worldvertexa = w[worldvertex - 3];
					vec3& worldvertexb = w[worldvertex - 2];
					vec3& worldvertexc = w[worldvertex - 1];
					vec3 vertexz = {vertexa.z, vertexb.z, vertexc.z};
					double primitive1 = 1 / primitive[0][3][0];
					double primitive2 = 1 / primitive[p - 1][3][0];
					double primitive3 = 1 / primitive[p][3][0];
					vec2 v0 = {vertexc.x - vertexa.x, vertexc.y - vertexa.y};
					vec2 v1 = {vertexb.x - vertexa.x, vertexb.y - vertexa.y};
					double dot00 = v0 * v0;
					double dot01 = v0 * v1;
					double dot11 = v1 * v1;
					double inverse = 1 / (dot00 * dot11 - dot01 * dot01);
					vec3 vertexl = vertexa;
					vec3 vertexm = vertexb;
					vec3 vertexu = vertexc;
					if (vertexl.y > vertexm.y) std::swap(vertexl, vertexm);
					if (vertexm.y > vertexu.y) std::swap(vertexm, vertexu);
					if (vertexl.y > vertexm.y) std::swap(vertexl, vertexm);
					double lower = fmax(floor(vertexl.y) + 1, 0);
					double upper = fmin(floor(vertexu.y) + 1, height);
					double median = vertexm.y;
					double inverseml = 1 / (vertexm.y - vertexl.y);
					double inverseum = 1 / (vertexu.y - vertexm.y);
					double inverseul = 1 / (vertexu.y - vertexl.y);
					for (int y = lower; y < upper; ++y) {
						double left = y < median ?
							(vertexl.x * (vertexm.y - y) + vertexm.x * (y - vertexl.y)) * inverseml :
							(vertexm.x * (vertexu.y - y) + vertexu.x * (y - vertexm.y)) * inverseum;
						double right = (vertexl.x * (vertexu.y - y) + vertexu.x * (y - vertexl.y)) * inverseul;
						if (left > right) std::swap(left, right);
						left = fmax(floor(left) + 1, 0);
						right = fmin(floor(right) + 1, width);
						for (int x = left; x < right; ++x) {
							vec2 v2 = {x - vertexa.x + .5, y - vertexa.y + .5};
							double dot02 = v0 * v2;
							double dot12 = v1 * v2;
							double u = (dot11 * dot02 - dot01 * dot12) * inverse;
							double v = (dot00 * dot12 - dot01 * dot02) * inverse;
							vec3 barycenter = {1 - u - v, v, u};
							double z = vertexz * barycenter;
							int location = x + y * width;
							if (z > -1 && z < 1 && z < zbuffer[location] + eps) {
								zbuffer[location] = z;
								barycenter.x *= primitive1;
								barycenter.y *= primitive2;
								barycenter.z *= primitive3;
								barycenter /= barycenter.x + barycenter.y + barycenter.z;
								barycenter = fixed1 * barycenter.x + fixed2 * barycenter.y + fixed3 * barycenter.z;
								gbuffer[location].vaild = true;
								gbuffer[location].position = worldvertexa * barycenter.x +
									worldvertexb * barycenter.y + worldvertexc * barycenter.z;
								// fragment shader
								s.fragmentshader(modeldata, f, barycenter, gbuffer[location].normal, c[location]);
							}
						}
					}
				}
			}
		}
	}
	
	void vertexlighting(const std::vector<light*>& ls, std::vector<uint32_t>& c) {
		for (int i = width * height; i --> 0;) {
			if (!gbuffer[i].vaild) continue;
			vec3 radiance;
			for (auto& l : ls) radiance += l->lighting(gbuffer[i].position, gbuffer[i].normal);
			lighting(radiance, c[i]);
		}
	}
	
	void classicraytracing(const std::vector<light*>& ls, const std::vector<vec3>& w, std::vector<uint32_t>& c) {
		for (int i = width * height; i --> 0;) {
			if (!gbuffer[i].vaild) continue;
			vec3 radiance;
			for (auto& l : ls) radiance += l->cast(gbuffer[i].position, gbuffer[i].normal, w);
			lighting(radiance, c[i]);
		}
	}
	
	void pathtracing(const std::vector<light*>& ls, const std::vector<vec3>& w, std::vector<uint32_t>& c) {
		std::vector<arealight*> lights(ls.size());
		for (size_t i = ls.size(); i --> 0;) lights[i] = static_cast<arealight*>(ls[i]);
		for (int i = width * height; i --> 0;) {
			if (!gbuffer[i].vaild) continue;
			double probability = 0.6;
			double coefficient = 1;
			vec3 radiance;
			vec3 position = gbuffer[i].position;
			vec3 normal = gbuffer[i].normal;
			while (true) {
				for (auto& l : lights) radiance += l->cast(position, normal, w) * coefficient;
				if (probability <= random01()) break;
				vec3 direction = vec3::randomunit();
				if (direction * normal < 0) direction = -direction;
				coefficient *= direction * normal * brdf / probability;
				if (!arealight::emit(position + direction * eps, direction, w, position, normal)) break;
			}
			lighting(radiance, c[i]);
		}
	}
	
protected:
	vec3 position;
	vec3 direction;
	vec3 up;
	mat<4, 4> viewing;
	mat<4, 4> projection;
	std::vector<double> zbuffer;
	std::vector<ginfo> gbuffer;
};

#endif
