#ifndef OBJ_HPP
#define OBJ_HPP

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "gl.hpp"

class objmodel : public model {
public:
	std::vector<vec3> diffusergb;
	std::vector<int> facediffuse;
};

inline void readmtl(const std::string& f, objmodel& m, std::unordered_map<std::string, int>& im) {
	std::ifstream in;
	in.open(f, std::ifstream::in);
	std::string line;
	while (!in.eof()) {
		std::getline(in, line);
		std::istringstream stream(line.c_str());
		std::string keyword;
		stream >> keyword;
		if (keyword == "newmtl") {
			std::string name;
			stream >> name;
			im.insert({name, im.size()});
		} else if (keyword == "Kd") {
			vec3 rgb;
			stream >> rgb.x >> rgb.y >> rgb.z;
			m.diffusergb.push_back(rgb);
		}
	}
	in.close();
}

inline void readobj(const std::string& f, objmodel& m) {
	std::ifstream in;
	in.open(f, std::ifstream::in);
	std::string line;
	std::unordered_map<std::string, int> indexmap;
	int index = 0;
	while (!in.eof()) {
		std::getline(in, line);
		std::istringstream stream(line.c_str());
		std::string keyword;
		stream >> keyword;
		if (keyword == "v") {
			vec3 vertex;
			stream >> vertex.x >> vertex.y >> vertex.z;
			m.vertices.push_back(vertex);
		} else if (keyword == "vn") {
			vec3 normal;
			stream >> normal.x >> normal.y >> normal.z;
			m.normals.push_back(normal);
		} else if (keyword == "vt") {
			vec2 uv;
			stream >> uv.x >> uv.y;
			m.uvs.push_back(uv);
		} else if (keyword == "f") {
			char _;
			int vertex;
			int uv;
			int normal;
			for (int i = 0; i < 3; ++i) {
				stream >> vertex >> _ >> uv >> _ >> normal;
				m.facevertices.push_back(vertex - 1);
				m.faceuvs.push_back(uv - 1);
				m.facenormals.push_back(normal - 1);
			}
			m.facediffuse.push_back(index);
		} else if (keyword == "mtllib") {
			std::string file;
			stream >> file;
			readmtl(f.substr(0, f.rfind('/') + 1) + file, m, indexmap);
		} else if (keyword == "usemtl") {
			std::string name;
			stream >> name;
			index = indexmap[name];
		}
	}
	in.close();
}

#endif
