#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_(), texs_() {
	std::ifstream in;
	in.open(filename, std::ifstream::in);
	if (in.fail()) return;
	std::string line;
	while (!in.eof()) {
		std::getline(in, line);
		std::istringstream iss(line.c_str());
		char trash;
		if (!line.compare(0, 2, "v ")) {
			iss >> trash;
			Vec3f v;
			for (int i = 0; i<3; i++) iss >> v[i];
			verts_.push_back(v);
		}
		else if (!line.compare(0, 2, "f ")) {
			std::vector<std::vector<int>> f;
			int itrash, idx, idy;
			iss >> trash;
			while (iss >> idx >> trash >> idy >> trash >> itrash) {
				idx--; // in wavefront obj all indices start at 1, not zero
				idy--;
				f.push_back({ idx, idy });
			}
			faces_.push_back(f);
		}
		else if (!line.compare(0, 4, "vt  ")) {
			iss >> trash >> trash;
			Vec3f t;
			for (int i = 0; i < 3; i++) iss >> t[i];
			texs_.push_back(t);
		}
	}
	std::cerr << "# v# " << verts_.size() << " f# " << faces_.size() << " t# " << texs_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
	return (int)verts_.size();
}

int Model::nfaces() {
	return (int)faces_.size();
}

int Model::ntexs() {
	return (int)texs_.size();
}

std::vector<std::vector<int>> Model::face(int idx) {
	return faces_[idx];
}

Vec3f Model::vert(int i) {
	return verts_[i];
}

Vec3f Model::tex(int i) {
	return texs_[i];
}