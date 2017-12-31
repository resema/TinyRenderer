#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include "geometry.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<std::vector<std::vector<int>>> faces_;
	std::vector<Vec3f> texs_;

public:
	Model(const char* filename);
	~Model();

	int nverts();
	int nfaces();
	int ntexs();
	Vec3f vert(int i);
	std::vector<std::vector<int>> face(int idx);
	Vec3f tex(int i);
};

#endif // !MODEL_H

