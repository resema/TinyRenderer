#include <vector>
#include <cmath>
#include <algorithm>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
TGAImage* uvmap = NULL;

Model *model = NULL;
const int width = 800;
const int height = 800;

//
//  Tutorial 1

//void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(p0.x - p1.x)<std::abs(p0.y - p1.y)) {
//		std::swap(p0.x, p0.y);
//		std::swap(p1.x, p1.y);
//		steep = true;
//	}
//	if (p0.x>p1.x) {
//		std::swap(p0, p1);
//	}
//
//	for (int x = p0.x; x <= p1.x; x++) {
//		float t = (x - p0.x) / (float)(p1.x - p0.x);
//		int y = p0.y*(1. - t) + p1.y*t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}

// 
//  Tutorial 2

//void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
//	if (t0.y == t1.y && t0.y == t2.y) return;
//	// sort the vertices t0,t1,t2 lower to upper, bubblesort....
//	if (t0.y>t1.y) std::swap(t0, t1);
//	if (t0.y>t2.y) std::swap(t0, t2);
//	if (t1.y>t2.y) std::swap(t1, t2);
//
//	int total_height = t2.y - t0.y;
//	for (int i = 0; i < total_height; i++) {
//		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
//		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
//		float alpha = (float)i / total_height;
//		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
//		Vec2i A = t0 + (t2 - t0)*alpha;
//		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;
//		if (A.x > B.x) std::swap(A, B);
//		
//		for (int j = A.x; j <= B.x; j++) {
//			image.set(j, t0.y + i, color);
//		}
//	}
//}

//
//  Tutorial 3

//void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
//		std::swap(p0.x, p0.y);
//		std::swap(p1.x, p1.y);
//		steep = true;
//	}
//	if (p0.x > p1.x) {
//		std::swap(p0, p1);
//	}
//
//	for (int x = p0.x; x <= p1.x; x++) {
//		float t = (x - p0.x) / (float)(p1.x - p0.x);
//		int y = p0.y * (1. - t) + p1.y * t + .5;
//		if (steep) {
//			image.set(y, x, color);
//		} 
//		else {
//			image.set(x, y, color);
//		}
//	}
//}

//void rasterize(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color, int ybuffer[]) {
//	if (p0.x > p1.x) {
//		std::swap(p0, p1);
//	}
//	for (int x = p0.x; x <= p1.x; x++) {
//		float t = (x - p0.x) / (float)(p1.x - p0.x);
//		int y = p0.y * (1. - t) + p1.y * t + .5;
//		if (ybuffer[x] < y) {
//			ybuffer[x] = y;
//			image.set(x, 0, color);
//		}
//	}
//}

//int main(int argc, char** argv) {
//	{ // just dumping the 2d scene (yay we have enough dimensions!)
//		TGAImage scene(width, height, TGAImage::RGB);
//
//		// scene "2d mesh"
//		line(Vec2i(20, 34), Vec2i(744, 400), scene, red);
//		line(Vec2i(120, 434), Vec2i(444, 400), scene, green);
//		line(Vec2i(330, 463), Vec2i(594, 200), scene, blue);
//
//		// screen line
//		line(Vec2i(10, 10), Vec2i(790, 10), scene, white);
//
//		scene.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//		scene.write_tga_file("scene.tga");
//	}
//
//	{
//		TGAImage render(width, 16, TGAImage::RGB);
//		int ybuffer[width];
//		for (int i = 0; i<width; i++) {
//			ybuffer[i] = std::numeric_limits<int>::min();
//		}
//		rasterize(Vec2i(20, 34), Vec2i(744, 400), render, red, ybuffer);
//		rasterize(Vec2i(120, 434), Vec2i(444, 400), render, green, ybuffer);
//		rasterize(Vec2i(330, 463), Vec2i(594, 200), render, blue, ybuffer);
//
//		// 1-pixel wide image is bad for eyes, lets widen it
//		for (int i = 0; i<width; i++) {
//			for (int j = 1; j<16; j++) {
//				render.set(i, j, render.get(i, 0));
//			}
//		}
//		render.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//		render.write_tga_file("render.tga");
//	}
//	return 0;
//}

//
// Multi threading approach

Vec3f barycentric(Vec2i *pts, Vec2i P) {
	Vec3f u = cross<float>(Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]), Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]));
	if (std::abs(u[2]) < 1) return Vec3f(-1, 1, 1);
	return Vec3f(1.f - (u.x + u.y) / u.z, u.y/u.z, u.x / u.z);
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
	Vec3f s[2];
	for (int i = 2; i--; ) {
		s[i][0] = C[i] - A[i];
		s[i][1] = B[i] - A[i];
		s[i][2] = A[i] - P[i];
	}
	Vec3f u = cross(s[0], s[1]);
	if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

//void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
//	Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
//	Vec2i bboxmax(0, 0);
//	Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 2; j++) {
//			bboxmin[j] = std::max<int>(0, std::min<int>(bboxmin[j], pts[i][j]));
//			bboxmax[j] = std::min<int>(clamp[j], std::max<int>(bboxmax[j], pts[i][j]));
//		}
//	}
//	Vec2i P;
//	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
//		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
//			Vec3f bc_screen = barycentric(pts, P);
//			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
//			image.set(P.x, P.y, color);
//		}
//	}
//}

//void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
//	if (t0.y == t1.y && t0.y == t2.y) return; // I don't care about degenerate triangles
//	if (t0.y > t1.y) std::swap(t0, t1);
//	if (t0.y > t2.y) std::swap(t0, t2);
//	if (t1.y > t2.y) std::swap(t1, t2);
//	int total_height = t2.y - t0.y;
//	for (int i = 0; i < total_height; i++) {
//		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
//		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
//		float alpha = (float)i / total_height;
//		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
//		Vec2i A =				t0 + (t2 - t0)*alpha;
//		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;
//		if (A.x > B.x) std::swap(A, B);
//		for (int j = A.x; j <= B.x; j++) {
//			image.set(j, t0.y + i, color);
//		}
//	}
//}

void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) {
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i<3; i++) {
		for (int j = 0; j<2; j++) {
			bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
			P.z = 0;
			for (int i = 0; i<3; i++) P.z += pts[i][2] * bc_screen[i];
			if (zbuffer[int(P.x + P.y*width)] < P.z) {
				zbuffer[int(P.x + P.y*width)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
	}
}

void triangle(Vec3f *pts, Vec3f *uvs, float *zbuffer, TGAImage &image, TGAImage &uvmap) {
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f texbboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f texbboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	int uv_width = uvmap.get_width();
	int uv_height = uvmap.get_height();
	Vec2f texclamp(uv_width-1, uv_height);

	for (int i = 0; i<3; i++) {
		for (int j = 0; j<2; j++) {
			bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
			texbboxmin[j] = std::max(0.f, std::min(texbboxmin[j], uvs[i][j] * uv_height));
			texbboxmax[j] = std::min(texclamp[j], std::max(texbboxmax[j], uvs[i][j] * uv_width));
		}
	}
	Vec2f uv_step;
	uv_step.x = (texbboxmax.x - texbboxmin.x) / (bboxmax.x - bboxmin.x);
	uv_step.y = (texbboxmax.y - texbboxmin.y) / (bboxmax.y - bboxmin.y);

	Vec3f P;
	Vec3f T;
	for (P.x = bboxmin.x, T.x = texbboxmin.x; P.x <= bboxmax.x; P.x++, T.x+=uv_step.x) {
		for (P.y = bboxmin.y, T.y = texbboxmin.y; P.y <= bboxmax.y; P.y++, T.y+=uv_step.y) {
			Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
			P.z = 0;
			TGAColor color = uvmap.get(T.x, T.y);
			for (int i = 0; i<3; i++) P.z += pts[i][2] * bc_screen[i];
			if (zbuffer[int(P.x + P.y*width)] < P.z) {
				zbuffer[int(P.x + P.y*width)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
	}
}

Vec3f world2screen(Vec3f v) {
	return Vec3f(int((v.x + 1.)*width / 2. + .5), int((v.y + 1.)*height / 2. + .5), v.z);
}

//
//  Tutorial 3: Include z-Buffer
int main(int argc, char** argv) {
	if (2 == argc) {
		model = new Model(argv[1]);
	}
	else {
		model = new Model("obj/african_head.obj");
	}

	uvmap = new TGAImage();
	uvmap->read_tga_file("tex/african_head_diffuse.tga");
	uvmap->flip_vertically();

	float *zbuffer = new float[width*height];
	for (int i = width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

	TGAImage image(width, height, TGAImage::RGB);
	for (int i = 0; i<model->nfaces(); i++) {
		std::vector<std::vector<int>> face = model->face(i);
		Vec3f world_coords[3];
		Vec3f pts[3];
		Vec3f texs[3];
		for (int j = 0; j < 3; j++) {
			world_coords[j] = model->vert(face[j][0]);
			pts[j] = world2screen(world_coords[j]);
			texs[j] = model->tex(face[j][1]);
		}
		triangle(pts, texs, zbuffer, image, *uvmap);
	}

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	delete model;

	system("pause");
	return 0;
}

//
//  Tutorial 2: Draw face

//int main(int argc, char** argv) {
//	if (2 == argc) {
//		model = new Model(argv[1]);
//	}
//	else {
//		model = new Model("obj/african_head.obj");
//	}
//	TGAImage image(width, height, TGAImage::RGB);
//	Vec3f light_dir(0, 0, -1);
//	//Vec2i pts[3] = { Vec2i(10,10),Vec2i(100,30), Vec2i(190,160) };
//	//triangle(pts, frame, TGAColor(200, 100, 150, 0));
//	for (int i = 0; i < model->nfaces(); i++) {
//		std::vector<int> face = model->face(i);
//		Vec2i screen_coords[3];
//		Vec3f world_coords[3];
//		for (int j = 0; j < 3; j++) {
//			Vec3f v = model->vert(face[j]);
//			screen_coords[j] = Vec2i((v.x + 1.)*width / 2., (v.y + 1.)*height / 2.);
//			world_coords[j] = v;
//		}
//		Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
//		n.normalize();
//		float intensity = n*light_dir;
//		if (intensity > 0) {
//			/*triangle(screen_coords[0], screen_coords[1], screen_coords[2], image,
//				TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));*/
//			triangle(screen_coords, image,
//				TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
//		}
//		//triangle(screen_coords, frame, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
//	}
//	
//	image.flip_vertically(); // to place the origin in the bottom left corner of the image
//	image.write_tga_file("framebuffer.tga");
//
//	return 0;
//}

//
//  Tutorial 1: Draw face verts

//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//		steep = true;
//	}
//	if (x0 > x1) { // make it left-to-right
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//	}
//
//	for (int x = x0; x <= x1; x++) {
//		float t = (x - x0) / (float)(x1 - x0);
//		int y = y0*(1. - t) + y1*t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}

//int main(int argc, char** argv) {
//	if (2 == argc) {
//		model = new Model(argv[1]);
//	}
//	else {
//		model = new Model("obj/african_head.obj");
//	}
//
//	TGAImage image(width, height, TGAImage::RGB);
//	for (int i = 0; i<model->nfaces(); i++) {
//		std::vector<int> face = model->face(i);
//		for (int j = 0; j<3; j++) {
//			Vec3f v0 = model->vert(face[j]);
//			Vec3f v1 = model->vert(face[(j + 1) % 3]);
//			int x0 = (v0.x + 1.)*width / 2.;
//			int y0 = (v0.y + 1.)*height / 2.;
//			int x1 = (v1.x + 1.)*width / 2.;
//			int y1 = (v1.y + 1.)*height / 2.;
//			line(x0, y0, x1, y1, image, white);
//		}
//	}
//
//	image.flip_vertically(); // to have the origin at the left bottom corner of the image
//	image.write_tga_file("output.tga");
//	delete model;
//	return 0;
//}