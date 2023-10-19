#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <iostream>
#include <ostream>
#include <unistd.h>    

#include "Bezier.cpp"
#include "Util.cpp"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;


void build_linspace(MatrixXd &linspace, const MatrixXd &V){
    for (size_t i = 0; i < linspace.rows(); i++){
      linspace(i, 0) = V.col(0).minCoeff() + ((V.col(0).maxCoeff() - V.col(0).minCoeff()) / (linspace.rows() - 1)) * i;
    }
}
// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier){
    Util curve_renderer(viewer); // 3D renderer for drawing points, polylines, curves, ...
    std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
    if (key == 'B'){
    }
    else if (key == 'S'){
      // subdivide
    }
    else if (key == 'U'){
      // compute tangents
    }
    else if (key == 'N'){
      // complute normals
    }
    else
      return false;

    //viewer.data(0).clear();
    //viewer.data(0).set_mesh(V1, F1);

    return false;
}


void test1(Util &curve_renderer){
    std::cout << "basic rendering of a Bezier curve of degree " << (V1.rows() - 1) << std::endl;
    Bezier bezier;
    int resolution = 500;
    std::cout<<bezier.de_casteljau(V1,0.76358)<<std::endl;
    MatrixXd plot = bezier.plot_curve(V1, resolution);
    curve_renderer.draw_curve(plot);
    curve_renderer.draw_control_polygon(V1);
}

void test2(Util &curve_renderer, double time_step){
    Bezier bezier;
    
    curve_renderer.draw_control_polygon(V1);
    MatrixXd V2 = V1;
    while(V2.rows()>1){
      V2 = bezier.de_casteljau_intermediate(V2,time_step);
      curve_renderer.draw_edges(V2);
    }

    //now draw on the same data channel the curve

    MatrixXd pts_for_plot = bezier.plot_curve(V1,100,time_step);

    curve_renderer.draw_edges(pts_for_plot);
}

void test3(Util &curve_renderer){
	std::cout << "splitting a Bezier curve of degree " << (V1.rows() - 1) << std::endl;
	Bezier bezier;
	int resolution = 500;
	RowVector3d green(0.1, 0.9, 0.1);
	RowVector3d yellow(1.0, 1.0, 0.1);

	vector<MatrixXd> curves = bezier.subdivide(V1, 0.5); // subdivide the curve for t=0.5
	MatrixXd b0 = curves[0];
	MatrixXd b1 = curves[1];

	MatrixXd p(1, 3), q(1, 3);
	p(0, 0) = 1.0; p(0, 1) = 0.0; p(0, 2) = 1.0;
	q(0, 0) = 1.0; q(0, 1) = 0.0; q(0, 2) = -0.8;

	MatrixXd plot = bezier.plot_curve(V1, resolution); //the input curve B(t)
	curve_renderer.draw_curve(plot);

	b0 = curve_renderer.translate_points(b0, p);
	MatrixXd plot0 = bezier.plot_curve(b0, resolution); // draw the first sub-curve
	curve_renderer.draw_colored_curve(plot0, green);
	curve_renderer.draw_control_polygon(b0);

	b1 = curve_renderer.translate_points(b1, q);
	MatrixXd plot1 = bezier.plot_curve(b1, resolution); // draw the second sub-curve
	curve_renderer.draw_colored_curve(plot1, yellow);
	curve_renderer.draw_control_polygon(b1);
	
}

void test4(Util &curve_renderer){
	std::cout << "rendering a Bezier curve with recursive subdivision" << std::endl;
	Bezier bezier;
	int levels = 2;

	MatrixXd plot = bezier.subdivision_plot(V1, levels);
	curve_renderer.draw_curve(plot);
	curve_renderer.draw_control_polygon(V1);
}

void test5(Util &curve_renderer){
	std::cout << "Computing tangents and normals " << std::endl;
	Bezier bezier;
	int resolution = 500;

	MatrixXd plot = bezier.plot_curve(V1, resolution);
	curve_renderer.draw_curve(plot);
	curve_renderer.draw_control_polygon(V1);

	for (double step = 0.1; step < 1.; step = step + 0.1){
		MatrixXd a = bezier.de_casteljau(V1, step);
		MatrixXd tangent = bezier.compute_tangent(V1, step);
		MatrixXd normal = bezier.compute_normal(V1, step);
		curve_renderer.draw_tangents(a, tangent * 4);
		curve_renderer.draw_normals(a, normal * 4);
	}
}

void test6(Util &curve_renderer){
	std::cout << "3D rendering of a Bezier curves " << std::endl;
	Bezier bezier;
	int resolution = 500; // for rendering the bezier curve
	int w=8; // width of the mesh
	int h=100; // height of the mesh
	double delta=1.0/h;

	MatrixXd plot = bezier.plot_curve(V1, resolution);
	curve_renderer.draw_curve(plot);
	//curve_renderer.draw_control_polygon(V1);

	MatrixXd loops[h]; // store all loop of vertices
	double t=0.0;
	for (int i=0; i<h; i++){
		MatrixXd loop = bezier.compute_loop_of_vertices(V1, t, w, 2.0*t);
		loops[i]=loop;
		t=t+delta;
		//curve_renderer.draw_close_curve(loop); // show only the loops of vertices
	}
	curve_renderer.draw_cylinder_grid(loops, h); // uncomment this line to render a thin cylinder
}



int main(int argc, char *argv[]){
	if(argc<2) {
		std::cout << "Error: input file required (.ply)" << std::endl;
		return 0;
	}
	std::cout << "reading input file: " << argv[1] << std::endl;
	igl::readPLY(argv[1], V1, F1);
	
	double time = 0.3;

	//  print the number of mesh elements
	std::cout << "Points: " << V1.rows() << std::endl;

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	Util curve_renderer(viewer);      // 3D renderer for drawing points, polylines, curves, ...


	

	//test1(curve_renderer);
	//test2(curve_renderer, time);
	//test3(curve_renderer);
	test4(curve_renderer);
	// test5(curve_renderer);
	// test6(curve_renderer);


	// //------------------- Uncomment this part for the animation -------------------
	// viewer.core().is_animating = true;  // Enable animation
	// double t = 0.0;  // Initialize time
	// double time_step = 0.01; // experiment with this parameter
		
	// viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool {
	// 	// Use a delay to control the speed of the color change
	// 	std::this_thread::sleep_for(std::chrono::milliseconds(150));
	// 	// Calculate new vertex colors based on your time-dependent function
	// 	viewer.data(0).clear();
	// 	test2(curve_renderer, t);
	// 	t += time_step;
	// 	if (t > 1.0) {
	// 		t = 0.0;  // Reset time to loop the animation
	// 	}
	// 	return false;
	// };

	viewer.core(0).align_camera_center(V1, F1);
	viewer.launch(); // run the viewer

}
