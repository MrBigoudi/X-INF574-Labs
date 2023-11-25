/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class SphereGeneration
{

public:
	/** 
	 * Initialize the data structures
	 **/
	SphereGeneration(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
	{
		he = &mesh;
		V = &V_original;
		F = &F_original; // NOT NEEDED if using the half-edge data structure
		int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
		int n = V_original.rows();		   // number of vertices in the original mesh

		// TO BE COMPLETED (initialize arrays V1 and F1)
		V1 = new MatrixXd(MatrixXd::Zero(2*n, V->cols()));
		F1 = new MatrixXi(MatrixXi::Zero(4*F->rows(), F->cols()));
	}

	/** 
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
	void subdivide()
	{
		std::cout << "Performing one round subdivision" << endl;
		int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
		int n = he->sizeOfVertices();	  // number of vertices in the original mesh
		int f = he->sizeOfFaces();		   // number of faces in the original mesh

		// TO BE COMPLETED
		// an array containing the index of the midpoint vertex for each half edge
		int halfEdgeMidPointsIdx[2*e];
		for(int i=0; i<2*e; i++){
			halfEdgeMidPointsIdx[i] = -1;
		}
		// first step: perform a copy of the original points
		HalfedgeDS* newHe = new HalfedgeDS(2*n, 2*(3*e), 4*f);
		for(int i=0; i<n; i++){
			V1->row(i) = V->row(i);
		}

		// second step: compute new midpoint vertices and assign a number, between 0..e-1, to all halfedges
		for(int i=0; i<e; i++){
			if(halfEdgeMidPointsIdx[i] == -1){
				MatrixXd midPoint = computeEdgePoint(i);
				int newIdx = n+1;
				V1->row(newIdx) = midPoint;
				halfEdgeMidPointsIdx[i] = newIdx;
				halfEdgeMidPointsIdx[he->getOpposite(i)] = newIdx;
			}
		}


		// third step: set the face/vertex incidence relations
		int edgeCounter = 0;
		for(int i=0; i<f; i++){
			// get the current face
			int curFace = i;

			// get the 3 vertices
			int h01 = he->getEdgeInFace(curFace); // edge from v0 to v1
			int h12 = he->getNext(h01);
			int h20 = he->getNext(h12);
			int v0 = he->getTarget(h20);
			int v1 = he->getTarget(h01);
			int v2 = he->getTarget(h12);

			// get the 3 new points
			int mid01 = halfEdgeMidPointsIdx[h01];
			int mid12 = halfEdgeMidPointsIdx[h12];
			int mid20 = halfEdgeMidPointsIdx[h20];

			// get the 9 new edges
			int h_0_mid01, h_mid01_0 = 0;
			int h_mid01_1, h_1_mid01 = 0;
			int h_1_mid12, h_mid12_1 = 0;
			int h_mid12_2, h_2_mid12 = 0;
			int h_2_mid20, h_mid20_2 = 0;
			int h_mid20_0, h_0_mid20 = 0;
			// pair h_0_mid01, h_mid01_0
			newHe->setNewEdge(mid01, h_0_mid01, h_mid01_0, edgeCounter);
			// pair h_1_mid01, h_mid01_1
			newHe->setNewEdge(v1, h_mid01_1, h_1_mid01, edgeCounter);
			// pair h_1_mid12, h_mid12_1
			newHe->setNewEdge(mid12, h_1_mid12, h_mid12_1, edgeCounter);
			// pair h_2_mid12, h_mid12_2
			newHe->setNewEdge(v2, h_mid12_2, h_2_mid12, edgeCounter);
			// pair h_2_mid20, h_mid20_2
			newHe->setNewEdge(mid20, h_2_mid20, h_mid20_2, edgeCounter);
			// pair h_0_mid20, h_mid20_0
			newHe->setNewEdge(v0, h_mid20_0, h_0_mid20, edgeCounter);
			
			int h_mid01_mid12 = edgeCounter;
			int h_mid12_mid01 = edgeCounter+1;
			newHe->setOpposite(h_mid01_mid12, h_mid12_mid01);
			newHe->setOpposite(h_mid12_mid01, h_mid01_mid12);
			int h_mid12_mid20 = edgeCounter+2;
			int h_mid20_mid12 = edgeCounter+3;
			newHe->setOpposite(h_mid12_mid20, h_mid20_mid12);
			newHe->setOpposite(h_mid20_mid12, h_mid12_mid20);
			int h_mid20_mid01 = edgeCounter+4;
			int h_mid01_mid20 = edgeCounter+5;
			newHe->setOpposite(h_mid20_mid01, h_mid01_mid20);
			newHe->setOpposite(h_mid01_mid20, h_mid20_mid01);
			edgeCounter+= 6;


			// create the relations
			// set the nexts
		}

	}

	/** 
	 * Return the number of half-edges
	 **/
	MatrixXd getVertexCoordinates()
	{
		return *V1;
	}

	/** 
	 * Return the number of faces
	 **/
	MatrixXi getFaces()
	{
		return *F1;
	}

	/** 
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
	void print(int verbosity)
	{
		cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

		if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
		{
			for (int i = 0; i < nVertices; i++)
			{
				cout << "v" << i << ": " << V1->row(i) << endl;
			}

			std::cout << "new faces: " << nFaces << endl;
			for (int i = 0; i < nFaces; i++)
			{
				cout << "f" << i << ": " << F1->row(i) << endl;
			}
		}
	}

private:
	/**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
	MatrixXd computeEdgePoint(int h)
	{
		// get u, v
		int uIdx = he->getTarget(h);
		int vIdx = (he->getTarget(he->getOpposite(h)));
		MatrixXd u = V->row(uIdx);
		MatrixXd v = V->row(vIdx);
		// get half
		MatrixXd edgePoint = (u+v) / 2.0f;
		// reproject
		edgePoint.normalize();
		return edgePoint;
	}

	/** Half-edge representation of the original input mesh */
	HalfedgeDS *he;
	MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure

	int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
	MatrixXd *V1;		   // vertex coordinates of the new subdivided mesh
	MatrixXi *F1;		   // faces of the new subdivided mesh
};
