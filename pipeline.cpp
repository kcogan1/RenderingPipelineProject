#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <stdio.h>
#include <errno.h>
#include <vector>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <cfloat>
#define MAX DBL_MAX
using namespace std;

//STRUCTS
struct viewPoint {
	Eigen::Vector3d from;
	Eigen::Vector3d at;
	Eigen::Vector3d up;
	float angle;
	float hither;
	int resolution[2];
};

struct color {
	double fillColor[3];
	double Kd;
	double Ks;
	double shine;
	double T;
	double refraction;
};

struct polygon {
	Eigen::Matrix3d points;
	color colorInfo;
	Eigen::Matrix3d normals;
	Eigen::Matrix4d pixels;
	Eigen::Matrix3d colors;
};

//GLOBAL VARS
float bgColor[3];
viewPoint camera;
vector <polygon> worldShapes;
vector <color> colorVector;
Eigen::Vector3d light1;
Eigen::Vector3d light2;
unsigned char pixel[512][512][3]; //HARD CODED RESOLUTION
double zBuffer[512][512];

//FUNCTION DEFINITIONS
void parseFile();
void intersections(polygon*, double, double, double, double);
double shade(polygon*, int, int);
void processVertex(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
void raster();
bool blend(int, int, double);

void parseFile()
{
	//indicates what info is being parsed
	char indicator;
	char tempIndicator;
	//Char of arbitrary length to store the line read in
	char input[200];

	//temp vars to store the info being parsed
	char tempStr[20];
	int tempInt1;
	int tempInt2;
	float temp1;
	float temp2;
	float temp3;
	float temp4;
	float temp5;
	float temp6;
	float temp7;
	float temp8;

	int colorCounter = 0;

	//opens file
	FILE* file;
	file = fopen("teapot-3.nff", "r");
	//file = fopen("tetra-3.nff", "r");

	//if file is empty or does not exist
	if (file == NULL)
	{
		cout << "File is NULL" << endl;
		cout << errno << endl;
	}

	else
	{
		//gets first character of the line
		indicator = fgetc(file);

		//loops through file until end of file
		while (tempIndicator != EOF)
		{
			//BACKGROUND COLOR
			if (indicator == 'b')
			{
				//reads and parses line
				fgets(input, 50, file);
				sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

				//stores values
				bgColor[0] = temp1;
				bgColor[1] = temp2;
				bgColor[2] = temp3;

				//cout << "bgColor:" << temp1 << " " << temp2 << " " << temp3 << endl;
			}

			//VIEWPOINT
			else if (indicator == 'v')
			{
				//FROM
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.from(0) = temp1;
				camera.from(1) = temp2;
				camera.from(2) = temp3;
				//cout << "From: " << camera.from[0] << " " << camera.from[1] << " " << camera.from[2] << endl;

				//AT
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.at(0) = temp1;
				camera.at(1) = temp2;
				camera.at(2) = temp3;
				//cout << "at: " << camera.at[0] << " " << camera.at[1] << " " << camera.at[2] << endl;

				//UP
				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);
				camera.up(0) = temp1;
				camera.up(1) = temp2;
				camera.up(2) = temp3;
				//cout << "UP: " << camera.up[0] << " " << camera.up[1] << " " << camera.up[2] << endl;

				//ANGLE
				fgets(input, 50, file);
				sscanf(input, "%s %g", tempStr, &temp1);
				camera.angle = temp1;
				//cout << "Angle: " << camera.angle << endl;

				//HITHER
				fgets(input, 50, file);
				sscanf(input, "%s %g", tempStr, &temp1);
				camera.hither = temp1;
				//cout << "Hither: " << camera.hither << endl;

				//RESOLUTION
				fgets(input, 50, file);
				sscanf(input, "%s %d %d", tempStr, &tempInt1, &tempInt2);
				camera.resolution[0] = tempInt1;
				camera.resolution[1] = tempInt2;
				//cout << "resolution: " << camera.resolution[0] << " " << camera.resolution[1] << endl;
			}

			//LIGHT
			else if (indicator == 'l')
			{
				fgets(input, 50, file);
				sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

				light1[0] = temp1;
				light1[1] = temp2;
				light1[2] = temp3;

				//cout << "light1: " << light1[0] << " " << light1[1] << " " << light1[2] << endl;

				fgets(input, 50, file);
				sscanf(input, "%s %g %g %g", tempStr, &temp1, &temp2, &temp3);

				light2[0] = temp1;
				light2[1] = temp2;
				light2[2] = temp3;
				//cout << "light2: " << light2[0] << " " << light2[1] << " " << light2[2] << endl;
			}

			//FILL COLOR
			else if (indicator == 'f')
			{
				color* tempColor = new color;
				//reads and parses the line
				fgets(input, 50, file);
				sscanf(input, "%g %g %g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6, &temp7, &temp8);

				//stores values
				tempColor->fillColor[0] = temp1;
				tempColor->fillColor[1] = temp2;
				tempColor->fillColor[2] = temp3;

				tempColor->Kd = temp4;
				tempColor->Ks = temp5;
				tempColor->shine = temp6;
				tempColor->T = temp7;
				tempColor->refraction = temp8;

				colorVector.push_back(*tempColor);
				colorCounter++;
				
				//cout << "Fill Color:" << fillColor[0] << " " << fillColor[1] << " " << fillColor[2] << endl;
			}

			//POLYGON
			else if (indicator == 'p')
			{
				fgets(input, 200, file);
				sscanf(input, "%d", &tempInt1);

				//if the nff give the normals as well
				if (tempIndicator == 'p')
				{
					//if the polygon has 3 verticies
					if (tempInt1 == 3)
					{
						//creates new polygon object
						polygon* tempPoly = new polygon;

						//loops through the 3 verticies and stores the info
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);

							tempPoly->points(i, 0) = temp1;
							tempPoly->points(i, 1) = temp2;
							tempPoly->points(i, 2) = temp3;

							tempPoly->normals(i, 0) = temp4;
							tempPoly->normals(i, 1) = temp5;
							tempPoly->normals(i, 2) = temp6;
						}

						//adds the new polygon to a vector called "shapes"
						tempPoly->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;
						worldShapes.push_back(*tempPoly);
					}

					//if the polygon has 4 verticies
					else if (tempInt1 == 4)
					{
						//creates new polygon object
						polygon* tempPoly1 = new polygon;
						polygon* tempPoly2 = new polygon;

						int counter = 0;
						for (int i = 0; i < 3; i++)
						{
							fgets(input, 200, file);
							sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);

							tempPoly1->points(i, 0) = temp1;
							tempPoly1->points(i, 1) = temp2;
							tempPoly1->points(i, 2) = temp3;

							tempPoly1->normals(i, 0) = temp4;
							tempPoly1->normals(i, 1) = temp5;
							tempPoly1->normals(i, 2) = temp6;

							if (i == 0 || i == 2 || i == 3)
							{
								tempPoly2->points(counter, 0) = temp1;
								tempPoly2->points(counter, 1) = temp2;
								tempPoly2->points(counter, 2) = temp3;

								tempPoly2->normals(counter, 0) = temp4;
								tempPoly2->normals(counter, 1) = temp5;
								tempPoly2->normals(counter, 2) = temp6;

								counter++;
							}
						}

						//after the loop, it has to get the fourth vertex becuase the for loop will never hit it.
						fgets(input, 200, file);
						sscanf(input, "%g %g %g %g %g %g", &temp1, &temp2, &temp3, &temp4, &temp5, &temp6);

						tempPoly2->points(2, 0) = temp1;
						tempPoly2->points(2, 1) = temp2;
						tempPoly2->points(2, 2) = temp3;

						tempPoly2->normals(2, 0) = temp4;
						tempPoly2->normals(2, 1) = temp5;
						tempPoly2->normals(2, 2) = temp6;

						//adds the new polygon to a vector called "shapes"
						tempPoly1->colorInfo = colorVector[(colorCounter - 1)];
						tempPoly2->colorInfo = colorVector[(colorCounter - 1)];
						//cout << "color counter = " << colorCounter << endl;
						worldShapes.push_back(*tempPoly1);
						worldShapes.push_back(*tempPoly2);
					}
				}

				//if the polygon has 3 verticies
				else if (tempInt1 == 3)
				{
					//creates new polygon object
					polygon* tempPoly = new polygon;

					//loops through the 3 verticies and stores the info
					for (int i = 0; i < 3; i++)
					{
						fgets(input, 200, file);
						sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

						tempPoly->points(i, 0) = temp1;
						tempPoly->points(i, 1) = temp2;
						tempPoly->points(i, 2) = temp3;
					}
					//adds the new polygon to a vector called "shapes"
					tempPoly->colorInfo = colorVector[(colorCounter - 1)];
					//cout << "color counter = " << colorCounter << endl;
					tempPoly->normals(0, 0) = -1;
					tempPoly->normals(0, 1) = -1;
					tempPoly->normals(0, 2) = -1;
					worldShapes.push_back(*tempPoly);
				}

				//if the polygon has 4 verticies
				else if (tempInt1 == 4)
				{
					//creates new polygon object
					polygon* tempPoly1 = new polygon;
					polygon* tempPoly2 = new polygon;

					//loops through the 3 verticies and stores the info
					int counter = 0;
					for (int i = 0; i < 3; i++)
					{
						fgets(input, 200, file);
						sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

						tempPoly1->points(i, 0) = temp1;
						tempPoly1->points(i, 1) = temp2;
						tempPoly1->points(i, 2) = temp3;

						if (i == 0 || i == 2 || i == 3)
						{
							tempPoly2->points(counter, 0) = temp1;
							tempPoly2->points(counter, 1) = temp2;
							tempPoly2->points(counter, 2) = temp3;
							counter++;
						}
					}

					//after the loop, it has to get the fourth vertex becuase the for loop will never hit it.
					fgets(input, 200, file);
					sscanf(input, "%g %g %g", &temp1, &temp2, &temp3);

					tempPoly2->points(2, 0) = temp1;
					tempPoly2->points(2, 1) = temp2;
					tempPoly2->points(2, 2) = temp3;

//adds the new polygon to a vector called "shapes"
tempPoly1->colorInfo = colorVector[(colorCounter - 1)];
tempPoly2->colorInfo = colorVector[(colorCounter - 1)];
//cout << "color counter = " << colorCounter << endl;

tempPoly1->normals(0, 0) = -1;
tempPoly2->normals(0, 0) = -1;
tempPoly1->normals(0, 1) = -1;
tempPoly2->normals(0, 1) = -1;
tempPoly1->normals(0, 2) = -1;
tempPoly2->normals(0, 2) = -1;
worldShapes.push_back(*tempPoly1);
worldShapes.push_back(*tempPoly2);
				}
			}

			//handles lines in the NFF file that do not have an indicator letter
			else
			{
			fgets(input, 50, file);
			}

			//gets new indicator for the next iteration of the loop
			indicator = fgetc(file);
			tempIndicator = fgetc(file);
		}
	}

	if (file != NULL)
		fclose(file);
}

void intersections(polygon* polyInfo, double xMin, double xMax, double yMin, double yMax)
{
	//initialize numbers and x and y values
	double f01;
	double f12;
	double f20;
	double falpha;
	double fbeta;
	double fgamma;

	double fNeg12;
	double fNeg20;
	double fNeg01;

	double alpha;
	double beta;
	double gamma;

	double x0 = polyInfo->pixels(0, 0);
	double x1 = polyInfo->pixels(1, 0);
	double x2 = polyInfo->pixels(2, 0);
	double y0 = polyInfo->pixels(0, 1);
	double y1 = polyInfo->pixels(1, 1);
	double y2 = polyInfo->pixels(2, 1);

	//used for blender
	bool illuminate1;
	bool illuminate2;
	bool illuminate3;

	int xLoc;
	int yLoc;

	//clamp values between 0 and resolution
	if (xMin < 0)
		xMin = 0;
	if (xMax > camera.resolution[0])
		xMax = camera.resolution[0];
	if (yMin < 0)
		yMin = 0;
	if (yMax > camera.resolution[1])
		yMax = camera.resolution[1];

	//loop over bounding box for triangle
	for(int y = yMin; y < yMax; y++)
	{
		for(int x = xMin; x < xMax; x++)
		{
			//compute alpha, beta, and gamma
			f12 = (y1 - y2) * x + (x2 - x1) * y + (x1 * y2) - (x2 * y1);
			falpha = (y1 - y2) * x0 + (x2 - x1) * y0 + (x1 * y2) - (x2 * y1);
			alpha = f12 / falpha;
			fNeg12 = (y1 - y2) * -1 + (x2 - x1) * -1 + (x1 * y2) - (x2 * y1);

			f20 = (y2 - y0) * x + (x0 - x2) * y + (x2 * y0) - (x0 * y2);
			fbeta = (y2 - y0) * x1 + (x0 - x2) * y1 + (x2 * y0) - (x0 * y2);
			beta = f20 / fbeta;
			fNeg20 = (y2 - y0) * -1 + (x0 - x2) * -1 + (x2 * y0) - (x0 * y2);

			f01 = (y0 - y1) * x + (x1 - x0) * y + (x0 * y1) - (x1 * y0);
			fgamma = (y0 - y1) * x2 + (x1 - x0) * y2 + (x0 * y1) - (x1 * y0);
			gamma = f01 / fgamma;
			fNeg01 = (y0 - y1) * -1 + (x1 - x0) * -1 + (x0 * y1) - (x1 * y0);

			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				//handes triangle edges on pixels
				if ((alpha > 0 || falpha * fNeg12 > 0) && (beta > 0 || fbeta * fNeg20 > 0) && (gamma > 0 || fgamma * fNeg01 > 0))
				{
					//z buffer- blend using all 3 z values of the current trianlge to reduce black spots
					illuminate1 = blend(y, x, polyInfo->pixels(0, 2));
					illuminate2 = blend(y, x, polyInfo->pixels(1, 2));
					illuminate3 = blend(y, x, polyInfo->pixels(2, 2));

					if (illuminate1 == true || illuminate2 == true || illuminate3 == true)
					{
						//then, calculate shading for all verts
						//color at v0 * alpha
						double r1 = polyInfo->colors(0, 0) * alpha;
						double g1 = polyInfo->colors(0, 1) * alpha;
						double b1 = polyInfo->colors(0, 2) * alpha;

						//color at v1 * beta
						double r2 = polyInfo->colors(1, 0) * beta;
						double g2 = polyInfo->colors(1, 1) * beta;
						double b2 = polyInfo->colors(1, 2) * beta;

						//color at v2 * gamma
						double r3 = polyInfo->colors(2, 0) * gamma;
						double g3 = polyInfo->colors(2, 1) * gamma;
						double b3 = polyInfo->colors(2, 2) * gamma;

						double r = r1 + r2 + r3;
						double g = g1 + g2 + g3;
						double b = b1 + b2 + b3;

						//clamp values between 1 and 0
						if (r > 1)
							r = 1;
						if (g > 1)
							g = 1;
						if (b > 1)
							b = 1;

						//drae pixels
						pixel[y][x][0] = r * 255;
						pixel[y][x][1] = g * 255;
						pixel[y][x][2] = b * 255;
					}
				}

					//if (illuminate == false) just dont illuminate becuase the current fragment is not closer than past fragments.	
			}
		}
	}
}

double shade(polygon* polyInfo, int index, int vert)
{
	Eigen::Vector3d normal;
	double localColor;
	double diffuse;
	double specular;
	double lightIntensity = (1.0 / sqrt(2));
	Eigen::Vector3d h;
	Eigen::Vector3d location;

	//if the point is not a pp, calculate normal
	if (polyInfo->normals(0, 0) == -1 && polyInfo->normals(0, 1) == -1 && polyInfo->normals(0,2) == -1)
	{
		//cout << "not a polygon patch, calculating normal" << endl;

		//calculate surface normal of current polygon
		Eigen::Vector3d U; //v1 - v0
		Eigen::Vector3d V; //v2 - v0
		U(0) = polyInfo->points(1, 0) - polyInfo->points(0, 0);
		U(1) = polyInfo->points(1, 1) - polyInfo->points(0, 1);
		U(2) = polyInfo->points(1, 2) - polyInfo->points(0, 2);

		V(0) = polyInfo->points(2, 0) - polyInfo->points(0, 0);
		V(1) = polyInfo->points(2, 1) - polyInfo->points(0, 1);
		V(2) = polyInfo->points(2, 2) - polyInfo->points(0, 2);

		//normal is a vector
		normal = U.cross(V);
		normal.normalize();
		
		location(0) = polyInfo->points(vert, 0);
		location(1) = polyInfo->points(vert, 1);
		location(2) = polyInfo->points(vert, 2);

		Eigen::Vector3d l = light1 - location;
		Eigen::Vector3d direction = location - camera.from;
		//instead of using intersect locations since you dont have that, just shade each vertex and use that as the intersect location.
		l.normalize();

		//do for first light----------------------------------------------------------------------
		h = (direction + l) / (direction + l).norm();
		diffuse = max(0.0, (normal.dot(l)));
		specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);
		localColor = ((polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse) + polyInfo->colorInfo.Ks * specular) * lightIntensity;

		//do for second light----------------------------------------------------------------------
		l = light2 - location;
		l.normalize();
		h = (direction + l) / (direction + l).norm();
		diffuse = max(0.0, (normal.dot(l)));
		specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);
		localColor += ((polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse) + polyInfo->colorInfo.Ks * specular) * lightIntensity;
	}	

	//else, use given normal for the current vertex
	else
	{
		normal(0) = polyInfo->normals(vert, 0);
		normal(1) = polyInfo->normals(vert, 1);
		normal(2) = polyInfo->normals(vert, 2);
		normal.normalize();

		location(0) = polyInfo->points(vert, 0);
		location(1) = polyInfo->points(vert, 1);
		location(2) = polyInfo->points(vert, 2);

		Eigen::Vector3d l = light1 - location;
		Eigen::Vector3d direction = location - camera.from; //"ray" direction is just vertex - camera
		l.normalize();

		//do for first light----------------------------------------------------------------------
		h = (direction + l) / (direction + l).norm();
		diffuse = max(0.0, (normal.dot(l)));
		specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);
		localColor = ((polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse) + polyInfo->colorInfo.Ks * specular) * lightIntensity;

		//do for second light----------------------------------------------------------------------
		l = light2 - location;
		l.normalize();
		h = (direction + l) / (direction + l).norm();
		diffuse = max(0.0, (normal.dot(l)));
		specular = pow(max(0.0, normal.dot(h)), polyInfo->colorInfo.shine);
		localColor += ((polyInfo->colorInfo.Kd * polyInfo->colorInfo.fillColor[index] * diffuse) + polyInfo->colorInfo.Ks * specular) * lightIntensity;
	}

	//checks that the color is between 1 and 0
	if (localColor > 1.0)
		localColor = 1;
	
	if (localColor < 0.0)
		localColor = 0;

	return localColor;
	//must call this funtion for each vertex, and for each vertex must call it for r, g, and b.
}

void processVertex(Eigen::Vector3d u, Eigen::Vector3d v, Eigen::Vector3d w)
{
	//Initialize Matrix Mvp
	Eigen::Matrix4d Mvp;
	//row 1
	Mvp(0, 0) = camera.resolution[0] / 2.0;
	Mvp(0, 1) = 0;
	Mvp(0, 2) = 0;
	Mvp(0, 3) = (camera.resolution[0] - 1) / 2.0;
	//row 2
	Mvp(1, 0) = 0;
	Mvp(1, 1) = camera.resolution[1] / 2.0;
	Mvp(1, 2) = 0;
	Mvp(1, 3) = (camera.resolution[1] - 1) / 2.0;
	//row 3
	Mvp(2, 0) = 0;
	Mvp(2, 1) = 0;
	Mvp(2, 2) = 1;
	Mvp(2, 3) = 0;
	//row 4
	Mvp(3, 0) = 0;
	Mvp(3, 1) = 0;
	Mvp(3, 2) = 0;
	Mvp(3, 3) = 1;
	cout << "Mvp:" << endl << Mvp << endl << endl;

	//Initialize Matrix P and Morth
	Eigen::Matrix4d P;
	Eigen::Matrix4d Morth;
	double n =  -1.0 * camera.hither;
	double f = 1000 * n;

	//variables needed to calculate Morth
	double rads = (M_PI / 180.0) * (camera.angle);
	double t = -1 * n * tan(rads/2.0);
	double b = -1 * t;
	double r = (camera.resolution[0] / camera.resolution[1]) * t;
	double l = -1 * r;

	//row 1
	P(0, 0) = n;
	P(0, 1) = 0;
	P(0, 2) = 0;
	P(0, 3) = 0;
	//row 2
	P(1, 0) = 0;
	P(1, 1) = n;
	P(1, 2) = 0;
	P(1, 3) = 0;
	//row 3
	P(2, 0) = 0;
	P(2, 1) = 0;
	P(2, 2) = n + f;
	P(2, 3) = -1.0 * (f * n);
	//row 4
	P(3, 0) = 0;
	P(3, 1) = 0;
	P(3, 2) = 1;
	P(3, 3) = 0;

	//row 1
	Morth(0, 0) = 2.0/(r-l);
	Morth(0, 1) = 0;
	Morth(0, 2) = 0;
	Morth(0, 3) = -1.0 * ((r + l) / (r - l));
	//row 2
	Morth(1, 0) = 0;
	Morth(1, 1) = 2.0 / (t - b);
	Morth(1, 2) = 0;
	Morth(1, 3) = -1.0 * ((t + b) / (t - b));
	//row 3
	Morth(2, 0) = 0;
	Morth(2, 1) = 0;
	Morth(2, 2) = 2 / (n - f);
	Morth(2, 3) = -1.0 * ((n + f) / (n - f));
	//row 4
	Morth(3, 0) = 0;
	Morth(3, 1) = 0;
	Morth(3, 2) = 0;
	Morth(3, 3) = 1;
	//calculate Mperspective
	Eigen::Matrix4d Mper = Morth * P; //<-do not swap order of matrices in multiplication!!

	cout << "Mper:" << endl << Mper << endl << endl;

	Eigen::Matrix4d Mcam;
	Eigen::Matrix4d uvw;
	Eigen::Matrix4d e;
	//row 1
	uvw(0, 0) = u(0);
	uvw(0, 1) = u(1);
	uvw(0, 2) = u(2);
	uvw(0, 3) = 0;
	//row 2
	uvw(1, 0) = -1 * v(0);
	uvw(1, 1) = -1 * v(1);
	uvw(1, 2) = -1 * v(2);
	uvw(1, 3) = 0;
	//row 3
	uvw(2, 0) = w(0);
	uvw(2, 1) = w(1);
	uvw(2, 2) = w(2);
	uvw(2, 3) = 0;
	//row 4
	uvw(3, 0) = 0;
	uvw(3, 1) = 0;
	uvw(3, 2) = 0;
	uvw(3, 3) = 1;

	//row1
	e(0, 0) = 1;
	e(0, 1) = 0;
	e(0, 2) = 0;
	e(0, 3) = camera.from[0] * -1;
	//row 2
	e(1, 0) = 0;
	e(1, 1) = 1;
	e(1, 2) = 0;
	e(1, 3) = camera.from[1] * -1;
	//row 3
	e(2, 0) = 0;
	e(2, 1) = 0;
	e(2, 2) = 1;
	e(2, 3) = camera.from[2] * -1;
	//row 4
	e(3, 0) = 0;
	e(3, 1) = 0;
	e(3, 2) = 0;
	e(3, 3) = 1;
	Mcam = uvw * e;
	cout << "Mcam:" << endl << Mcam << endl << endl;

	//multiply all matrices together to get Matrix M
	Eigen::Matrix4d M = Mvp * Mper * Mcam;
	cout << "M:" << endl << M << endl;
	
	//initialize some vectors to do math
	Eigen::Vector4d temp;
	Eigen::Vector4d answer;
	temp(3) = 1;

	//loops through all polygons
	for (int k = 0; k < worldShapes.size(); k++)
	{
		//loops through all vertices of the current polygon
		for (int i = 0; i < 3; i++){ //row
			for (int j = 0; j < 3; j++) //column
			{
				//this will fill in xyz of temp, one row of points
				temp(j) = worldShapes[k].points(i, j);
			}

			//once temp is filled in, do the multiplication.
			answer = M * temp;

			//Homogeneous divide
			answer /= answer(3);

			//fills in matrix "pixels" with answer
			for (int b = 0; b < 4; b++)
			{
				worldShapes[k].pixels(i, b) = answer(b);
			}

			//then, calculate shading for all verts usting world coordinates
			double r = shade(&worldShapes[k], 0, i); //i is the vertex the loop is on currently
			double g = shade(&worldShapes[k], 1, i);
			double b = shade(&worldShapes[k], 2, i);

			//store new color
			worldShapes[k].colors(i, 0) = r;
			worldShapes[k].colors(i, 1) = g;
			worldShapes[k].colors(i, 2) = b;
		}
	}
}

void raster()
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	double xTemp;
	double yTemp;


	//loop over all polygons and computer the bounding box for each polygon
	for (int k = 0; k < worldShapes.size(); k++)
	{
		//xmin
		xTemp = worldShapes[k].pixels(0, 0);
		for (int x = 1; x < 3; x++)
		{
			if (worldShapes[k].pixels(x, 0) < xTemp)
				xTemp = worldShapes[k].pixels(x, 0);
		}
		xMin = floor(xTemp);

		//ymin
		yTemp = worldShapes[k].pixels(0, 1);
		for (int y = 1; y < 3; y++)
		{
			if (worldShapes[k].pixels(y, 1) < yTemp)
				yTemp = worldShapes[k].pixels(y, 1);
		}
		yMin = floor(yTemp);

		//xmax
		xTemp = worldShapes[k].pixels(0, 0);
		for (int x = 1; x < 3; x++)
		{
			if (worldShapes[k].pixels(x, 0) > xTemp)
				xTemp = worldShapes[k].pixels(x, 0);
		}
		xMax = ceil(xTemp);

		//ymax
		yTemp = worldShapes[k].pixels(0, 1);
		for (int y = 1; y < 3; y++)
		{
			if (worldShapes[k].pixels(y, 1) > yTemp)
				yTemp = worldShapes[k].pixels(y, 1);
		}
		yMax = ceil(yTemp);

		//find which pixels intersect the triangle
		intersections(&worldShapes[k], xMin, xMax, yMin, yMax);	
	}
}

bool blend(int y, int x, double currZ)
{
	//if pixel has not yet been illuminated, it is the first fragment in that pixel
	if (zBuffer[y][x] == NULL)
	{
		zBuffer[y][x] = currZ;
		//return true if we illluminate this pixel now
		return true;
	}

	//if this pixel has already been illuminated by another fragment, compare their z values
	else if (currZ > zBuffer[y][x])
	{
		zBuffer[y][x] = currZ;
		//return true if we illluminate this pixel now
		return true;
	}

	//return false if this fragment is farther than the fragment already in that pixel
	else
		return false;
}

int main()
{
	parseFile();

	//initialize pixels as 2 to keep track of unvisited pixels, which will later be the background color
	for (int y = 0; y < camera.resolution[1]; y++){
		for (int x = 0; x < camera.resolution[0]; x++){
			pixel[y][x][0] = 2;
		}
	}

	//w
	Eigen::Vector3d w;
	w = camera.from - camera.at;
	w /= w.norm();
	//u
	//made u negative in processVertex function to flip the image right side up
	Eigen::Vector3d u;
	u = (camera.up).cross(w);
	u.normalize();
	//v
	Eigen::Vector3d v;
	v = w.cross(u);

	//process vertex and raster, raster calls the intersect function
	processVertex(u, v, w);
	raster();

	//fill in background color, loop through all pixels
	for (int y = 0; y < camera.resolution[1]; y++){
		for (int x = 0; x < camera.resolution[0]; x++){

			//if the pixel is unvisited by the vertex and raster functions (using indicator of 2), make it the background color
			if (pixel[y][x][0] == 2)
			{
				pixel[y][x][0] = bgColor[0] * 255;
				pixel[y][x][1] = bgColor[1] * 255;
				pixel[y][x][2] = bgColor[2] * 255;
			}
		}
	}

	//writes pixels to ppm file
	FILE* k = fopen("C:\\Users\\bluej\\OneDrive\\Documents\\UMBC\\SENIOR+\\435\\proj1435\\hide.ppm", "wb");
	fprintf(k, "P6\n%d %d\n%d\n", 512, 512, 255); //HARD CODED RESOLIUTION
	fwrite(pixel, 1, 512 * 512 * 3, k);
	fclose(k);

	return 0;
}