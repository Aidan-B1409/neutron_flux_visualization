#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

Polyhedron* poly;
std::vector<PolyLine> lines;
std::vector<icVector3> points;
std::vector<_RGB> colors;

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		64
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/normalfuelgeometry.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}
	
	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->ntris; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Triangle* temp_q = this_poly->tlist[i];
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 3; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_tri == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Triangle* temp_q = this_poly->tlist[this_poly->selected_tri];

	glBegin(GL_POLYGON);
	for (j = 0; j < 3; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0,0.0);

	CHECK_GL_ERROR();
}

// credit - https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
_RGB hsv2rgb(HSV in)
{
	double      hh, p, q, t, ff;
	long        i;
	_RGB        out;

	if (in.s <= 0.0) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0 - in.s);
	q = in.v * (1.0 - (in.s * ff));
	t = in.v * (1.0 - (in.s * (1.0 - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}

void classify_vertices(double s) {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vtmp = poly->vlist[i];
		if (vtmp->s >= s) {
			vtmp->level = 1;
		}
		else {
			vtmp->level = 0;
		}
	}
}

void classify_and_raise_vertices(double s, Vertex* v_min, Vertex* v_max) {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vtmp = poly->vlist[i];
		double vscale = (vtmp->s - v_min->s) / (v_max->s - v_min->s);
		if (vtmp->s >= s) {
			vtmp->level = 1;
		}
		else {
			vtmp->level = 0;
		}
		vtmp->z = vscale;
	}
}


void calculate_intersections(double s) {
	for (int i = 0; i < poly->nedges; i++) {
		Edge* e_tmp = poly->elist[i];
		Vertex* v1 = e_tmp->verts[0];
		Vertex* v2 = e_tmp->verts[1];

		if (e_tmp->cross_point) {
			delete e_tmp->cross_point;
			e_tmp->cross_point = NULL;
		}

		if (v1->level != v2->level) {
			double alpha = (s - v1->s) / (v2->s - v1->s);
			double cross_x = (alpha * (v2->x - v1->x)) + v1->x;
			double cross_y = (alpha * (v2->y - v1->y)) + v1->y;
			double cross_z = (alpha * (v2->z - v1->z)) + v1->z;

			// TODO - free this
			icVector3* cross = new icVector3(cross_x, cross_y, cross_z);
			e_tmp->cross_point = cross;

		}
	}
}

_RGB get_bicolor_mapping(double s, Vertex* v_min, Vertex* v_max) {
	double r = (1 - 0) * ((s - v_min->s) / (v_max->s - v_min->s)) + 0;
	double g = 0;
	double b = (1 - 0) * ((v_max->s - s) / (v_max->s - v_min->s)) + 0;
	_RGB rgb = { r, g, b };
	return rgb;
}

PolyLine connect_crossings(double s) {
	PolyLine contour; 
	for (int i = 0; i < poly->ntris; i++) {
		Triangle* q_tmp = poly->tlist[i];

		std::vector<icVector3*> crossingPoints;

		for (int j = 0; j < 3; j++) {
			Edge* e_tmp = q_tmp->edges[j];
			if (e_tmp->cross_point) {
				crossingPoints.push_back(e_tmp->cross_point);
			}
		}

		switch (crossingPoints.size()) {
		case 2:
		{
			LineSegment line(*crossingPoints[0], *crossingPoints[1]);
			contour.push_back(line);

		}
		break;
		case 4:
		{

		}
		break;
		}
	}
	return contour;
}

Vertex* lookup_vertex_from_coords(double x, double y, Vertex* vertlist[]) {
	for (int i = 0; i < 3; i++) {
		Vertex* tmp = vertlist[i];
		if (tmp->x == x && tmp->y == y) {
			return tmp;
		}
	}
	return &Vertex(0, 0, 0);
}

void find_critical_points(Vertex* v_min, Vertex* v_max) {
	for (int i = 0; i < poly->ntris; i++) {
		Triangle* q_tmp = poly->tlist[i];
		double x1 = q_tmp->verts[0]->x,
			x2 = q_tmp->verts[0]->x,
			y1 = q_tmp->verts[0]->y,
			y2 = q_tmp->verts[0]->y,
			fx1y1 = 0,
			fx2y1 = 0,
			fx1y2 = 0,
			fx2y2 = 0,
			x0 = 0,
			y0 = 0,
			vscale = 0;

		for (int j = 0; j < 4; j++) {
			Vertex* v_tmp = q_tmp->verts[j];
			if (v_tmp->x < x1) {
				x1 = v_tmp->x;
			}
			if (v_tmp->x > x2) {
				x2 = v_tmp->x;
			}
			if (v_tmp->y < y1) {
				y1 = v_tmp->y;
			}
			if (v_tmp->y > y2) {
				y2 = v_tmp->y;
			}
		}

		fx1y1 = lookup_vertex_from_coords(x1, y1, q_tmp->verts)->s;
		fx2y1 = lookup_vertex_from_coords(x2, y1, q_tmp->verts)->s;
		fx1y2 = lookup_vertex_from_coords(x1, y2, q_tmp->verts)->s;
		fx2y2 = lookup_vertex_from_coords(x2, y2, q_tmp->verts)->s;


		x0 = ((x2 * fx1y1) - (x1 * fx2y1) - (x2 * fx1y2) + (x1 * fx2y2)) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);
		y0 = ((y1 * fx1y1) - (y2 * fx2y1) - (y1 * fx1y2) + (y1 * fx2y2)) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);


		// Dumb and bad approximation
		vscale = (fx1y1 - v_min->s) / (v_max->s - v_min->s);

		if (((x0 > x1) && (x0 < x2)) && ((y0 > y1) && (y0 < y2))) {
			icVector3* cp = new icVector3(x0, y0, vscale);
			q_tmp->critical = cp;
			points.push_back(*cp);

		}
	}
}

double hue_mapping(double s, Vertex* v_min, Vertex* v_max) {
	return (((s - v_min->s) / (v_max->s - v_min->s)) * 360);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	// clear out lines and points
	lines.clear();
	points.clear();

	switch (key) {
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':	// checkerboard display
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];
			for (int j = 0; j < 3; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;

	case '4':	// Drawing points and lines created by the dots_and_lines_example() function
		display_mode = 4;
		dots_and_lines_example(&points, &lines);
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		glutPostRedisplay();
		break;

	// bicolor contour
	case '6':
	{
		display_mode = 6;

		Vertex* v_max = poly->max_vertex();
		Vertex* v_min = poly->min_vertex();

		double n = 20;
		double n_dist = (v_max->s - v_min->s) / (double)(n);
		for (float i = v_min->s; i < v_max->s; i += n_dist) {
			classify_vertices(i);
			calculate_intersections(i);
			lines.push_back(connect_crossings(i));

		}

		for (float i = v_min->s; i < v_max->s; i += n_dist * 2) {
			colors.push_back(get_bicolor_mapping(i, v_min, v_max));
		}
		int len = colors.size();
		std::vector<_RGB> temp_color = colors;
		std::reverse(temp_color.begin(), temp_color.end());
		colors.insert(colors.end(), temp_color.begin(), temp_color.end());
	}
	glutPostRedisplay();
	break;

	// rainbow contour
	case '7':
	{
		display_mode = 7;

		Vertex* v_max = poly->max_vertex();
		Vertex* v_min = poly->min_vertex();

		double n = 20;
		double n_dist = (v_max->s - v_min->s) / (double)(n);
		for (float i = v_min->s; i < v_max->s; i += n_dist) {
			classify_vertices(i);
			calculate_intersections(i);
			lines.push_back(connect_crossings(i));

		}
		for (float i = v_min->s; i < v_max->s; i += n_dist * 2) {
			HSV hsv = { hue_mapping(i, v_min, v_max), 1, 1 };
			colors.push_back(hsv2rgb(hsv));
		}
		int len = colors.size();
		std::vector<_RGB> temp_color = colors;
		std::reverse(temp_color.begin(), temp_color.end());
		colors.insert(colors.end(), temp_color.begin(), temp_color.end());
	}

	glutPostRedisplay();
	break;

	// bicolor height field
	case '8' :
	{
		display_mode = 8;

		Vertex* v_max = poly->max_vertex();
		Vertex* v_min = poly->min_vertex();

		double n = 30;
		double n_dist = (v_max->s - v_min->s) / (double)(n);
		for (float i = v_min->s; i < v_max->s; i += n_dist) {
			classify_and_raise_vertices(i, v_min, v_max);
			calculate_intersections(i);
			lines.push_back(connect_crossings(i));

		}

		for (float i = v_min->s; i < v_max->s; i += n_dist * 2) {
			colors.push_back(get_bicolor_mapping(i, v_min, v_max));
		}
		int len = colors.size();
		std::vector<_RGB> temp_color = colors;
		std::reverse(temp_color.begin(), temp_color.end());
		colors.insert(colors.end(), temp_color.begin(), temp_color.end());
	}
	glutPostRedisplay();
	break;

	// rainbow height field 
	case '9':
	{
		display_mode = 9;

		Vertex* v_max = poly->max_vertex();
		Vertex* v_min = poly->min_vertex();

		double n = 20;
		double n_dist = (v_max->s - v_min->s) / (double)(n);
		for (float i = v_min->s; i < v_max->s; i += n_dist) {
			classify_and_raise_vertices(i, v_min, v_max);
			calculate_intersections(i);
			lines.push_back(connect_crossings(i));

		}
		for (float i = v_min->s; i < v_max->s; i += n_dist * 2) {
			HSV hsv = { hue_mapping(i, v_min, v_max), 1, 1 };
			colors.push_back(hsv2rgb(hsv));
		}
		int len = colors.size();
		std::vector<_RGB> temp_color = colors;
		std::reverse(temp_color.begin(), temp_color.end());
		colors.insert(colors.end(), temp_color.begin(), temp_color.end());
	}

	glutPostRedisplay();
	break;

	// rainbow height field with critical points
	case '0':
	{
		display_mode = 0;

		Vertex* v_max = poly->max_vertex();
		Vertex* v_min = poly->min_vertex();

		double n = 40;
		double n_dist = (v_max->s - v_min->s) / (double)(n);
		find_critical_points(v_min, v_max);
		for (float i = v_min->s; i < v_max->s; i += n_dist) {
			classify_and_raise_vertices(i, v_min, v_max);
			calculate_intersections(i);
			lines.push_back(connect_crossings(i));

		}
		for (float i = v_min->s; i < v_max->s; i += n_dist * 2) {
			HSV hsv = { hue_mapping(i, v_min, v_max), 1, 1 };
			colors.push_back(hsv2rgb(hsv));
		}
		int len = colors.size();
		std::vector<_RGB> temp_color = colors;
		std::reverse(temp_color.begin(), temp_color.end());
		colors.insert(colors.end(), temp_color.begin(), temp_color.end());
	}

	glutPostRedisplay();
	break;



	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_tri = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_tri);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	for (int i = 0; i < poly->ntris; i++)
	{
		Triangle* qtemp = poly->tlist[i];

		glBegin(GL_QUADS);
		for (int j = 0; j < 3; j++)
		{
			Vertex* vtemp = qtemp->verts[j];

			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;
			
			icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
			normalize(dp);
			dp *= dmax;

			double dx = -dp.x;
			double dy = -dp.y;

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}

	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->ntris; i++)
	{
		Triangle* qtemp = poly->tlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 3; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// checkerboard pattern display
	{
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 4: // points and lines drawing example
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			drawPolyLine(lines[k], 1.0, 1.0, 0.0, 0.0);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z);
		}
		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		displayIBFV();
		glutPostRedisplay();
	}
	break;

	case 6:
	case 7:
	case 8:
	case 9:
	case 0:
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = {0.0, 0.0, 0.0, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->ntris; i++) {
			Triangle* temp_q = poly->tlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			double test = colors[k].r;
			drawPolyLine(lines[k], 1.0, colors[k].r, colors[k].g, colors[k].b);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z, 0.15, 1.0, 0.0, 0.0);
		}
		break;
	}
	break;

	default:
	{
		// don't draw anything
	}

	}
}
