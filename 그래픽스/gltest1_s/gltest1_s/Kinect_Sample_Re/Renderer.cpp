#include "Renderer.h"
#include <cmath>
#include <ctime>
#define PI 3.1415926535897

void draw_center(void)
{
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f); 
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.2f, 0.0f, 0.0f);
	glEnd();
	glRasterPos3f(0.2f, 0.0f, 0.0f);
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'x');

	glBegin(GL_LINES);
	glColor3f(0.0f, 1.0f, 0.0f); 
	glVertex3f(0.0f, 0.2f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glEnd();
	glRasterPos3f(0.0f, 0.2f, 0.0f);
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'y');

	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 1.0f); 
	glVertex3f(0.0f, 0.0f, -0.2f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glEnd();
	glRasterPos3f(0.0f, 0.0f, -0.2f);
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'z');
}

void idle() {
	static GLuint previousClock = glutGet(GLUT_ELAPSED_TIME);
	static GLuint currentClock = glutGet(GLUT_ELAPSED_TIME);
	static GLfloat deltaT;

	currentClock = glutGet(GLUT_ELAPSED_TIME);
	deltaT = currentClock - previousClock;
	if (deltaT < 1000.0 / 20.0) { return; }
	else { previousClock = currentClock; }

	//char buff[256];
	//sprintf_s(buff, "Frame Rate = %f", 1000.0 / deltaT);
	//frameRate = buff;

	glutPostRedisplay();
}

void close()
{
	glDeleteTextures(1, &dispBindIndex);
	glutLeaveMainLoop();
	CloseHandle(hMutex);
}

void add_quats(float q1[4], float q2[4], float dest[4])
{
	static int count = 0;
	float t1[4], t2[4], t3[4];
	float tf[4];

	vcopy(q1, t1);
	vscale(t1, q2[3]);

	vcopy(q2, t2);
	vscale(t2, q1[3]);

	vcross(q2, q1, t3);
	vadd(t1, t2, tf);
	vadd(t3, tf, tf);
	tf[3] = q1[3] * q2[3] - vdot(q1, q2);

	dest[0] = tf[0];
	dest[1] = tf[1];
	dest[2] = tf[2];
	dest[3] = tf[3];

	if (++count > RENORMCOUNT) {
		count = 0;
		normalize_quat(dest);
	}
}

void reshape(int width, int height)
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(58, (double)width / height, 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
}

void motion(int x, int y)
{
	GLfloat spin_quat[4];
	float gain;
	gain = 2.0; /* trackball gain */

	if (drag_state == GLUT_DOWN)
	{
		if (button_state == GLUT_LEFT_BUTTON)
		{
			trackball(spin_quat,
				(gain * rot_x - 500) / 500,
				(500 - gain * rot_y) / 500,
				(gain * x - 500) / 500,
				(500 - gain * y) / 500);
			add_quats(spin_quat, quat, quat);
		}
		else if (button_state == GLUT_RIGHT_BUTTON)
		{
			t[0] -= (((float)trans_x - x) / 500);
			t[1] += (((float)trans_y - y) / 500);
		}
		else if (button_state == GLUT_MIDDLE_BUTTON)
			t[2] -= (((float)trans_z - y) / 500 * 4);
		else if (button_state == 3 || button_state == 4) // scroll
		{

		}
		//glutPostRedisplay();
	}

	rot_x = x;
	rot_y = y;

	trans_x = x;
	trans_y = y;
	trans_z = y;
}

void mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		if (button == GLUT_LEFT_BUTTON)
		{
			rot_x = x;
			rot_y = y;

			//t[0] = t[0] + 1;


		}
		else if (button == GLUT_RIGHT_BUTTON)
		{
			trans_x = x;
			trans_y = y;
		}
		else if (button == GLUT_MIDDLE_BUTTON)
		{
			//trcon = trcon + 1;
			trans_z = y;
		}
		else if (button == 3 || button == 4)
		{
			const float sign = (static_cast<float>(button)-3.5f) * 2.0f;
			t[2] -= sign * 500 * 0.00015f;
		}
	}

	drag_state = state;
	button_state = button;
}

void vzero(float* v)
{
	v[0] = 0.0f;
	v[1] = 0.0f;
	v[2] = 0.0f;
}

void vset(float* v, float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

void vsub(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] - src2[0];
	dst[1] = src1[1] - src2[1];
	dst[2] = src1[2] - src2[2];
}

void vcopy(const float *v1, float *v2)
{
	register int i;
	for (i = 0; i < 3; i++)
		v2[i] = v1[i];
}

void vcross(const float *v1, const float *v2, float *cross)
{
	float temp[3];

	temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
	vcopy(temp, cross);
}

float vlength(const float *v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void vscale(float *v, float div)
{
	v[0] *= div;
	v[1] *= div;
	v[2] *= div;
}

void vnormal(float *v)
{
	vscale(v, 1.0f / vlength(v));
}

float vdot(const float *v1, const float *v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void vadd(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
	dst[2] = src1[2] + src2[2];
}

void trackball(float q[4], float p1x, float p1y, float p2x, float p2y)
{
	float a[3]; /* Axis of rotation */
	float phi;  /* how much to rotate about axis */
	float p1[3], p2[3], d[3];
	float t;

	if (p1x == p2x && p1y == p2y) {
		/* Zero rotation */
		vzero(q);
		q[3] = 1.0;
		return;
	}

	/*
	 * First, figure out z-coordinates for projection of P1 and P2 to
	 * deformed sphere
	 */
	vset(p1, p1x, p1y, tb_project_to_sphere(TRACKBALLSIZE, p1x, p1y));
	vset(p2, p2x, p2y, tb_project_to_sphere(TRACKBALLSIZE, p2x, p2y));

	/*
	 *  Now, we want the cross product of P1 and P2
	 */
	vcross(p2, p1, a);

	/*
	 *  Figure out how much to rotate around that axis.
	 */
	vsub(p1, p2, d);
	t = vlength(d) / (2.0f*TRACKBALLSIZE);

	/*
	 * Avoid problems with out-of-control values...
	 */
	if (t > 1.0) t = 1.0;
	if (t < -1.0) t = -1.0;
	phi = 2.0f * asin(t);

	axis_to_quat(a, phi, q);
}

void axis_to_quat(float a[3], float phi, float q[4])
{
	vnormal(a);
	vcopy(a, q);
	vscale(q, sin(phi / 2.0f));
	q[3] = cos(phi / 2.0f);
}

float tb_project_to_sphere(float r, float x, float y)
{
	float d, t, z;

	d = sqrt(x*x + y*y);
	if (d < r * 0.70710678118654752440f) {    /* Inside sphere */
		z = sqrt(r*r - d*d);
	}
	else {           /* On hyperbola */
		t = r / 1.41421356237309504880f;
		z = t*t / d;
	}
	return z;
}

void normalize_quat(float q[4])
{
	int i;
	float mag;

	mag = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	for (i = 0; i < 4; i++) q[i] /= mag;
}

void build_rotmatrix(float m[4][4], float q[4])
{
	m[0][0] = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
	m[0][1] = 2.0f * (q[0] * q[1] - q[2] * q[3]);
	m[0][2] = 2.0f * (q[2] * q[0] + q[1] * q[3]);
	m[0][3] = 0.0f;

	m[1][0] = 2.0f * (q[0] * q[1] + q[2] * q[3]);
	m[1][1] = 1.0f - 2.0f * (q[2] * q[2] + q[0] * q[0]);
	m[1][2] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
	m[1][3] = 0.0f;

	m[2][0] = 2.0f * (q[2] * q[0] - q[1] * q[3]);
	m[2][1] = 2.0f * (q[1] * q[2] + q[0] * q[3]);
	m[2][2] = 1.0f - 2.0f * (q[1] * q[1] + q[0] * q[0]);
	m[2][3] = 0.0f;

	m[3][0] = 0.0f;
	m[3][1] = 0.0f;
	m[3][2] = 0.0f;
	m[3][3] = 1.0f;
}

void InitializeWindow(int argc, char* argv[])
{
	// initialize glut settings
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowSize(1000 / 2, 1000 / 2);

	glutInitWindowPosition(0, 0);

	dispWindowIndex = glutCreateWindow("3D Model");

	trackball(quat, 90.0, 0.0, 0.0, 0.0);

	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutSpecialFunc(special);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutCloseFunc(close);
	//GLuint image = load   ("./my_texture.bmp");
	
	//glBindTexture(1,)

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	// bind textures
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);

	reshape(1000, 1000);

	/*glGenTextures(1, &dispBindIndex);
	glBindTexture(GL_TEXTURE_2D, dispBindIndex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);*/
}


void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, 0.1, 200);		
	glTranslatef(t[0], t[1], t[2] - 1.0f);
	glScalef(1, 1, 1);	
	GLfloat m[4][4],m1[4][4], m2[4][4], m3[4][4];
	build_rotmatrix(m, quat);
	gluLookAt(0, 0.2, 2.0, 0, 0, 0, 0, 1.0, 0);

	GLfloat r, g, b;
	glMultMatrixf(&m[0][0]);

	//floor
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 276, 276, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	for (float fl = -1.0; fl < 1.0; fl = fl + 0.2)
	{
		for (float fr = -1.0; fr < 1.0; fr = fr + 0.2)
		{
			glTexCoord2d(0.0, 0.0);
			glVertex3f(-0.2 + fl, 0.0, -0.2 + fr);
			glTexCoord2d(1.0, 0.0);
			glVertex3f(-0.2 + fl, 0.0, 0.2 + fr);
			glTexCoord2d(1.0, 1.0);
			glVertex3f(0.2 + fl, 0.0, 0.2 + fr);
			glTexCoord2d(0.0, 1.0);
			glVertex3f(0.2 + fl, 0.0, -0.2 + fr);
		}

	}
	glEnd();


	//wall
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 512, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, wall);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_2D, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	for (float fl = -1.0; fl < 1.0; fl = fl + 0.2)
	{
		for (float fr = -1.0; fr < 1.0; fr = fr + 0.2)
		{
			glTexCoord2d(0.0, 0.0);
			glVertex3f(-0.2, 0.0, 0.2);
			glTexCoord2d(1.0, 0.0);
			glVertex3f(-0.2, 0.4, 0.2);
			glTexCoord2d(1.0, 1.0);
			glVertex3f(0.2, 0.4, 0.2);
			glTexCoord2d(0.0, 1.0);
			glVertex3f(0.2, 0.0, 0.2);

			glTexCoord2d(1.0, 1.0);
			glVertex3f(-0.2, 0.4, -0.2);			
			glTexCoord2d(1.0, 0.0);
			glVertex3f(-0.2, 0.4, 0.2);			
			glTexCoord2d(0.0, 0.0);
			glVertex3f(-0.2, 0.0, 0.2);
			glTexCoord2d(0.0, 1.0);
			glVertex3f(-0.2, 0.0, -0.2);

			glTexCoord2d(0.0, 1.0);
			glVertex3f(0.2, 0.0, -0.2);
			glTexCoord2d(0.0, 0.0);
			glVertex3f(0.2, 0.0, 0.2);
			glTexCoord2d(1.0, 0.0);
			glVertex3f(0.2, 0.4, 0.2);
			glTexCoord2d(1.0, 1.0);
			glVertex3f(0.2, 0.4, -0.2);

		}

	}
	glEnd();


	// 침대 texture 설정
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 256,256, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels3);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 침대
	for (int jj = 0; jj <q_index3; jj = jj + 1)
	{
		
		glTexCoord2d(0, 0);
		glVertex3f((cos(PI / 2)*vertex3[quad3[jj].V1].X - sin(PI / 2)*vertex3[quad3[jj].V1].Z) + 0.08, vertex3[quad3[jj].V1].Y + 0.045, (sin(PI / 2)*vertex3[quad3[jj].V1].X + cos(PI / 2)*vertex3[quad3[jj].V1].Z) + 0.13);
		glTexCoord2d(1, 0);
		glVertex3f((cos(PI / 2)*vertex3[quad3[jj].V2].X - sin(PI / 2)*vertex3[quad3[jj].V2].Z) + 0.08, vertex3[quad3[jj].V2].Y + 0.045, (sin(PI / 2)*vertex3[quad3[jj].V2].X + cos(PI / 2)*vertex3[quad3[jj].V2].Z) + 0.13);
		glTexCoord2d(1, 1);
		glVertex3f((cos(PI / 2)*vertex3[quad3[jj].V3].X - sin(PI / 2)*vertex3[quad3[jj].V3].Z) + 0.08, vertex3[quad3[jj].V3].Y + 0.045, (sin(PI / 2)*vertex3[quad3[jj].V3].X + cos(PI / 2)*vertex3[quad3[jj].V3].Z) + 0.13);
		glTexCoord2d(0, 1);
		glVertex3f((cos(PI / 2)*vertex3[quad3[jj].V4].X - sin(PI / 2)*vertex3[quad3[jj].V4].Z) + 0.08, vertex3[quad3[jj].V4].Y + 0.045, (sin(PI / 2)*vertex3[quad3[jj].V4].X + cos(PI / 2)*vertex3[quad3[jj].V4].Z) + 0.13);
		
	}
	glEnd();

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 침대
	for (register int j = 0; j < t_index3; j = j + 1)
	{
		glTexCoord2d(0, 0);
		glVertex3f((cos(PI / 2)*vertex3[triangle3[j].V1].X - sin(PI / 2)*vertex3[triangle3[j].V1].Z) + 0.08, vertex3[triangle3[j].V1].Y + 0.04, (sin(PI / 2)*vertex3[triangle3[j].V1].X + cos(PI / 2)*vertex3[triangle3[j].V1].Z) + 0.13);
		glTexCoord2d(1, 0);
		glVertex3f((cos(PI / 2)*vertex3[triangle3[j].V2].X - sin(PI / 2)*vertex3[triangle3[j].V2].Z) + 0.08, vertex3[triangle3[j].V2].Y + 0.04, (sin(PI / 2)*vertex3[triangle3[j].V2].X + cos(PI / 2)*vertex3[triangle3[j].V2].Z) + 0.13);
		glTexCoord2d(1, 1);
		glVertex3f((cos(PI / 2)*vertex3[triangle3[j].V3].X - sin(PI / 2)*vertex3[triangle3[j].V3].Z) + 0.08, vertex3[triangle3[j].V3].Y + 0.04, (sin(PI / 2)*vertex3[triangle3[j].V3].X + cos(PI / 2)*vertex3[triangle3[j].V3].Z) + 0.13);

	}
	glEnd();

	// 의자 texture 설정

	glTexImage2D(GL_TEXTURE_2D, 0, 3, 571, 572, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels4);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 :의자
	for (int jj = 0; jj <q_index4; jj = jj + 1)
	{
		
		glTexCoord2d(0.0,0.0);
		glVertex3f(vertex4[quad4[jj].V1].X + 0.08, vertex4[quad4[jj].V1].Z , vertex4[quad4[jj].V1].Y-0.1);
		glTexCoord2d(0.0,1.0);
		glVertex3f(vertex4[quad4[jj].V2].X + 0.08, vertex4[quad4[jj].V2].Z , vertex4[quad4[jj].V2].Y - 0.1);
		glTexCoord2d(1.0,1.0);
		glVertex3f(vertex4[quad4[jj].V3].X + 0.08, vertex4[quad4[jj].V3].Z , vertex4[quad4[jj].V3].Y - 0.1);
		glTexCoord2d(1.0,0.0);
		glVertex3f(vertex4[quad4[jj].V4].X + 0.08, vertex4[quad4[jj].V4].Z , vertex4[quad4[jj].V4].Y - 0.1);
		
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 의자
	for (register int j = 0; j < t_index4; j = j + 1)
	{
	 glTexCoord2d(0.0, 0.0);
	 glVertex3f(vertex4[triangle4[j].V1].X + 0.08, vertex4[triangle4[j].V1].Z , vertex4[triangle4[j].V1].Y - 0.1);
	 glTexCoord2d(0.0, 1.0);
	 glVertex3f(vertex4[triangle4[j].V2].X + 0.08, vertex4[triangle4[j].V2].Z , vertex4[triangle4[j].V2].Y - 0.1);
	 glTexCoord2d(1.0, 1.0);
	 glVertex3f(vertex4[triangle4[j].V3].X + 0.08, vertex4[triangle4[j].V3].Z , vertex4[triangle4[j].V3].Y - 0.1);
	
	}
	glEnd();

	//TV texture 설정

	glTexImage2D(GL_TEXTURE_2D, 0, 3, 276, 276, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels5);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : TV
	for (int jj = 0; jj <q_index5; jj = jj + 1)
	{
		
		glTexCoord2d(vertex_color5[quad5[jj].V1 -1].X, vertex_color5[quad5[jj].G1 -1].Y);
		glVertex3f((cos(PI / 2)*vertex5[quad5[jj].V1].X - sin(PI / 2)*(-vertex5[quad5[jj].V1].Z)) - 0.17, vertex5[quad5[jj].V1].Y + 0.0875, (sin(PI / 2)*vertex5[quad5[jj].V1].X + cos(PI / 2) *(-vertex5[quad5[jj].V1].Z)) + 0.1);
		glTexCoord2d(vertex_color5[quad5[jj].V2 - 1].X, vertex_color5[quad5[jj].G2 - 1].Y);
		glVertex3f((cos(PI / 2)*vertex5[quad5[jj].V2].X - sin(PI / 2)*(-vertex5[quad5[jj].V2].Z)) - 0.17, vertex5[quad5[jj].V2].Y + 0.0875, (sin(PI / 2)*vertex5[quad5[jj].V2].X + cos(PI / 2) *(-vertex5[quad5[jj].V2].Z)) + 0.1);
		glTexCoord2d(vertex_color5[quad3[jj].V3 - 1].X, vertex_color5[quad3[jj].G3 - 1].Y);
		glVertex3f((cos(PI / 2)*vertex5[quad5[jj].V3].X - sin(PI / 2)*(-vertex5[quad5[jj].V3].Z)) - 0.17, vertex5[quad5[jj].V3].Y + 0.0875, (sin(PI / 2)*vertex5[quad5[jj].V3].X + cos(PI / 2) *(-vertex5[quad5[jj].V3].Z)) + 0.1);
		glTexCoord2d(vertex_color5[quad3[jj].V4 - 1].X, vertex_color5[quad3[jj].G4 - 1].Y);
		glVertex3f((cos(PI / 2)*vertex5[quad5[jj].V4].X - sin(PI / 2)*(-vertex5[quad5[jj].V4].Z)) - 0.17, vertex5[quad5[jj].V4].Y + 0.0875, (sin(PI / 2)*vertex5[quad5[jj].V4].X + cos(PI / 2) *(-vertex5[quad5[jj].V4].Z)) + 0.1);
		}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : TV
	for (register int j = 0; j < t_index5; j = j + 1)
	{
	 glTexCoord2d(vertex_color5[triangle5[j].V1 - 1].X, vertex_color5[triangle5[j].G1 - 1].Y);
	 glVertex3f(vertex5[triangle5[j].V1].X, vertex5[triangle5[j].V1].Y, -vertex5[triangle5[j].V1].Z);
	 glTexCoord2d(vertex_color5[triangle5[j].V2 - 1].X, vertex_color5[triangle5[j].G2 - 1].Y);
	 glVertex3f(vertex5[triangle5[j].V2].X, vertex5[triangle5[j].V2].Y, -vertex5[triangle5[j].V2].Z);
	 glTexCoord2d(vertex_color5[triangle5[j].V3 - 1].X, vertex_color5[triangle5[j].G3 - 1].Y);
	 glVertex3f(vertex5[triangle5[j].V3].X, vertex5[triangle5[j].V3].Y, -vertex5[triangle5[j].V3].Z);
	 
	}
	glEnd();

	//서랍 texture 설정
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 1600, 1200, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels6);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 서랍
	for (int jj = 0; jj <q_index6; jj = jj + 1)
	{
		
		glTexCoord2d(0, 0);
		glVertex3f(vertex6[quad6[jj].V1].X - 0.165, vertex6[quad6[jj].V1].Y, vertex6[quad6[jj].V1].Z*1.3 + 0.06);
		glTexCoord2d(1, 0);
		glVertex3f(vertex6[quad6[jj].V2].X - 0.165, vertex6[quad6[jj].V2].Y, vertex6[quad6[jj].V2].Z*1.3 + 0.06);
		glTexCoord2d(1, 1);
		glVertex3f(vertex6[quad6[jj].V3].X - 0.165, vertex6[quad6[jj].V3].Y, vertex6[quad6[jj].V3].Z*1.3 + 0.06);
		glTexCoord2d(0, 1);
		glVertex3f(vertex6[quad6[jj].V4].X - 0.165, vertex6[quad6[jj].V4].Y, vertex6[quad6[jj].V4].Z*1.3 + 0.06);
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 서랍
	for (register int j = 0; j < t_index6; j = j + 1)
	{
		glTexCoord2d(0, 0);
		glVertex3f(vertex6[triangle6[j].V1].X, vertex6[triangle6[j].V1].Y, vertex6[triangle6[j].V1].Z);
		glTexCoord2d(1, 0);
		glVertex3f(vertex6[triangle6[j].V2].X, vertex6[triangle6[j].V2].Y, vertex6[triangle6[j].V2].Z);
		glTexCoord2d(1, 1);
		glVertex3f(vertex6[triangle6[j].V3].X, vertex6[triangle6[j].V3].Y, vertex6[triangle6[j].V3].Z);
	}
	glEnd();



	//책상 texture 설정
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 512, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels7);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 책상
	for (int jj = 0; jj <q_index7; jj = jj + 1)
	{
		
		glTexCoord2d(1, 0);
		glVertex3f(vertex7[quad7[jj].V1].X + 0.08, vertex7[quad7[jj].V1].Z, -vertex7[quad7[jj].V1].Y - 0.15);
		glTexCoord2d(1, 1);
		glVertex3f(vertex7[quad7[jj].V2].X + 0.08, vertex7[quad7[jj].V2].Z, -vertex7[quad7[jj].V2].Y - 0.15);
		glTexCoord2d(0, 1);
		glVertex3f(vertex7[quad7[jj].V3].X + 0.08, vertex7[quad7[jj].V3].Z, -vertex7[quad7[jj].V3].Y - 0.15);
		glTexCoord2d(0, 0);
		glVertex3f(vertex7[quad7[jj].V4].X + 0.08, vertex7[quad7[jj].V4].Z, -vertex7[quad7[jj].V4].Y - 0.15);
		
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 책상
	for (register int j = 0; j < t_index7; j = j + 1)
	{
		glNormal3f(vertex_normal5[j].X, vertex_normal5[j].Y, vertex_normal5[j].Z);
		glTexCoord2d(0, 0);
		glVertex3f(vertex7[triangle7[j].V1].X, vertex7[triangle7[j].V1].Z, -vertex7[triangle7[j].V1].Y-0.15);
		glTexCoord2d(0, 1);
		glVertex3f(vertex7[triangle7[j].V2].X, vertex7[triangle7[j].V2].Z, -vertex7[triangle7[j].V2].Y - 0.15);
		glTexCoord2d(0.5,0.5);
		glVertex3f(vertex7[triangle7[j].V3].X, vertex7[triangle7[j].V3].Z, -vertex7[triangle7[j].V3].Y - 0.15);
	}
	glEnd();


	//모니터 texture 설정
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 276, 276, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels5);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 모니터
	for (int jj = 0; jj <q_index5; jj = jj + 1)
	{
		
		glTexCoord2d(0, 0);
		glVertex3f(vertex5[quad5[jj].V1].X / 1.5 + 0.08, vertex5[quad5[jj].V1].Y / 1.5 + 0.135, vertex5[quad5[jj].V1].Z / 1.5 - 0.16);
		glTexCoord2d(1, 0);
		glVertex3f(vertex5[quad5[jj].V2].X / 1.5 + 0.08, vertex5[quad5[jj].V2].Y / 1.5 + 0.135, vertex5[quad5[jj].V2].Z / 1.5 - 0.16);
		glTexCoord2d(1, 1);
		glVertex3f(vertex5[quad5[jj].V3].X / 1.5 + 0.08, vertex5[quad5[jj].V3].Y / 1.5 + 0.135, vertex5[quad5[jj].V3].Z / 1.5 - 0.16);
		glTexCoord2d(0, 1);
		glVertex3f(vertex5[quad5[jj].V4].X / 1.5 + 0.08, vertex5[quad5[jj].V4].Y / 1.5 + 0.135, vertex5[quad5[jj].V4].Z / 1.5 - 0.16);
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 모니터
	for (register int j = 0; j < t_index5; j = j + 1)
	{ 
		glTexCoord2d(0, 0);
		glVertex3f(vertex5[triangle5[j].V1].X / 1.5, vertex5[triangle5[j].V1].Y / 1.5, -vertex5[triangle5[j].V1].Z / 1.5);
		glTexCoord2d(1, 0);
		glVertex3f(vertex5[triangle5[j].V2].X / 1.5, vertex5[triangle5[j].V2].Y / 1.5, -vertex5[triangle5[j].V2].Z / 1.5);
		glTexCoord2d(1, 1);
		glVertex3f(vertex5[triangle5[j].V3].X / 1.5, vertex5[triangle5[j].V3].Y / 1.5, -vertex5[triangle5[j].V3].Z / 1.5);

	}
	glEnd();

	// 벌들을 움직이기 위한 설정
	srand((unsigned int)time(NULL)); 
	trcon = trcon - (rand() % 20 - 10); //벌1에 사용할 랜덤 변수로 초기화

	float cosval = cosf(trcon / 50.0);
	float sinval = sinf(trcon / 50.0);
	
	trcon2 = trcon2 - (rand() % 20 - 10); //벌1에 사용할 랜덤 변수로 초기화

	float cosval2 = cosf(trcon2 / 50.0);
	float sinval2 = sinf(trcon2 / 50.0);

	trcon3 = trcon3 - (rand() % 20 - 10); //벌1에 사용할 랜덤 변수로 초기화

	float cosval3 = cosf(trcon3 / 50.0);
	float sinval3 = sinf(trcon3 / 50.0);

	//벌1 행렬
	m1[0][0] = cosval;
	m1[0][1] = 0 - sinval;
	m1[0][2] = cosval;
	m1[0][3] = 0.0f;

	m1[1][0] = sinval;
	m1[1][1] = cosval;
	m1[1][2] = -sinval;
	m1[1][3] = 0.0f;

	m1[2][0] = -sinval2;
	m1[2][1] = sinval3;
	m1[2][2] = cosval;
	m1[2][3] = 0.0f;

	m1[3][0] = 0.0f;
	m1[3][1] = 0.0f;
	m1[3][2] = 0.0f;
	m1[3][3] = 1.0f;

	
	//벌2 행렬
	m2[0][0] = cosval2;
	m2[0][1] = 0 - sinval2;
	m2[0][2] = cosval2;
	m2[0][3] = 0.0f;

	m2[1][0] = sinval2;
	m2[1][1] = cosval2;
	m2[1][2] = -sinval2;
	m2[1][3] = 0.0f;

	m2[2][0] = -sinval3;
	m2[2][1] = sinval;
	m2[2][2] = cosval2;
	m2[2][3] = 0.0f;

	m2[3][0] = 0.0f;
	m2[3][1] = 0.0f;
	m2[3][2] = 0.0f;
	m2[3][3] = 1.0f;

	//벌3 행렬
	m3[0][0] = cosval3;
	m3[0][1] = 0 - sinval3;
	m3[0][2] = cosval3;
	m3[0][3] = 0.0f;

	m3[1][0] = sinval3;
	m3[1][1] = cosval3;
	m3[1][2] = -sinval3;
	m3[1][3] = 0.0f;

	m3[2][0] = -sinval;
	m3[2][1] = sinval2;
	m3[2][2] = cosval3;
	m3[2][3] = 0.0f;

	m3[3][0] = 0.0f;
	m3[3][1] = 0.0f;
	m3[3][2] = 0.0f;
	m3[3][3] = 1.0f;

	//벌들의 texture 설정
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 276, 276, 0, GL_RGB, GL_UNSIGNED_BYTE, mytexels5);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	float xx1, yy1, zz1;
	float xx2, yy2, zz2;
	float xx3, yy3, zz3;

	//벌1
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 벌
	for (int jj = 0; jj <q_index8; jj = jj + 1)
	{
		
		glTexCoord2d(0, 0);
		xx1 = m1[0][0] * (vertex8[quad8[jj].V1].X) + m1[0][1] * (vertex8[quad8[jj].V1].Y + 0.08) + m1[0][2] * vertex8[quad8[jj].V1].Z;
		yy1 = m1[1][0] * (vertex8[quad8[jj].V1].X) + m1[1][1] * (vertex8[quad8[jj].V1].Y + 0.08) + m1[1][2] * vertex8[quad8[jj].V1].Z;
		zz1 = m1[2][0] * (vertex8[quad8[jj].V1].X) + m1[2][1] * (vertex8[quad8[jj].V1].Y + 0.08) + m1[2][2] * vertex8[quad8[jj].V1].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);

		glTexCoord2d(1, 0);
		xx1 = m1[0][0] * (vertex8[quad8[jj].V2].X) + m1[0][1] * (vertex8[quad8[jj].V2].Y + 0.08) + m1[0][2] * vertex8[quad8[jj].V2].Z;
		yy1 = m1[1][0] * (vertex8[quad8[jj].V2].X) + m1[1][1] * (vertex8[quad8[jj].V2].Y + 0.08) + m1[1][2] * vertex8[quad8[jj].V2].Z;
		zz1 = m1[2][0] * (vertex8[quad8[jj].V2].X) + m1[2][1] * (vertex8[quad8[jj].V2].Y + 0.08) + m1[2][2] * vertex8[quad8[jj].V2].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);

		glTexCoord2d(1, 1);
		xx1 = m1[0][0] * (vertex8[quad8[jj].V3].X) + m1[0][1] * (vertex8[quad8[jj].V3].Y + 0.08) + m1[0][2] * vertex8[quad8[jj].V3].Z;
		yy1 = m1[1][0] * (vertex8[quad8[jj].V3].X) + m1[1][1] * (vertex8[quad8[jj].V3].Y + 0.08) + m1[1][2] * vertex8[quad8[jj].V3].Z;
		zz1 = m1[2][0] * (vertex8[quad8[jj].V3].X) + m1[2][1] * (vertex8[quad8[jj].V3].Y + 0.08) + m1[2][2] * vertex8[quad8[jj].V3].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);
		
		glTexCoord2d(0, 1);
		xx1 = m1[0][0] * (vertex8[quad8[jj].V4].X) + m1[0][1] * (vertex8[quad8[jj].V4].Y + 0.08) + m1[0][2] * vertex8[quad8[jj].V4].Z;
		yy1 = m1[1][0] * (vertex8[quad8[jj].V4].X) + m1[1][1] * (vertex8[quad8[jj].V4].Y + 0.08) + m1[1][2] * vertex8[quad8[jj].V4].Z;
		zz1 = m1[2][0] * (vertex8[quad8[jj].V4].X) + m1[2][1] * (vertex8[quad8[jj].V4].Y + 0.08) + m1[2][2] * vertex8[quad8[jj].V4].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);
		
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 벌
	for (register int j = 0; j < t_index8; j = j + 1)
	{
		glTexCoord2d(0, 0);
		xx1 = m1[0][0] * (vertex8[triangle8[j].V1].X) + m1[0][1] * (vertex8[triangle8[j].V1].Y + 0.08) + m1[0][2] * vertex8[triangle8[j].V1].Z;
		yy1 = m1[1][0] * (vertex8[triangle8[j].V1].X) + m1[1][1] * (vertex8[triangle8[j].V1].Y + 0.08) + m1[1][2] * vertex8[triangle8[j].V1].Z;
		zz1 = m1[2][0] * (vertex8[triangle8[j].V1].X) + m1[2][1] * (vertex8[triangle8[j].V1].Y + 0.08) + m1[2][2] * vertex8[triangle8[j].V1].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);

		glTexCoord2d(1, 0);
		xx1 = m1[0][0] * (vertex8[triangle8[j].V2].X) + m1[0][1] * (vertex8[triangle8[j].V2].Y + 0.08) + m1[0][2] * vertex8[triangle8[j].V2].Z;
		yy1 = m1[1][0] * (vertex8[triangle8[j].V2].X) + m1[1][1] * (vertex8[triangle8[j].V2].Y + 0.08) + m1[1][2] * vertex8[triangle8[j].V2].Z;
		zz1 = m1[2][0] * (vertex8[triangle8[j].V2].X) + m1[2][1] * (vertex8[triangle8[j].V2].Y + 0.08) + m1[2][2] * vertex8[triangle8[j].V2].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);

		glTexCoord2d(1, 1);
		xx1 = m1[0][0] * (vertex8[triangle8[j].V3].X) + m1[0][1] * (vertex8[triangle8[j].V3].Y + 0.08) + m1[0][2] * vertex8[triangle8[j].V3].Z;
		yy1 = m1[1][0] * (vertex8[triangle8[j].V3].X) + m1[1][1] * (vertex8[triangle8[j].V3].Y + 0.08) + m1[1][2] * vertex8[triangle8[j].V3].Z;
		zz1 = m1[2][0] * (vertex8[triangle8[j].V3].X) + m1[2][1] * (vertex8[triangle8[j].V3].Y + 0.08) + m1[2][2] * vertex8[triangle8[j].V3].Z;
		glVertex3f(xx1, yy1 + 0.25, zz1);
		
	}
	glEnd();

	
	//벌2
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 벌2
	for (int jj = 0; jj <q_index8; jj = jj + 1)
	{

		glTexCoord2d(0, 0);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V1].X) + m2[0][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V1].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V1].X) + m2[1][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V1].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V1].X) + m2[2][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V1].Z;
		glVertex3f(xx1 + 0.03, yy1 + 0.13, zz1);

		glTexCoord2d(1, 0);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V2].X) + m2[0][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V2].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V2].X) + m2[1][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V2].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V2].X) + m2[2][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V2].Z;
		glVertex3f(xx1 + 0.03, yy1 + 0.13, zz1);

		glTexCoord2d(1, 1);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V3].X) + m2[0][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V3].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V3].X) + m2[1][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V3].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V3].X) + m2[2][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V3].Z;
		glVertex3f(xx1 + 0.03, yy1 + 0.13, zz1);

		glTexCoord2d(0, 1);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V4].X) + m2[0][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V4].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V4].X) + m2[1][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V4].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V4].X) + m2[2][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V4].Z;
		glVertex3f(xx1 + 0.02, yy1 + 0.13, zz1);
		
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 벌2
	for (register int j = 0; j < t_index8; j = j + 1)
	{
		glTexCoord2d(0, 0);
		xx1 = m2[0][0] * (vertex8[triangle8[j].V1].X + 0.15) + m2[0][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m2[0][2] * vertex8[triangle8[j].V1].Z;
		yy1 = m2[1][0] * (vertex8[triangle8[j].V1].X + 0.15) + m2[1][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m2[1][2] * vertex8[triangle8[j].V1].Z;
		zz1 = m2[2][0] * (vertex8[triangle8[j].V1].X + 0.15) + m2[2][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m2[2][2] * vertex8[triangle8[j].V1].Z;
		glVertex3f(xx1 + 0.02, yy1 + 0.13, zz1);

		glTexCoord2d(1, 0);
		xx1 = m2[0][0] * (vertex8[triangle8[j].V2].X + 0.15) + m2[0][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m2[0][2] * vertex8[triangle8[j].V2].Z;
		yy1 = m2[1][0] * (vertex8[triangle8[j].V2].X + 0.15) + m2[1][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m2[1][2] * vertex8[triangle8[j].V2].Z;
		zz1 = m2[2][0] * (vertex8[triangle8[j].V2].X + 0.15) + m2[2][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m2[2][2] * vertex8[triangle8[j].V2].Z;
		glVertex3f(xx1 + 0.02, yy1 + 0.13, zz1);

		glTexCoord2d(1, 1);
		xx1 = m2[0][0] * (vertex8[triangle8[j].V3].X + 0.15) + m2[0][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m2[0][2] * vertex8[triangle8[j].V3].Z;
		yy1 = m2[1][0] * (vertex8[triangle8[j].V3].X + 0.15) + m2[1][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m2[1][2] * vertex8[triangle8[j].V3].Z;
		zz1 = m2[2][0] * (vertex8[triangle8[j].V3].X + 0.15) + m2[2][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m2[2][2] * vertex8[triangle8[j].V3].Z;
		glVertex3f(xx1 + 0.02, yy1 + 0.13, zz1);
		
	}
	glEnd();


	//벌3
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//QUAD 그리기 : 벌3
	for (int jj = 0; jj <q_index8; jj = jj + 1)
	{

		glTexCoord2d(0, 0);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V1].X) + m2[0][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V1].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V1].X) + m2[1][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V1].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V1].X) + m2[2][1] * (vertex8[quad8[jj].V1].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V1].Z;
		glVertex3f(xx1 - 0.03, yy1 + 0.2, zz1);

		glTexCoord2d(1, 0);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V2].X) + m2[0][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V2].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V2].X) + m2[1][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V2].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V2].X) + m2[2][1] * (vertex8[quad8[jj].V2].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V2].Z;
		glVertex3f(xx1 - 0.03, yy1 + 0.2, zz1);

		glTexCoord2d(1, 1);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V3].X) + m2[0][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V3].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V3].X) + m2[1][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V3].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V3].X) + m2[2][1] * (vertex8[quad8[jj].V3].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V3].Z;
		glVertex3f(xx1 - 0.03, yy1 + 0.2, zz1);

		glTexCoord2d(0, 1);
		xx1 = m2[0][0] * (vertex8[quad8[jj].V4].X) + m2[0][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[0][2] * vertex8[quad8[jj].V4].Z;
		yy1 = m2[1][0] * (vertex8[quad8[jj].V4].X) + m2[1][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[1][2] * vertex8[quad8[jj].V4].Z;
		zz1 = m2[2][0] * (vertex8[quad8[jj].V4].X) + m2[2][1] * (vertex8[quad8[jj].V4].Y + 0.03) + m2[2][2] * vertex8[quad8[jj].V4].Z;
		glVertex3f(xx1 - 0.02, yy1 + 0.2, zz1);
		
	}
	glEnd();


	glEnable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLES);//TRIANGLE 그리기 : 벌3
	for (register int j = 0; j < t_index8; j = j + 1)
	{
		glTexCoord2d(0, 0);
		xx1 = m3[0][0] * (vertex8[triangle8[j].V1].X + 0.15) + m3[0][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m3[0][2] * vertex8[triangle8[j].V1].Z;
		yy1 = m3[1][0] * (vertex8[triangle8[j].V1].X + 0.15) + m3[1][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m3[1][2] * vertex8[triangle8[j].V1].Z;
		zz1 = m3[2][0] * (vertex8[triangle8[j].V1].X + 0.15) + m3[2][1] * (vertex8[triangle8[j].V1].Y + 0.03) + m3[2][2] * vertex8[triangle8[j].V1].Z;
		glVertex3f(xx1 - 0.02, yy1 + 0.2, zz1);

		glTexCoord2d(1, 0);
		xx1 = m3[0][0] * (vertex8[triangle8[j].V2].X + 0.15) + m3[0][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m3[0][2] * vertex8[triangle8[j].V2].Z;
		yy1 = m3[1][0] * (vertex8[triangle8[j].V2].X + 0.15) + m3[1][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m3[1][2] * vertex8[triangle8[j].V2].Z;
		zz1 = m3[2][0] * (vertex8[triangle8[j].V2].X + 0.15) + m3[2][1] * (vertex8[triangle8[j].V2].Y + 0.03) + m3[2][2] * vertex8[triangle8[j].V2].Z;
		glVertex3f(xx1 - 0.02, yy1 + 0.2, zz1);

		glTexCoord2d(1, 1);
		xx1 = m3[0][0] * (vertex8[triangle8[j].V3].X + 0.15) + m3[0][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m3[0][2] * vertex8[triangle8[j].V3].Z;
		yy1 = m3[1][0] * (vertex8[triangle8[j].V3].X + 0.15) + m3[1][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m3[1][2] * vertex8[triangle8[j].V3].Z;
		zz1 = m3[2][0] * (vertex8[triangle8[j].V3].X + 0.15) + m3[2][1] * (vertex8[triangle8[j].V3].Y + 0.03) + m3[2][2] * vertex8[triangle8[j].V3].Z;
		glVertex3f(xx1 - 0.02, yy1 + 0.2, zz1);
		
	}
	glEnd();

	glutSwapBuffers();
}

int main(int argc, char* argv[])
{
	vertex = new Vertex[100000];
	vertex_color = new Vertex[100000];
	//mymesh = new Meshmodel[100000];
	
	int i, j, k = 0;
	FILE* f = fopen("carpet.bmp", "rb");
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
	int width = *(int*)&info[18];
	int height = *(int*)&info[22];
	int size = 3 * width * height;
	unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
	fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
	fclose(f);
	for (i = 0; i < width; i++)
		for (j = 0; j < height; j++)
		{
			mytexels[j][i][0] = data[k * 3 + 2];
			mytexels[j][i][1] = data[k * 3 + 1];
			mytexels[j][i][2] = data[k * 3];
			k++;
		}
	k = 0;
	FILE* f2 = fopen("mymap.bmp", "rb");
	unsigned char info2[54];
	fread(info2, sizeof(unsigned char), 54, f2); // read the 54-byte header
	int width2 = *(int*)&info2[18];
	int height2 = *(int*)&info2[22];
	int size2 = 3 * width2 * height2;
	unsigned char* data2 = new unsigned char[size2]; // allocate 3 bytes per pixel
	fread(data2, sizeof(unsigned char), size2, f2); // read the rest of the data at once
	fclose(f2);
	for (i = 0; i < width2; i++)
		for (j = 0; j < height2; j++)
		{
			wall[j][i][0] = data2[k * 3 + 2];
			wall[j][i][1] = data2[k * 3 + 1];
			wall[j][i][2] = data2[k * 3];
			k++;
		}


	//침대
	vertex3 = new Vertex[100000];
	vertex_color3 = new Vertex[100000];
	vertex_normal3 = new Vertex[100000];
	quad3 = new MMesh[100000];
	triangle3 = new MMesh[100000];

	char dsa5[128];
	char dsa6[128];
	i, j, k = 0;
	f = fopen("Bedside_Table_D_default_1_1.bmp", "rb");
	unsigned char info3[54];
	fread(info3, sizeof(unsigned char), 54, f); // read the 54-byte header
												// extract image height and width from header
	int width3 = *(int*)&info3[18];
	int height3 = *(int*)&info3[22];

	int size3 = 3 * width3 * height3;
	unsigned char* data3 = new unsigned char[size3]; // allocate 3 bytes per pixel
	fread(data3, sizeof(unsigned char), size3, f); // read the rest of the data at once
	fclose(f);
	for (i = 0; i < width3; i++)
		for (j = 0; j < height3; j++)
		{
			for (int t = 0; t <3; t++) {
				mytexels3[i][j][t] = data3[k * 3+ (2-t)];
			}
			k++;
		}

	FILE* fp = fopen("Bed.obj", "r");
	int count = 0;
	int num = 0;
	float x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;

	while (!feof(fp)) {
		fscanf(fp, "%s", dsa5);
		while (strcmp(dsa5, "v") != 0)
		{
			fscanf(fp, "%s", dsa5);
		}
		//v탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex3[0].X = x / 15;
		vertex3[0].Y = y / 15;
		vertex3[0].Z = z / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa6, &x, &y, &z);

			if (strcmp(dsa6, "v") == 0)
			{
				vertex3[j].X = x / 15;
				vertex3[j].Y = y / 15;
				vertex3[j].Z = z / 15;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa5);
		while (strcmp(dsa5, "vn") != 0)
		{
			fscanf(fp, "%s", dsa5);
		}
		// vn 탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex_normal3[0].X = x / 15;
		vertex_normal3[0].Y = y / 15;
		vertex_normal3[0].Z = z / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa6, &x, &y, &z);

			if (strcmp(dsa6, "vn") == 0)
			{
				vertex_normal3[j].X = x / 15;
				vertex_normal3[j].Y = y / 15;
				vertex_normal3[j].Z = z / 15;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa5);
		while (strcmp(dsa5, "vt") != 0)
		{
			fscanf(fp, "%s", dsa5);
		}
		// vt탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex_color3[0].X = x / 15;
		vertex_color3[0].Y = y / 15;
		vertex_color3[0].Z = z / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa6, &x, &y, &z);

			if (strcmp(dsa6, "vt") == 0)
			{
				vertex_color3[j].X = x / 15;
				vertex_color3[j].Y = y / 15;
				vertex_color3[j].Z = z / 15;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa5);
		while (strcmp(dsa5, "f") != 0)
		{
			fscanf(fp, "%s", dsa5);
		}
		// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

		count = fscanf(fp, "%f/%f/%f %f/%f/%f %f/%f/%f", &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);

		fscanf(fp, "%s", dsa5);
		if (dsa5 == "\n") { // 3
			triangle3[t_index3].V1 = x - 1;
			triangle3[t_index3].V2 = x2 - 1;
			triangle3[t_index3].V3 = x3 - 1;
			triangle3[t_index3].T1 = y - 1;
			triangle3[t_index3].T2 = y2 - 1;
			triangle3[t_index3].T3 = y3 - 1;
			triangle3[t_index3].G1 = z - 1;
			triangle3[t_index3].G2 = z2 - 1;
			triangle3[t_index3].G3 = z3 - 1;
			t_index3++;
		}
		else { // 4
			quad3[q_index3].V1 = x - 1;
			quad3[q_index3].V2 = x2 - 1;
			quad3[q_index3].V3 = x3 - 1;
			quad3[q_index3].T1 = y - 1;
			quad3[q_index3].T2 = y2 - 1;
			quad3[q_index3].T3 = y3 - 1;
			quad3[q_index3].G1 = z - 1;
			quad3[q_index3].G2 = z2 - 1;
			quad3[q_index3].G3 = z3 - 1;

			char* ptr = strtok(dsa5, "/");
			quad3[q_index3].V4 = atof(ptr) - 1;
			ptr = strtok(NULL, "/");
			quad3[q_index3].T4 = atof(ptr) - 1;
			ptr = strtok(NULL, "/n");
			quad3[q_index3].G4 = atof(dsa5) - 1;

			q_index3++;
		}

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f/%f/%f %f/%f/%f %f/%f/%f", dsa6, &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);


			if (strcmp(dsa6, "f") == 0)
			{
				fscanf(fp, "%s", dsa6);
				if (strstr(dsa6, "/") == NULL) { // 3
					triangle3[t_index3].V1 = x - 1;
					triangle3[t_index3].V2 = x2 - 1;
					triangle3[t_index3].V3 = x3 - 1;
					triangle3[t_index3].T1 = y - 1;
					triangle3[t_index3].T2 = y2 - 1;
					triangle3[t_index3].T3 = y3 - 1;
					triangle3[t_index3].G1 = z - 1;
					triangle3[t_index3].G2 = z2 - 1;
					triangle3[t_index3].G3 = z3 - 1;

					t_index3++;
					fseek(fp, -1, SEEK_CUR);
				}
				else { // 4
					quad3[q_index3].V1 = x - 1;
					quad3[q_index3].V2 = x2 - 1;
					quad3[q_index3].V3 = x3 - 1;
					quad3[q_index3].T1 = y - 1;
					quad3[q_index3].T2 = y2 - 1;
					quad3[q_index3].T3 = y3 - 1;
					quad3[q_index3].G1 = z - 1;
					quad3[q_index3].G2 = z2 - 1;
					quad3[q_index3].G3 = z3 - 1;

					char* ptr = strtok(dsa6, "/");
					quad3[q_index3].V4 = atof(ptr) - 1;
					ptr = strtok(NULL, "/");
					quad3[q_index3].T4 = atof(ptr) - 1;
					ptr = strtok(NULL, "/n");
					quad3[q_index3].G4 = atof(dsa5) - 1;

					q_index3++;
				}

			}
			else if (strcmp(dsa6, "s") == 0) {
				continue;
			}
			else {
				//break;
			}
		}

		break;
	}

	//의자
	vertex4 = new Vertex[100000];
	vertex_color4 = new Vertex[100000];
	vertex_normal4 = new Vertex[100000];
	quad4 = new MMesh[100000];
	triangle4 = new MMesh[100000];

	char dsa7[128];
	char dsa8[128];
	i, j, k = 0;
	f = fopen("mur.bmp", "rb");
	unsigned char info4[54];
	fread(info4, sizeof(unsigned char), 54, f); // read the 54-byte header
												// extract image height and width from header
	int width4 = *(int*)&info4[18];
	int height4 = *(int*)&info4[22];

	int size4 = 3 * width4 * height4;
	unsigned char* data4 = new unsigned char[size4]; // allocate 3 bytes per pixel
	fread(data4, sizeof(unsigned char), size4, f); // read the rest of the data at once
	fclose(f);
	for (i = 0; i < width4; i++)
		for (j = 0; j < height4; j++)
		{
			for (int t = 0; t <3; t++) {
				mytexels4[i][j][t] = data4[k * 3 + (2-t)];
			}
			k++;
		}



	fp = fopen("Chair.obj", "r");
	count = 0;
	num = 0;
	x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;

	while (!feof(fp)) {
		fscanf(fp, "%s", dsa7);
		while (strcmp(dsa7, "v") != 0)
		{
			fscanf(fp, "%s", dsa7);
		}
		//v탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex4[0].X = x / 800;
		vertex4[0].Y = y / 800;
		vertex4[0].Z = z / 800;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa8, &x, &y, &z);

			if (strcmp(dsa8, "v") == 0)
			{
				vertex4[j].X = x / 800;
				vertex4[j].Y = y / 800;
				vertex4[j].Z = z / 800;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa7);
		while (strcmp(dsa7, "vn") != 0)
		{
			fscanf(fp, "%s", dsa7);
		}
		// vn 탐색


		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex_normal4[0].X = x / 800;
		vertex_normal4[0].Y = y / 800;
		vertex_normal4[0].Z = z / 800;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa8, &x, &y, &z);

			if (strcmp(dsa8, "vn") == 0)
			{
				vertex_normal4[j].X = x / 800;
				vertex_normal4[j].Y = y / 800;
				vertex_normal4[j].Z = z / 800;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa7);
		while (strcmp(dsa7, "vt") != 0)
		{
			fscanf(fp, "%s", dsa7);
		}
		// vt탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex_color4[0].X = x / 800;
		vertex_color4[0].Y = y / 800;
		vertex_color4[0].Z = z / 800;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa8, &x, &y, &z);

			if (strcmp(dsa8, "vt") == 0)
			{
				vertex_color4[j].X = x / 800;
				vertex_color4[j].Y = y / 800;
				vertex_color4[j].Z = z / 800;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa7);
		while (strcmp(dsa7, "f") != 0)
		{
			fscanf(fp, "%s", dsa7);
		}
		// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

		count = fscanf(fp, "%f/%f/%f %f/%f/%f %f/%f/%f", &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);

		fscanf(fp, "%s", dsa7);
		if (dsa7 == "\n") { // 3
			triangle4[t_index4].V1 = x - 1;
			triangle4[t_index4].V2 = x2 - 1;
			triangle4[t_index4].V3 = x3 - 1;
			triangle4[t_index4].T1 = y - 1;
			triangle4[t_index4].T2 = y2 - 1;
			triangle4[t_index4].T3 = y3 - 1;
			triangle4[t_index4].G1 = z - 1;
			triangle4[t_index4].G2 = z2 - 1;
			triangle4[t_index4].G3 = z3 - 1;
			t_index4++;
		}
		else { // 4
			quad4[q_index4].V1 = x - 1;
			quad4[q_index4].V2 = x2 - 1;
			quad4[q_index4].V3 = x3 - 1;
			quad4[q_index4].T1 = y - 1;
			quad4[q_index4].T2 = y2 - 1;
			quad4[q_index4].T3 = y3 - 1;
			quad4[q_index4].G1 = z - 1;
			quad4[q_index4].G2 = z2 - 1;
			quad4[q_index4].G3 = z3 - 1;

			char* ptr = strtok(dsa7, "/");
			quad4[q_index4].V4 = atof(ptr) - 1;
			ptr = strtok(NULL, "/");
			quad4[q_index4].T4 = atof(ptr) - 1;
			ptr = strtok(NULL, "/n");
			quad4[q_index4].G4 = atof(dsa7) - 1;

			q_index4++;
		}

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f/%f/%f %f/%f/%f %f/%f/%f", dsa8, &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);


			if (strcmp(dsa8, "f") == 0)
			{
				fscanf(fp, "%s", dsa8);
				if (strstr(dsa8, "/") == NULL) { // 3
					triangle4[t_index4].V1 = x - 1;
					triangle4[t_index4].V2 = x2 - 1;
					triangle4[t_index4].V3 = x3 - 1;
					triangle4[t_index4].T1 = y - 1;
					triangle4[t_index4].T2 = y2 - 1;
					triangle4[t_index4].T3 = y3 - 1;
					triangle4[t_index4].G1 = z - 1;
					triangle4[t_index4].G2 = z2 - 1;
					triangle4[t_index4].G3 = z3 - 1;

					t_index4++;
					fseek(fp, -1, SEEK_CUR);
				}
				else { // 4
					quad4[q_index4].V1 = x - 1;
					quad4[q_index4].V2 = x2 - 1;
					quad4[q_index4].V3 = x3 - 1;
					quad4[q_index4].T1 = y - 1;
					quad4[q_index4].T2 = y2 - 1;
					quad4[q_index4].T3 = y3 - 1;
					quad4[q_index4].G1 = z - 1;
					quad4[q_index4].G2 = z2 - 1;
					quad4[q_index4].G3 = z3 - 1;

					char* ptr = strtok(dsa8, "/");
					quad4[q_index4].V4 = atof(ptr) - 1;
					ptr = strtok(NULL, "/");
					quad4[q_index4].T4 = atof(ptr) - 1;
					ptr = strtok(NULL, "/n");
					quad4[q_index4].G4 = atof(dsa7) - 1;

					q_index4++;
				}

			}
			else if (strcmp(dsa8, "s") == 0) {
				continue;
			}
			else {
				//break;
			}
		}

		break;
	}


	//TV
	vertex5 = new Vertex[100000];
	vertex_color5 = new Vertex[100000];
	vertex_normal5 = new Vertex[100000];
	quad5 = new MMesh[100000];
	triangle5 = new MMesh[100000];

	char dsa9[300];
	char dsa10[300];
	i, j, k = 0;
	f = fopen("Texture-LCD3.bmp", "rb");
	unsigned char info5[54];
	fread(info5, sizeof(unsigned char), 54, f); // read the 54-byte header
												// extract image height and width from header
	int width5 = *(int*)&info5[18];
	int height5 = *(int*)&info5[22];

	int size5 = 3 * width5 * height5;
	unsigned char* data5 = new unsigned char[size5]; // allocate 3 bytes per pixel
	fread(data5, sizeof(unsigned char), size5, f); // read the rest of the data at once
	fclose(f);
	for (i = 0; i < width5; i++)
		for (j = 0; j < height5; j++)
		{
			for (int t = 0; t <3; t++) {
				mytexels5[i][j][t] = data5[k * 3 + (2-t)];
			}
			k++;
		}



	fp = fopen("LCD TV.obj", "r");
	count = 0;
	num = 0;
	x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;

	while (!feof(fp)) {
		fscanf(fp, "%s", dsa9);
		while (strcmp(dsa9, "v") != 0)
		{
			fscanf(fp, "%s", dsa9);
		}
		//v탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex5[0].X = x / 15;
		vertex5[0].Y = y / 15;
		vertex5[0].Z = z / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa10, &x, &y, &z);

			if (strcmp(dsa10, "v") == 0)
			{
				vertex5[j].X = x / 15;
				vertex5[j].Y = y / 15;
				vertex5[j].Z = z / 15;
			}
			else
				break;
		}
		
		fscanf(fp, "%s", dsa9);
		while (strcmp(dsa9, "vn") != 0)
		{
			fscanf(fp, "%s", dsa9);
		}
		// vn 탐색


		count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
		vertex_normal5[0].X = x / 15;
		vertex_normal5[0].Y = y / 15;
		vertex_normal5[0].Z = z / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa10, &x, &y, &z);

			if (strcmp(dsa10, "vn") == 0)
			{
				vertex_normal5[j].X = x / 15;
				vertex_normal5[j].Y = y / 15;
				vertex_normal5[j].Z = z / 15;
			}
			else
				break;
		}
		
		fscanf(fp, "%s", dsa9);
		while (strcmp(dsa9, "vt") != 0)
		{
			fscanf(fp, "%s", dsa9);
		}
		// vt탐색

		count = fscanf(fp, "%f %f %f /n", &x, &y);
		vertex_color5[0].X = x / 15;
		vertex_color5[0].Y = y / 15;

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f %f %f /n", dsa10, &x, &y);

			if (strcmp(dsa10, "vt") == 0)
			{
				vertex_color5[j].X = x / 15;
				vertex_color5[j].Y = y / 15;
			}
			else
				break;
		}

		fscanf(fp, "%s", dsa9);
		while (strcmp(dsa9, "f") != 0)
		{
			fscanf(fp, "%s", dsa9);
		}
		// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

		count = fscanf(fp, "%f//%f %f//%f %f//%f", &x,  &z, &x2,  &z2, &x3,  &z3);

		fscanf(fp, "%s", dsa9);
		if (dsa9 == "\n") { // 3
			triangle5[t_index5].V1 = x - 1;
			triangle5[t_index5].V2 = x2 - 1;
			triangle5[t_index5].V3 = x3 - 1;
			triangle5[t_index5].G1 = z - 1;
			triangle5[t_index5].G2 = z2 - 1;
			triangle5[t_index5].G3 = z3 - 1;
			t_index5++;
		}
		else { // 4
			quad5[q_index5].V1 = x - 1;
			quad5[q_index5].V2 = x2 - 1;
			quad5[q_index5].V3 = x3 - 1;
			quad5[q_index5].G1 = z - 1;
			quad5[q_index5].G2 = z2 - 1;
			quad5[q_index5].G3 = z3 - 1;

			char* ptr = strtok(dsa9, "/");
			quad5[q_index5].V4 = atof(ptr) - 1;
			ptr = strtok(NULL, "/");
			quad5[q_index5].G4 = atof(dsa9) - 1;

			q_index5++;
		}

		for (int j = 1; j < 100000; j = j + 1)
		{
			count = fscanf(fp, "%s %f//%f %f//%f %f//%f", dsa10, &x,  &z, &x2,  &z2, &x3,  &z3);


			if (strcmp(dsa10, "f") == 0)
			{
				fscanf(fp, "%s", dsa10);
				if (strstr(dsa10, "/") == NULL) { // 3
					triangle5[t_index5].V1 = x - 1;
					triangle5[t_index5].V2 = x2 - 1;
					triangle5[t_index5].V3 = x3 - 1;
					triangle5[t_index5].G1 = z - 1;
					triangle5[t_index5].G2 = z2 - 1;
					triangle5[t_index5].G3 = z3 - 1;

					t_index5++;
					fseek(fp, -1, SEEK_CUR);
				}
				else { // 4
					quad5[q_index5].V1 = x - 1;
					quad5[q_index5].V2 = x2 - 1;
					quad5[q_index5].V3 = x3 - 1;
					quad5[q_index5].G1 = z - 1;
					quad5[q_index5].G2 = z2 - 1;
					quad5[q_index5].G3 = z3 - 1;

					char* ptr = strtok(dsa10, "/");
					quad5[q_index5].V4 = atof(ptr) - 1;
					ptr = strtok(NULL, "/");
					quad5[q_index5].G4 = atof(dsa9) - 1;

					q_index5++;
				}

			}
			else if (strcmp(dsa10, "s") == 0) {
				continue;
			}
			else {
				//break;
			}
		}

		break;
	}

	//서랍
	vertex6 = new Vertex[100000];
	vertex_color6 = new Vertex[100000];
	vertex_normal6 = new Vertex[100000];
	quad6 = new MMesh[100000];
	triangle6 = new MMesh[100000];

	char dsa11[300];
	char dsa12[300];
	i, j, k = 0;
	f = fopen("Metal_Texture_PLUS_Metal_Grid.bmp", "rb");
	unsigned char info6[54];
	fread(info6, sizeof(unsigned char), 54, f); // read the 54-byte header
												// extract image height and width from header
	int width6 = *(int*)&info6[18];
	int height6 = *(int*)&info6[22];

	int size6 = 3 * width6 * height6;
	unsigned char* data6 = new unsigned char[size6]; // allocate 3 bytes per pixel
	fread(data6, sizeof(unsigned char), size6, f); // read the rest of the data at once
	fclose(f);
	for (i = 0; i < width6; i++)
		for (j = 0; j < height6; j++)
		{
			for (int t = 0; t <3; t++) {
				mytexels6[i][j][t] = data6[k * 3 +(2- t)];
			}
			k++;
		}



	fp = fopen("final-tabel2.obj", "r");
	count = 0;
	num = 0;
	x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;
		while (!feof(fp)) {
			fscanf(fp, "%s", dsa11);
			while (strcmp(dsa11, "v") != 0)
			{
				fscanf(fp, "%s", dsa11);
			}
			//v탐색

			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex6[0].X = x / 150;
			vertex6[0].Y = y / 150;
			vertex6[0].Z = z / 150;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa12, &x, &y, &z);

				if (strcmp(dsa12, "v") == 0)
				{
					vertex6[j].X = x / 150;
					vertex6[j].Y = y / 150;
					vertex6[j].Z = z / 150;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa11);
			while (strcmp(dsa11, "vn") != 0)
			{
				fscanf(fp, "%s", dsa11);
			}
			// vn 탐색


			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex_normal6[0].X = x / 150;
			vertex_normal6[0].Y = y / 150;
			vertex_normal6[0].Z = z / 150;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa12, &x, &y, &z);

				if (strcmp(dsa12, "vn") == 0)
				{
					vertex_normal6[j].X = x / 150;
					vertex_normal6[j].Y = y / 150;
					vertex_normal6[j].Z = z / 150;
				}
				else
					break;
			}
			
			fscanf(fp, "%s", dsa11);
			while (strcmp(dsa11, "f") != 0)
			{
				fscanf(fp, "%s", dsa11);
			}
			// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

			count = fscanf(fp, "%f//%f %f//%f %f//%f", &x, &z, &x2, &z2, &x3, &z3);

			fscanf(fp, "%s", dsa11);
			if (dsa11 == "\n") { // 3
				triangle6[t_index6].V1 = x - 1;
				triangle6[t_index6].V2 = x2 - 1;
				triangle6[t_index6].V3 = x3 - 1;
				triangle6[t_index6].G1 = z - 1;
				triangle6[t_index6].G2 = z2 - 1;
				triangle6[t_index6].G3 = z3 - 1;
				t_index6++;
			}
			else { // 4
				quad6[q_index6].V1 = x - 1;
				quad6[q_index6].V2 = x2 - 1;
				quad6[q_index6].V3 = x3 - 1;
				quad6[q_index6].G1 = z - 1;
				quad6[q_index6].G2 = z2 - 1;
				quad6[q_index6].G3 = z3 - 1;

				char* ptr = strtok(dsa11, "/");
				quad6[q_index6].V4 = atof(ptr) - 1;
				ptr = strtok(NULL, "/");
				quad6[q_index6].G4 = atof(dsa11) - 1;

				q_index6++;
			}

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f//%f %f//%f %f//%f", dsa12, &x, &z, &x2, &z2, &x3, &z3);


				if (strcmp(dsa12, "f") == 0)
				{
					fscanf(fp, "%s", dsa12);
					if (strstr(dsa12, "/") == NULL) { // 3
						triangle6[t_index6].V1 = x - 1;
						triangle6[t_index6].V2 = x2 - 1;
						triangle6[t_index6].V3 = x3 - 1;
						triangle6[t_index6].G1 = z - 1;
						triangle6[t_index6].G2 = z2 - 1;
						triangle6[t_index6].G3 = z3 - 1;

						t_index6++;
						fseek(fp, -1, SEEK_CUR);
					}
					else { // 4
						quad6[q_index6].V1 = x - 1;
						quad6[q_index6].V2 = x2 - 1;
						quad6[q_index6].V3 = x3 - 1;
						quad6[q_index6].G1 = z - 1;
						quad6[q_index6].G2 = z2 - 1;
						quad6[q_index6].G3 = z3 - 1;

						char* ptr = strtok(dsa12, "/");
						quad6[q_index6].V4 = atof(ptr) - 1;
						ptr = strtok(NULL, "/");
						quad6[q_index6].G4 = atof(dsa11) - 1;

						q_index6++;
					}

				}
				else if (strcmp(dsa12, "s") == 0) {
					continue;
				}
				else {
					//break;
				}
			}

			break;
		}
		
		//책상
		vertex7 = new Vertex[100000];
		vertex_color7 = new Vertex[100000];
		vertex_normal7 = new Vertex[100000];
		quad7 = new MMesh[100000];
		triangle7 = new MMesh[100000];


		char dsa13[128];
		char dsa14[128];
		i, j, k = 0;
		f = fopen("sol.bmp", "rb");
		unsigned char info7[54];
		fread(info7, sizeof(unsigned char), 54, f); // read the 54-byte header
													// extract image height and width from header
		int width7 = *(int*)&info7[18];
		int height7 = *(int*)&info7[22];

		int size7 = 3 * width7 * height7;
		unsigned char* data7 = new unsigned char[size7]; // allocate 3 bytes per pixel
		fread(data7, sizeof(unsigned char), size7, f); // read the rest of the data at once
		fclose(f);
		for (i = 0; i < width7; i++)
			for (j = 0; j < height7; j++)
			{
				mytexels7[i][j][0] = data7[k * 3 + 2];
				mytexels7[i][j][1] = data7[k * 3 + 1];
				mytexels7[i][j][2] = data7[k * 3];
				k++;
			}



		fp = fopen("table.obj", "r");
		count = 0;
		num = 0;
		x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;

		while (!feof(fp)) {
			fscanf(fp, "%s", dsa13);
			while (strcmp(dsa13, "v") != 0)
			{
				fscanf(fp, "%s", dsa13);
			}
			//v탐색

			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex7[0].X = x / 650;
			vertex7[0].Y = y / 650;
			vertex7[0].Z = z / 650;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa14, &x, &y, &z);

				if (strcmp(dsa14, "v") == 0)
				{
					vertex7[j].X = x / 650;
					vertex7[j].Y = y / 650;
					vertex7[j].Z = z / 650;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa13);
			while (strcmp(dsa13, "vn") != 0)
			{
				fscanf(fp, "%s", dsa13);
			}
			// vn 탐색


			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex_normal7[0].X = x / 650;
			vertex_normal7[0].Y = y / 650;
			vertex_normal7[0].Z = z / 650;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa14, &x, &y, &z);

				if (strcmp(dsa14, "vn") == 0)
				{
					vertex_normal7[j].X = x / 650;
					vertex_normal7[j].Y = y / 650;
					vertex_normal7[j].Z = z / 650;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa13);
			while (strcmp(dsa13, "vt") != 0)
			{
				fscanf(fp, "%s", dsa13);
			}
			// vt탐색

			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex_color7[0].X = x / 650;
			vertex_color7[0].Y = y / 650;
			vertex_color7[0].Z = z / 650;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa14, &x, &y, &z);

				if (strcmp(dsa14, "vt") == 0)
				{
					vertex_color7[j].X = x / 650;
					vertex_color7[j].Y = y / 650;
					vertex_color7[j].Z = z / 650;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa13);
			while (strcmp(dsa13, "f") != 0)
			{
				fscanf(fp, "%s", dsa13);
			}
			// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

			count = fscanf(fp, "%f/%f/%f %f/%f/%f %f/%f/%f", &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);

			fscanf(fp, "%s", dsa13);
			if (dsa13 == "\n") { // 3
				triangle7[t_index7].V1 = x - 1;
				triangle7[t_index7].V2 = x2 - 1;
				triangle7[t_index7].V3 = x3 - 1;
				triangle7[t_index7].T1 = y - 1;
				triangle7[t_index7].T2 = y2 - 1;
				triangle7[t_index7].T3 = y3 - 1;
				triangle7[t_index7].G1 = z - 1;
				triangle7[t_index7].G2 = z2 - 1;
				triangle7[t_index7].G3 = z3 - 1;
				t_index7++;
			}
			else { // 4
				quad7[q_index7].V1 = x - 1;
				quad7[q_index7].V2 = x2 - 1;
				quad7[q_index7].V3 = x3 - 1;
				quad7[q_index7].T1 = y - 1;
				quad7[q_index7].T2 = y2 - 1;
				quad7[q_index7].T3 = y3 - 1;
				quad7[q_index7].G1 = z - 1;
				quad7[q_index7].G2 = z2 - 1;
				quad7[q_index7].G3 = z3 - 1;

				char* ptr = strtok(dsa13, "/");
				quad7[q_index7].V4 = atof(ptr) - 1;
				ptr = strtok(NULL, "/");
				quad7[q_index7].T4 = atof(ptr) - 1;
				ptr = strtok(NULL, "/n");
				quad7[q_index7].G4 = atof(dsa13) - 1;

				q_index7++;
			}

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f/%f/%f %f/%f/%f %f/%f/%f", dsa14, &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);


				if (strcmp(dsa14, "f") == 0)
				{
					fscanf(fp, "%s", dsa14);
					if (strstr(dsa14, "/") == NULL) { // 3
						triangle7[t_index7].V1 = x - 1;
						triangle7[t_index7].V2 = x2 - 1;
						triangle7[t_index7].V3 = x3 - 1;
						triangle7[t_index7].T1 = y - 1;
						triangle7[t_index7].T2 = y2 - 1;
						triangle7[t_index7].T3 = y3 - 1;
						triangle7[t_index7].G1 = z - 1;
						triangle7[t_index7].G2 = z2 - 1;
						triangle7[t_index7].G3 = z3 - 1;

						t_index7++;
						fseek(fp, -1, SEEK_CUR);
					}
					else { // 4
						quad7[q_index7].V1 = x - 1;
						quad7[q_index7].V2 = x2 - 1;
						quad7[q_index7].V3 = x3 - 1;
						quad7[q_index7].T1 = y - 1;
						quad7[q_index7].T2 = y2 - 1;
						quad7[q_index7].T3 = y3 - 1;
						quad7[q_index7].G1 = z - 1;
						quad7[q_index7].G2 = z2 - 1;
						quad7[q_index7].G3 = z3 - 1;

						char* ptr = strtok(dsa14, "/");
						quad7[q_index7].V4 = atof(ptr) - 1;
						ptr = strtok(NULL, "/");
						quad7[q_index7].T4 = atof(ptr) - 1;
						ptr = strtok(NULL, "/n");
						quad7[q_index7].G4 = atof(dsa13) - 1;

						q_index7++;
					}

				}
				else if (strcmp(dsa14, "s") == 0) {
					continue;
				}
				else {
					//break;
				}
			}

			break;
		}

		//벌
		vertex8 = new Vertex[100000];
		vertex_color8 = new Vertex[100000];
		vertex_normal8 = new Vertex[100000];
		quad8 = new MMesh[100000];
		triangle8 = new MMesh[100000];


		char dsa15[128];
		char dsa16[128];
		i, j, k = 0;
		f = fopen("Metal_Texture_PLUS_Metal_Grid.bmp", "rb");
		unsigned char info8[54];
		fread(info8, sizeof(unsigned char), 54, f); // read the 54-byte header
													// extract image height and width from header
		int width8 = *(int*)&info8[18];
		int height8 = *(int*)&info8[22];

		int size8 = 30 * width8 * height8;
		unsigned char* data8 = new unsigned char[size8]; // allocate 3 bytes per pixel
		fread(data8, sizeof(unsigned char), size8, f); // read the rest of the data at once
		fclose(f);
		for (i = 0; i < width8; i++)
			for (j = 0; j < height8; j++)
			{
				mytexels8[i][j][0] = data8[k * 30 + 2];
				mytexels8[i][j][1] = data8[k * 30 + 1];
				mytexels8[i][j][2] = data8[k * 30];
				k++;
			}



		fp = fopen("BEE.obj", "r");
		count = 0;
		num = 0;
		x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4 = 0;

		while (!feof(fp)) {
			fscanf(fp, "%s", dsa15);
			while (strcmp(dsa15, "v") != 0)
			{
				fscanf(fp, "%s", dsa15);
			}
			//v탐색

			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex8[0].X = x *1.2;
			vertex8[0].Y = y *1.2;
			vertex8[0].Z = z *1.2;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa16, &x, &y, &z);

				if (strcmp(dsa16, "v") == 0)
				{
					vertex8[j].X = x *1.2;
					vertex8[j].Y = y *1.2;
					vertex8[j].Z = z *1.2;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa15);
			while (strcmp(dsa15, "vt") != 0)
			{
				fscanf(fp, "%s", dsa15);
			}
			// vt탐색

			count = fscanf(fp, "%f %f /n", &x, &y);
			vertex_color8[0].X = x *1.2;
			vertex_color8[0].Y = y *1.2;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f /n", dsa16, &x, &y);

				if (strcmp(dsa16, "vt") == 0)
				{
					vertex_color8[j].X = x *1.2;
					vertex_color8[j].Y = y *1.2;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa15);
			while (strcmp(dsa15, "vn") != 0)
			{
				fscanf(fp, "%s", dsa15);
			}
			// vn 탐색


			count = fscanf(fp, "%f %f %f /n", &x, &y, &z);
			vertex_normal8[0].X = x *1.2;
			vertex_normal8[0].Y = y *1.2;
			vertex_normal8[0].Z = z *1.2;

			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f %f %f /n", dsa16, &x, &y, &z);

				if (strcmp(dsa16, "vn") == 0)
				{
					vertex_normal8[j].X = x *1.2;
					vertex_normal8[j].Y = y *1.2;
					vertex_normal8[j].Z = z *1.2;
				}
				else
					break;
			}

			fscanf(fp, "%s", dsa15);
			while (strcmp(dsa15, "f") != 0)
			{
				fscanf(fp, "%s", dsa15);
			}
			// f탐색 : f에 대한 값이 3개인지 4개인지 구분해야됨

			count = fscanf(fp, "%f/%f/%f %f/%f/%f %f/%f/%f", &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);

			fscanf(fp, "%s", dsa15);
			if (dsa15 == "\n") { // 3
				triangle8[t_index8].V1 = x - 1;
				triangle8[t_index8].V2 = x2 - 1;
				triangle8[t_index8].V3 = x3 - 1;
				triangle8[t_index8].T1 = y - 1;
				triangle8[t_index8].T2 = y2 - 1;
				triangle8[t_index8].T3 = y3 - 1;
				triangle8[t_index8].G1 = z - 1;
				triangle8[t_index8].G2 = z2 - 1;
				triangle8[t_index8].G3 = z3 - 1;
				t_index8++;
			}
			
			for (int j = 1; j < 100000; j = j + 1)
			{
				count = fscanf(fp, "%s %f/%f/%f %f/%f/%f %f/%f/%f", dsa16, &x, &y, &z, &x2, &y2, &z2, &x3, &y3, &z3);


				if (strcmp(dsa16, "f") == 0)
				{
					fscanf(fp, "%s", dsa16);
					if (strstr(dsa16, "/") == NULL) { // 3
						triangle8[t_index8].V1 = x - 1;
						triangle8[t_index8].V2 = x2 - 1;
						triangle8[t_index8].V3 = x3 - 1;
						triangle8[t_index8].T1 = y - 1;
						triangle8[t_index8].T2 = y2 - 1;
						triangle8[t_index8].T3 = y3 - 1;
						triangle8[t_index8].G1 = z - 1;
						triangle8[t_index8].G2 = z2 - 1;
						triangle8[t_index8].G3 = z3 - 1;

						t_index8++;
						fseek(fp, -1, SEEK_CUR);
					}
				}
				else if (strcmp(dsa16, "s") == 0) {
					continue;
				}
				else if (strcmp(dsa16, "s") == 0 && x == 785) {
					break;
				}


			}

			break;
		}

	InitializeWindow(argc, argv);

	display();

	glutMainLoop();
	delete[] vertex;
	delete[] vertex_color;
	
	delete[] vertex8;
	delete[] vertex_color8;

	delete[] vertex3;
	delete[] vertex_color3;

	delete[] vertex4;
	delete[] vertex_color4;

	delete[] vertex5;
	delete[] vertex_color5;

	delete[] vertex6;
	delete[] vertex_color6;

	delete[] vertex7;
	delete[] vertex_color7;

	
	return 0;
}