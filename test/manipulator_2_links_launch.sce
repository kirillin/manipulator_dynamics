// world values
g = 9.82;

// manipulator values
m_1 = 8;
m_2 = 5;
a_1 = 0.5;
a_2 = 0.4;
delta_1 = 60 * %pi / 180;
delta_2 = 40 * %pi / 180;
// for dynamics
x_c1 = -0.5 * a_1;
y_c1 = 0;
I_zz1 = m_1 * a_1^2 / 3;
x_c2 = -0.5 * a_2;
z_c2 = 0;
I_yy2 = m_2 * a_2^2 / 3;

// motors values
I_a = [0.1, 0; 0, 0.1];

// friction values
f_v = [0.5, 0; 0, 0.5];
f_c = [0.1, 0; 0, 0.1];
f_off = [0.0; 0.0];


// simulation values
T = 0.01;


// task angles values
q1 = [0, 1, 1.5, -0.4, 0.2];
q2 = [0, 1, 1.5, -0.4, 0.2];
goals = [q1; q2];
N = max(size(goals));


// import and model appropriate scheme
path = get_absolute_file_path('manipulator_2_links_launch.sce');
scheme_name = path + 'manipulator_2_links.zcos';
importXcosDiagram(scheme_name);
xcos_simulate(scs_m, 4);


// graphs plotting
scf();
time = angles.time;
subplot(2,2,1);
plot2d(time, angles.values);
subplot(2,2,2);
plot2d(time, speeds.values);
subplot(2,2,3);
plot2d(time, accels.values);
subplot(2,2,4);
plot2d(time, torques.values);


// saving results to text file
exp_data_1 = [angles.values, speeds.values, torques.values, time];
fprintfMat(path + "no_accels.txt", exp_data_1);
exp_data_2 = [angles.values, speeds.values, accels.values, torques.values, time];
fprintfMat(path + "with_accels.txt", exp_data_2);
