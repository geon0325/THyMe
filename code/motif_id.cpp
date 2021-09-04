#include <iostream>
#include <ctime> 

using namespace std;

int id_to_index_static[128] = {
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	17, 19, 18, 20, 19, 21, 20, 22,
	0, 0, 0, 0, 0, 0, 0, 0,
	17, 18, 19, 20, 19, 20, 21, 22,
	17, 19, 19, 21, 18, 20, 20, 22,
	23, 24, 24, 25, 24, 25, 25, 26,
	30, 27, 27, 1, 27, 1, 1, 2,
	28, 3, 3, 4, 29, 5, 5, 6,
	28, 29, 3, 5, 3, 5, 4, 6,
	7, 9, 8, 10, 9, 11, 10, 12,
	28, 3, 29, 5, 3, 4, 5, 6,
	7, 8, 9, 10, 9, 10, 11, 12,
	7, 9, 9, 11, 8, 10, 10, 12,
	13, 14, 14, 15, 14, 15, 15, 16
};

int id_to_index_temporal[128] = {
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	55, 61, 58, 67, 62, 73, 68, 76,
	0, 0, 0, 0, 0, 0, 0, 0,
	56, 59, 63, 69, 64, 70, 74, 77,
	57, 65, 66, 75, 60, 71, 72, 78,
	79, 80, 81, 83, 82, 84, 85, 86, 
	96, 87, 88, 1, 89, 2, 3, 4,
	90, 5, 6, 11, 93, 14, 15, 20,
	91, 94, 7, 16, 8, 17, 12, 21,
	23, 29, 26, 35, 30, 41, 36, 44,
	92, 9, 95, 18, 10, 13, 19, 22,
	24, 27, 31, 37, 32, 38, 42, 45, 
	25, 33, 34, 43, 28, 39, 40, 46,
	47, 48, 49, 51, 50, 52, 53, 54
};

int get_motif_index_static(int deg_a, int deg_b, int deg_c, int C_ab, int C_bc, int C_ca, int g_abc){
	int a = deg_a - (C_ab + C_ca) + g_abc;
	int b = deg_b - (C_bc + C_ab) + g_abc;
	int c = deg_c - (C_ca + C_bc) + g_abc;
	int d = C_ab - g_abc;
	int e = C_bc - g_abc;
	int f = C_ca - g_abc;
	int g = g_abc;

	int motif_id = ((a > 0) << 0) + ((b > 0) << 1) + ((c > 0) << 2) + ((d > 0) << 3) + ((e > 0) << 4) + ((f > 0) << 5) + ((g > 0) << 6);
	return id_to_index_static[motif_id] - 1;
}

int get_motif_index(int deg_a, int deg_b, int deg_c, int C_ab, int C_bc, int C_ca, int g_abc,
			double time_a, double time_b, double time_c){
	int a = deg_a - (C_ab + C_ca) + g_abc;
	int b = deg_b - (C_bc + C_ab) + g_abc;
	int c = deg_c - (C_ca + C_bc) + g_abc;
	int d = C_ab - g_abc;
	int e = C_bc - g_abc;
	int f = C_ca - g_abc;
	int g = g_abc;

	int temp = time_a;
	double rand_a = (double)rand() / RAND_MAX / 1e7;
	double rand_b = (double)rand() / RAND_MAX / 1e7;
	double rand_c = (double)rand() / RAND_MAX / 1e7;

	time_a += rand_a;
	time_b += rand_b;
	time_c += rand_c;

	int time_id = (((time_a < time_b) or (time_a == time_b and rand_a < rand_b)) << 2) 
		    + (((time_b < time_c) or (time_b == time_c and rand_b < rand_c)) << 1) 
		    + (((time_c < time_a) or (time_c == time_a and rand_c < rand_a))  << 0);
	
	int motif_id = 0;
	switch(time_id){
		case 1:
			// time_c < time_b < time_a
			motif_id = ((c > 0) << 0) + ((b > 0) << 1) + ((a > 0) << 2) + ((e > 0) << 3) + ((d > 0) << 4) + ((f > 0) << 5) + ((g > 0) << 6);
			break;
		case 2:
			// time_b < time_a < time_c
			motif_id = ((b > 0) << 0) + ((a > 0) << 1) + ((c > 0) << 2) + ((d > 0) << 3) + ((f > 0) << 4) + ((e > 0) << 5) + ((g > 0) << 6);
			break;
		case 3:
			// time_b < time_c < time_a
			motif_id = ((b > 0) << 0) + ((c > 0) << 1) + ((a > 0) << 2) + ((e > 0) << 3) + ((f > 0) << 4) + ((d > 0) << 5) + ((g > 0) << 6);
			break;
		case 4:
			// time_a < time_c < time_b
			motif_id = ((a > 0) << 0) + ((c > 0) << 1) + ((b > 0) << 2) + ((f > 0) << 3) + ((e > 0) << 4) + ((d > 0) << 5) + ((g > 0) << 6);
			break;
		case 5:
			// time_c < time_a < time_b
			motif_id = ((c > 0) << 0) + ((a > 0) << 1) + ((b > 0) << 2) + ((f > 0) << 3) + ((d > 0) << 4) + ((e > 0) << 5) + ((g > 0) << 6);
			break;
		case 6:
			// time_a < time_b < time_c
			motif_id = ((a > 0) << 0) + ((b > 0) << 1) + ((c > 0) << 2) + ((d > 0) << 3) + ((e > 0) << 4) + ((f > 0) << 5) + ((g > 0) << 6);
			break;
	}

	return id_to_index_temporal[motif_id] - 1;
}
