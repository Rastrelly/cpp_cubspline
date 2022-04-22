#include <cmath>
#include <iostream>

const int n = 10;
const double a = 0;
const double b = 1;
const double h = 0.1; //h = (b-a)/n
const double param = 0.3;
const double arr_p_inv[] = { 0.0000001, 0.000276, 0.000021, 0.00014, 0.00000001, 0.0001083, 0.00018, 0.9, 0.002425, 0.00013, 0.000000001 };
const double arr_rand[] = { 0.0001221, -0.002012, 0.000231, -0.0001234, -0.0051021, -0.0055103, 0.0001528, -0.00521, -0.0051342, -0.008132, 0.000001 };

using namespace std;

double f(double x)
{
	double res = cos(7 * x)*exp(pow(x, 3)) + pow(x, 2);
	return res;
}

double dfdx(double x)
{
	double res = (-7 * sin(7 * x) + 3 * pow(x, 2)*cos(7 * x))*exp(pow(x, 3)) + 2 * x;
	return res;
}

double d2fdx2(double x)
{
	double sin_x = sin(7 * x);
	double cos_x = cos(7 * x);
	double res = (-49 * cos_x - 42 * pow(x, 2)*sin_x + 6 * x*cos_x + 9 * pow(x, 4)*cos_x)*exp(pow(x, 3)) + 2;
	return res;
}

int def_elem_interval(double x)
{
	for (int i = 1; i <= n; i++)
		if ((x >= a + (i - 1)*h) && (x <= a + i * h))
			return i;
	return -1;
}

double* get_moments_cub()
{
	double* moments = (double*)malloc((n + 1) * sizeof(double));
	//a,b,c,k,g
	double A = h / 6;
	double B = 2 * h / 3;
	double C = h / 6;
	double K = 0;
	double* F_knots = (double*)malloc((n + 1) * sizeof(double));
	for (int i = 0; i < n + 1; i++)
		F_knots[i] = f(a + i * h);
	double* G = (double*)malloc((n + 1) * sizeof(double));
	G[0] = 0;
	for (int i = 1; i < n; i++)
		G[i] = (F_knots[i + 1] - 2 * F_knots[i] + F_knots[i - 1]) / h;
	G[n] = 0;
	//direct
	double* alpha = (double*)malloc((n + 1) * sizeof(double));
	double* beta = (double*)malloc((n + 1) * sizeof(double));
	alpha[0] = 0;
	beta[0] = 0;
	alpha[1] = K;
	beta[1] = G[0];
	for (int i = 1; i <= n - 1; i++) {
		alpha[i + 1] = -C / (A*alpha[i] + B);
		beta[i + 1] = (G[i] - A * beta[i]) / (A*alpha[i] + B);
	}
	moments[n] = (G[n] + K * beta[n]) / (1 - K * alpha[n]);
	for (int i = n - 1; i >= 0; i--)
		moments[i] = alpha[i + 1] * moments[i + 1] + beta[i + 1];
	return moments;
}

double* disc_get_moments_cub()
{
	double* moments = (double*)malloc((n + 1) * sizeof(double));
	//a,b,c,k,g
	double A = h / 6;
	double B = 2 * h / 3;
	double C = h / 6;
	double K = 0;
	double* F_knots = (double*)malloc((n + 1) * sizeof(double));
	for (int i = 0; i < n + 1; i++)
		F_knots[i] = f(a + i * h) + arr_rand[i];
	double* G = (double*)malloc((n + 1) * sizeof(double));
	G[0] = 0;
	for (int i = 1; i < n; i++)
		G[i] = (F_knots[i + 1] - 2 * F_knots[i] + F_knots[i - 1]) / h;
	G[n] = 0;
	//direct
	double* alpha = (double*)malloc((n + 1) * sizeof(double));
	double* beta = (double*)malloc((n + 1) * sizeof(double));
	alpha[0] = 0;
	beta[0] = 0;
	alpha[1] = K;
	beta[1] = G[0];
	for (int i = 1; i <= n - 1; i++) {
		alpha[i + 1] = -C / (A*alpha[i] + B);
		beta[i + 1] = (G[i] - A * beta[i]) / (A*alpha[i] + B);
	}
	moments[n] = (G[n] + K * beta[n]) / (1 - K * alpha[n]);
	for (int i = n - 1; i >= 0; i--)
		moments[i] = alpha[i + 1] * moments[i + 1] + beta[i + 1];
	return moments;
}

double cubic_spl(double x)
{
	double* moments = get_moments_cub();
	int i = def_elem_interval(x);
	double m_1 = moments[i - 1];
	double m_2 = moments[i];
	double x_1 = a + (i - 1)*h;
	double x_2 = a + i * h;
	double f_1 = f(x_1);
	double f_2 = f(x_2);
	double res = m_1 * pow(x_2 - x, 3) / (6 * h) + m_2 * pow(x - x_1, 3) / (6 * h) + ((f_2 - f_1) / h - (m_2 - m_1)*h / 6)*(x - x_2) + f_2 - m_2 * pow(h, 2) / 6;
	return res;
}

double disc_cubic_spl(double x)
{
	double* moments = disc_get_moments_cub();
	int i = def_elem_interval(x);
	double m_1 = moments[i - 1];
	double m_2 = moments[i];
	double x_1 = a + (i - 1)*h;
	double x_2 = a + i * h;
	double f_1 = f(x_1) + arr_rand[i - 1];
	double f_2 = f(x_2) + arr_rand[i];
	double res = m_1 * pow(x_2 - x, 3) / (6 * h) + m_2 * pow(x - x_1, 3) / (6 * h) + ((f_2 - f_1) / h - (m_2 - m_1)*h / 6)*(x - x_2) + f_2 - m_2 * pow(h, 2) / 6;
	return res;
}

void print_cubic_spl()
{
	cout << "        Cubic spline" << endl;
	cout << "  x        f(x)       s(x)       |(f-s)(x)|" << endl;
	double max1 = 0;
	for (int i = 0; i < n; i++) {
		double x = a + i * h + param * h;
		double f_ = f(x);
		double s = cubic_spl(x);
		cout << x << "    " << f_ << "    " << s << "    " << fabs(f_ - s) << endl;
		if (fabs(f_ - s) > max1)
			max1 = fabs(f_ - s);
	}
	cout << "Max deviation equals " << max1 << endl;
	cout << "        Cubic spline in the case when f has a discrepancy" << endl;
	cout << "  x        f(x)       s(x)       |(f-s)(x)|" << endl;
	double max2 = 0;
	for (int i = 0; i < n; i++) {
		double x = a + i * h + param * h;
		double f_ = f(x);
		double s = disc_cubic_spl(x);
		cout << x << "    " << f_ << "    " << s << "    " << fabs(f_ - s) << endl;
		if (fabs(f_ - s) > max2)
			max2 = fabs(f_ - s);
	}
	cout << "Max deviation equals " << max2 << endl;
}

int getmax_abs_i(int k, double** arr)
{
	int max_i = k;
	double max_a = fabs(arr[k][k]);
	for (int i = k; i < n - 1; i++)
		for (int j = k; j < n - 1; j++)
			if (fabs(arr[i][j]) > max_a) {
				max_a = fabs(arr[i][j]);
				max_i = i;
			}
	return max_i;
}

int getmax_abs_j(int k, double** arr)
{
	int max_j = k;
	double max_a = fabs(arr[k][k]);
	for (int i = k; i < n - 1; i++)
		for (int j = k; j < n - 1; j++)
			if (fabs(arr[i][j]) > max_a) {
				max_a = fabs(arr[i][j]);
				max_j = j;
			}
	return max_j;
}

double* Gauss_method(double** arr_a, double* vect)
{
	int* order;
	order = (int*)malloc((n - 1) * sizeof(int));
	for (int i = 0; i < n - 1; i++)
		order[i] = i;
	double** arr_a_k_minus_1;
	arr_a_k_minus_1 = (double**)malloc((n - 1) * sizeof(double*));
	for (int i = 0; i < n - 1; i++)
		arr_a_k_minus_1[i] = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n - 1; i++)
		for (int j = 0; j <= n - 1; j++)
			arr_a_k_minus_1[i][j] = (j == (n - 1)) ? vect[i] : arr_a[i][j];
	for (int k = 0; k < n - 1; k++) {
		int ik = getmax_abs_i(k, arr_a_k_minus_1);
		int jk = getmax_abs_j(k, arr_a_k_minus_1);
		if (ik != k) {
			for (int j = k; j <= n - 1; j++) {
				double temp = arr_a_k_minus_1[k][j];
				arr_a_k_minus_1[k][j] = arr_a_k_minus_1[ik][j];
				arr_a_k_minus_1[ik][j] = temp;
			}
		}
		if (jk != k) {
			for (int i = 0; i < n - 1; i++) {
				double temp = arr_a_k_minus_1[i][k];
				arr_a_k_minus_1[i][k] = arr_a_k_minus_1[i][jk];
				arr_a_k_minus_1[i][jk] = temp;
			}
			int temp = order[k];
			order[k] = order[jk];
			order[jk] = temp;
		}
		double diag = 1 / arr_a_k_minus_1[k][k];
		for (int j = k; j <= n; j++) {
			arr_a_k_minus_1[k][j] *= diag;
		}
		for (int i = 0; i < n - 1; i++)
			if (i != k) {
				double temp = arr_a_k_minus_1[i][k];
				for (int j = k; j <= n - 1; j++)
					arr_a_k_minus_1[i][j] -= temp * arr_a_k_minus_1[k][j];
			}
	}
	double* x;
	x = (double*)malloc((n - 1) * sizeof(double));
	for (int i = 0; i < n - 1; i++) {
		x[order[i]] = arr_a_k_minus_1[i][n - 1];
	}
	return x;
}

double** build_matr()
{
	double** matr;
	matr = (double**)malloc((n - 1) * sizeof(double*));
	for (int i = 0; i < n - 1; i++)
		matr[i] = (double*)malloc((n - 1) * sizeof(double));
	for (int i = 0; i < n - 1; i++)
		for (int j = 0; j < n - 1; j++) {
			matr[i][j] = 0;
			if (i == j)
				matr[i][j] += 2 * h / 3;
			if ((i - j == 1) || (j - i == 1))
				matr[i][j] += h / 6;
			double A = 0, B = 0, C = 0;
			if (i == j - 2)
				A = 1 / h;
			if (i == j - 1) {
				A = -2 / h;
				B = 1 / h;
			}
			if (i == j) {
				A = 1 / h;
				B = -2 / h;
				C = 1 / h;
			}
			if (i == j + 1) {
				B = 1 / h;
				C = -2 / h;
			}
			if (i == j + 2)
				C = 1 / h;
			matr[i][j] += (arr_p_inv[j] * A - 2 * arr_p_inv[j + 1] * B + arr_p_inv[j + 2] * C) / h;
		}
	return matr;
}

double* build_hf()
{
	double* f_ = (double*)malloc((n + 1) * sizeof(double));
	for (int i = 0; i <= n; i++)
		f_[i] = f(a + i * h) + arr_rand[i];
	double* hf = (double*)malloc((n - 1) * sizeof(double));
	for (int i = 0; i < n - 1; i++)
		hf[i] = (f_[i] - 2 * f_[i + 1] + f_[i + 2]) / h;
	return hf;
}

double* get_moments_smooth()
{
	double* moments = (double*)malloc((n - 1) * sizeof(double));
	double** matr = build_matr();
	double* hf = build_hf();
	moments = Gauss_method(matr, hf);
	return moments;
}

double* get_min_G(double* moments)
{
	double* y = (double*)malloc((n + 1) * sizeof(double));
	for (int i = 0; i <= n; i++) {
		y[i] = f(a + i * h) + pow(h, 10) + arr_p_inv[i] * moments[i] / h + arr_rand[i];
		if (i > 0)
			y[i] -= 2 * arr_p_inv[i] * moments[i - 1] / h;
		if (i > 1)
			y[i] += arr_p_inv[i] * moments[i - 2] / h;
	}
	return y;
}

double smoothing_spl(double x, double* moments, double* y)
{
	int i = def_elem_interval(x);
	double m_1 = moments[i - 1];
	double m_2 = moments[i];
	double x_1 = a + (i - 1)*h;
	double x_2 = a + i * h;
	double y_1 = y[i - 1];
	double y_2 = y[i];
	double res = m_1 * pow(x_2 - x, 3) / (6 * h) + m_2 * pow(x - x_1, 3) / (6 * h) + (y_1 - pow(h, 2)*m_1 / 6)*(x_2 - x) / h + (y_2 - pow(h, 2)*m_2 / 6)*(x - x_1) / h;
	return res;
}

double smoothing_spl(double x)
{
	double* moments = get_moments_smooth();
	double* y = get_min_G(moments);
	double res = smoothing_spl(x, moments, y);
	return res;
}

void print_smooth_spl()
{
	cout << "        Smoothing spline" << endl;
	cout << "  x        f(x)       s(x)       |(f-s)(x)|" << endl;
	double max = 0;
	for (int i = 0; i < n; i++) {
		double x = a + i * h + param * h;
		double f_ = f(x);
		double s = smoothing_spl(x);
		cout << x << "    " << f_ << "    " << s << "    " << fabs(f_ - s) << endl;
		if (fabs(f_ - s) > max)
			max = fabs(f_ - s);
	}
	cout << "Max deviation equals " << max << endl;
}

int main()
{
	print_cubic_spl();
	print_smooth_spl();
	return 0;
}