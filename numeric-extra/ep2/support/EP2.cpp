#include <iostream>
#include <iomanip>
#include <math.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)<(b)?(b):(a))

class cMatrix {
  private:
	double *A;
	int m, n;
  protected:
  public:
    const int row_size() const;
    const int column_size() const;

	cMatrix(const int m, const int n);
	cMatrix(const int m, const int n, const double A[]);

	~cMatrix();

	cMatrix& set(const double A[]);
	cMatrix& eye();
	cMatrix& zero();
	cMatrix transpose();
    double&  operator[](int i);
    double   operator[](int i) const;
	cMatrix  operator*(const cMatrix& B);
	cMatrix  operator*(const double s);
	cMatrix  operator+(const cMatrix& B);
	cMatrix  operator-(const cMatrix& B);
	cMatrix& operator=(const cMatrix& B);
	void householderDecomposition(cMatrix& Q, cMatrix& R);
	void householderBidiagonalization(cMatrix& Q, cMatrix& R, cMatrix& S);
	void output();
};


const int cMatrix::row_size() const {
    return this->m;
}

const int cMatrix::column_size() const {
    return this->n;
}

void cMatrix::output() {
	for (int j = 0; j < this->m; j++) {
		for (int i = 0; i < this->n; i++) {
			std::cout << std::setw(16) << A[j * this->n + i] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

cMatrix::cMatrix(const int m, const int n) : m(m), n(n) {
	A = new double[m * n];
	this->eye();
}

cMatrix::cMatrix(const int m, const int n, const double A[]) : m(m), n(n) {
	this->A = new double[m * n];
	for (int i = 0; i < m * n; i++) this->A[i] = A[i];
}

cMatrix::~cMatrix() {
	delete [] A;
}

cMatrix& cMatrix::set(const double A[]) {
	for (int i = 0; i < m * n; i++) this->A[i] = A[i];
	return *this;
}

cMatrix& cMatrix::eye() {
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			A[j * this->n + i] = i == j ? 1.0 : 0.0;
	return *this;
}

cMatrix& cMatrix::zero() {
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			A[j * this->n + i] = 0.0;
	return *this;
}

cMatrix cMatrix::transpose() {
	cMatrix R(this->n, this->m);
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			R[i * this->m + j] = A[j * this->n + i];
	return R;
}

double& cMatrix::operator[](int i) {
    return this->A[i];
}

double cMatrix::operator[](int i) const {
    return this->A[i];
}

cMatrix cMatrix::operator*(const cMatrix& B) {
	cMatrix R(this->m, B.column_size());
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < B.column_size(); i++) {
			R[j * B.column_size() + i] = 0.0;
			for (int k = 0; k < this->n; k++)
				R[j * B.column_size() + i] +=
                    A[j * this->n + k] * B[k * B.column_size() + i];
		}
	return R;
}

cMatrix cMatrix::operator*(const double s) {
	cMatrix R(this->m, this->n);
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			R[j * this->n + i] = A[j * this->n + i] * s;
	return R;
}

cMatrix cMatrix::operator+(const cMatrix& B) {
	cMatrix R(this->m, this->n);
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			R.A[j * this->n + i] = A[j * this->n + i] + B.A[j * this->n + i];
	return R;
}

cMatrix cMatrix::operator-(const cMatrix& B) {
	cMatrix R(this->m, this->n);
	for (int j = 0; j < this->m; j++)
		for (int i = 0; i < this->n; i++)
			R[j * this->n + i] = A[j * this->n + i] - B[j * this->n + i];
	return R;
}

cMatrix& cMatrix::operator=(const cMatrix& B) {
	if (this->m * this->n != B.row_size() * B.column_size()) {
		delete [] this->A;
		this->A = new double[B.row_size() * B.column_size()];
	}
	this->m = B.row_size();
	this->n = B.column_size();

	for (int i = 0; i < this->m * this->n; i++) A[i] = B[i];
	return *this;
}

void cMatrix::householderDecomposition(cMatrix& Q, cMatrix& R) {
	double mag, alpha;
	cMatrix u(this->m, 1), v(this->m, 1);
	cMatrix P(this->m, this->m), I(this->m, this->m);

	Q = cMatrix(this->m, this->m);
	R = *this;

	for (int i = 0; i < this->n; i++) {
		u.zero(); v.zero();

		mag = 0.0;
		for (int j = i; j < this->m; j++) {
			u[j] = R[j * this->n + i];
			mag += u[j] * u[j];
		}
		mag = sqrt(mag);

		alpha = u[i] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i; j < this->m; j++) {
			v[j] = j == i ? u[j] + alpha : u[j];
			mag += v[j] * v[j];
		}
		mag = sqrt(mag);

		if (mag < 0.0000000001) continue;

		for (int j = i; j < this->m; j++) v[j] /= mag;

		P = I - (v * v.transpose()) * 2.0;

		R = P * R;
		Q = Q * P;
	}
}


double Inner_Product(double u[], double v[], int n)
{
   double inner_product = 0.0;

   for (n--; n >= 0; n--) inner_product +=  u[n] * v[n];

   return inner_product;
}

void cMatrix::householder_tridiagonalization() {
    cMatrix P(this->n, this->n);
    P.eye();


}

int main(int argc, char* argv[]) {

    double A_i[] = {2, -1 , 1, 3,   -1, 1, 4, 2,   1, 4, 2, -1,   3, 2, -1, 1};
	//double temp[] = {2, 1, 3, 5,    -1, 0, 7, 1,    0, -1, -1, 3,   -3, 7, 4, 3,   1, 6, 4, -3};
	cMatrix A__(4, 4, A_i);
//	double temp[] = {4,3,0,2, 2,1,2,1, 4,4,0,3};
//	cMatrix A__(3,4,temp);
	cMatrix Q__(4,4), R__(4,4), S__(4,4);

	A__.output();

	A__.householder_tridiagonalization();

	A__.output();

	//(Q__*R__*S__).output();

	//(Q__*Q__.transpose()).output();
	//(S__*S__.transpose()).output();

	return 0;
}