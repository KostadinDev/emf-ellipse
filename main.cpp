#include <iostream>
#include<iomanip>
#include <cmath>
#include<limits>


#include<math.h>

using namespace std;

const long double epsilon = std::numeric_limits<long double>::epsilon();

/**
 * Returns the the area of an ellipse.
 * @param a length of vertical half-axis
 * @param b length of horizontal half-axis
 * @return Area of an ellipse
 */
long double area(long double a, long double b) {
    return M_PI * a * b;
}


long double Exact(long double a, long double b, int n) {
    if (a < b) {
        long double temp = a;
        a = b;
        b = temp;
    }
    long double e = sqrt(1 - b * b / (a * a));
    long double s, len;
    long double fact = 1;
    int i;
    len = 1;
    for (i = 1; i < n; i++) {
        fact *= (2 * i * (2 * i - 1));
        fact /= i;
        fact /= i;
        fact /= 4;

        s = fact * fact;

        s *= -pow(e, (2.0 * i));
        s /= (2.0 * i - 1);
        len += s;
    }
    len *= 2 * M_PI * a;
    return (len);
}

long double circumference(long double a, long double b) {
    return M_PI * (a + b) * (3 * ((a - b) * (a - b)) /
                             ((a + b) * (a + b) * (sqrt(-3 * ((a - b) * (a - b)) / ((a + b) * (a + b)) + 4) + 10)) + 1);
}

long double b_f_naive(long double a_0, long double b_0, long double a_f) {
    return sqrt((circumference(a_0, b_0) * circumference(a_0, b_0) - 2 * M_PI * M_PI * a_f * a_f)) / (sqrt(2) * M_PI);
}

long double
dFlux(long double B_0, long double B_f, long double a_0, long double a_f, long double b_0, long double b_f) {
    return B_f * area(a_f, b_f) - B_0 * area(a_0, b_0);
}

long double dT(long double t_0, long double t_f) {
    return t_f - t_0;
}

long double emf(long double dFlux, long double dT) {
    return dFlux / dT;
}

long double binarySearchApprox(long double l, long double r, long double C, long double a_f, int max_iterations) {
    long double left_approximation;
    long double right_approximation;
    long double mid_approximation;
    long double old_m = 0;
    int count_iterations = 0;
    while (1) {
        long double m = (r + l) / 2;

        // Check if x is present at mid
        if (a_f > l)
            left_approximation = circumference(a_f, l);
        else
            left_approximation = circumference(l, a_f);
        if (a_f > l)
            right_approximation = circumference(a_f, r);
        else
            right_approximation = circumference(r, a_f);
        if (a_f > l)
            mid_approximation = circumference(a_f, m);
        else
            mid_approximation = circumference(m, a_f);
        if (abs(mid_approximation - C) <= epsilon) {
            cout << "Error < epsilon" << mid_approximation << " " << C << " " << m << endl;
            return m;
        }
        if (abs(left_approximation - C) <= abs(right_approximation - C))
            r = m;
        else
            l = m;
        count_iterations++;
        if (count_iterations == max_iterations) {
            cout << "Max iterations reached" << endl;
            return m;
        }
        if (old_m == m) {
            cout << "Max change has been reached" << endl;
            return m;
        }
        old_m = m;
        cout << "Iteration " << count_iterations << " m: " << m << endl;
    }

    // if we reach here, then element was
    // not present
    return -1;
}

long double binarySearchExact(long double l, long double C, long double a_f, int max_iterations) {
    long double r = l - -215.0 / 1000.0 * l;
    long double left_approximation;
    long double right_approximation;
    long double mid_approximation;
    long double old_m = 0;
    int count_iterations = 0;
    while (1) {
        long double m = (r + l) / 2;

        // Check if x is present at mid
        if (a_f > l)
            left_approximation = Exact(a_f, l, 10000);
        else
            left_approximation = Exact(l, a_f, 10000);
        if (a_f > l)
            right_approximation = Exact(a_f, r, 10000);
        else
            right_approximation = Exact(r, a_f, 10000);
        if (a_f > l)
            mid_approximation = Exact(a_f, m, 10000);
        else
            mid_approximation = Exact(m, a_f, 10000);
        if (abs(mid_approximation - C) <= epsilon) {
            //cout << "Error < epsilon" << mid_approximation << " " << C << " " << m << endl;
            return m;
        }
        if (abs(left_approximation - C) <= abs(right_approximation - C))
            r = m;
        else
            l = m;
        count_iterations++;
        if (count_iterations == max_iterations) {
            //cout << "Max iterations reached" << endl;
            return m;
        }
        if (old_m == m) {
            //cout << "Max change has been reached" << endl;
            return m;
        }
        old_m = m;
        //cout << "Iteration " << count_iterations << " m: " << m << endl;
    }

    // if we reach here, then element was
    // not present
    return -1;

}

void table() {

    //                            a1,b1,af
    long double combinations[5][3] = {
            {10, 1, 5},
            {9,  2, 5},
            {8,  4, 5},
            {7,  5, 5},
            {6,  6, 5}
    };
    long double B_0 = 1;
    long double B_f = 1;
    long double t_0 = 0;
    long double t_f = 1;
    long double circumference;
    long int n = 100000;

    for (int i = 0; i < 5; i++) {
        long double a_0 = combinations[i][0];
        long double b_0 = combinations[i][1];
        long double a_f = combinations[i][2];
        if (a_0 >= b_0) {
            circumference = Exact(a_0, b_0, n);
        } else {
            circumference = Exact(a_0, b_0, n);
        }
        long double b_f = b_f_naive(a_0, b_0, a_f);
        b_f = binarySearchExact(b_f, circumference, a_f, n);
        long double dF = dFlux(B_0, B_f, a_0, a_f, b_0, b_f);
        long double E = emf(dF, dT(t_0, t_f));
        cout
                << "-----------------------------------------------------------------------------------------------------------------------"
                << i + 1 << endl;
        cout << "B_0: " << B_0 << endl;
        cout << "B_f: " << B_f << endl;
        cout << "dT: " << dT(t_0, t_f) << endl;
        cout << "a_0: " << a_0 << endl;
        cout << "b_0: " << b_0 << endl;
        cout << "a_f: " << a_f << endl;
        cout << "b_f: " << b_f << endl;
        cout << "C:   " << circumference << endl;
        cout << "A_0: " << area(a_0, b_0) << endl;
        cout << "A_f: " << area(a_f, b_f) << endl;
        cout << "EMF: " << E << endl;

    }
}

void test() {
    long double C = Exact(10, 1, 10000000);

    cout << "C: " << C << endl;
    long double b_final = b_f_naive(10, 1, 5);
    cout << "b_final approximation: " << b_final << endl;
    long double max_error = -215.0 / 1000.0 * b_final;
    cout << "b_final - max_error: " << b_final - max_error << endl;
    long double b_final_computed_exact = binarySearchExact(b_final, C, 5, 10000);
    long double b_final_computed_approx = binarySearchApprox(b_final, b_final - max_error, circumference(10, 1), 5,
                                                             10000);
    cout << endl << "b_f approx: " << b_final_computed_exact << endl;
    cout
            << "------------------------------------| EXACT |-------------------------------------------------------------------------"
            << endl;
    cout << endl << "Computed Circumference: "
         << Exact(b_final_computed_exact, 5, 10000000);
    cout << endl << "Actual Circumference:   " << C << endl;
    cout << "------------------------------| Approximation |---------------------------------------" << endl;

    cout << endl << "Computed Circumference: " << circumference(b_final_computed_approx, 5);
    cout << endl << "Actual Circumference:   " << circumference(10, 1) << endl;


}

int main() {
    cout << std::setprecision(64) << fixed;
    table();
    return 0;
}
